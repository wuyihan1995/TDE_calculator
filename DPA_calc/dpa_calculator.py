#!/usr/bin/env python3
"""
DPA calculator — analytical tool for alloy materials

Theory references:
  NRT model:   Norgett, Robinson & Torrens (1975). Nucl. Eng. Des., 33, 50-54.
  ARC-DPA:     Nordlund, Zinkle, Sand et al. (2018). Nat. Commun., 9, 1084.
  CB-DPA:      Chen, Bernard, Tommasi & De Saint Jean (2020). EPJ Web Conf., 239, 08003.
  Robinson partition: Robinson, M.T. (1970). ORNL-4552.
  Lindhard theory: Lindhard, Scharff & Schiott (1963). Kgl. Dan. Vid. Selsk. Mat.-Fys. Medd., 33, 1-42.

Input parameters and their physical meanings:
  E_PKA  – Primary Knock-on Atom kinetic energy (eV).
           Energy initially transferred from the incident particle (neutron/ion)
           to a lattice atom, setting it in motion through the crystal.
  Ed     – Threshold displacement energy (eV).
           Minimum kinetic energy a lattice atom must receive to be permanently
           displaced from its lattice site, creating a stable Frenkel pair (1 vacancy
           + 1 interstitial) that survives at least ~1 ps.
  N_atoms – Number of atoms in the simulation cell.
            Used to convert the total displacement count into the dimensionless
            "displacements per atom" (DPA) metric.
  model  – DPA theoretical model:
           "NRT" → Norgett-Robinson-Torrens (1975), the international standard.
                   Displacement efficiency g = 0.8 from binary collision simulations.
                   Overestimates DPA by ~3x (ignores in-cascade recombination).
           "ARC" → Nordlund et al. (2018), athermal recombination corrected.
                   Uses MD-fitted power-law efficiency: ξ=(1-c)·(0.8Ea/2Ed)^b+c.
                   Requires --b_arc and --c_arc (material-specific fitting parameters).
           "CB"  → Chen-Bernard (2020), an improved model without fitting parameters.
                   Includes an athermal-recombination efficiency ξ(Ea) derived from
                   the defect annihilation equation dn/dt = -kn^2.
  partition – How to compute damage energy Ea from PKA energy E_PKA:
           "direct"   → Ea = E_PKA (no correction; assumes electronic stopping
                         loss is negligible, valid for low-energy PKAs ~1-5 keV).
           "robinson" → Robinson analytic fit to the Lindhard energy partition
                         theory. Accounts for electronic stopping losses via
                         the reduced energy ε and the function g(ε).
                         Ea = E_PKA / (1 + k * g(ε)).
  target – Desired DPA dose level.
           Used to compute how many cascades are needed to reach this dose.
           Example: target=0.5 means "how many cascades for 0.5 displacements per atom".
  b_arc  – (ARC-DPA only) Power-law exponent from MD simulation fitting.
           Controls the decay rate of displacement efficiency with energy.
           Typical values: Fe=-0.568, Ni=-1.01, Cu=-0.68, Ag=-1.06.
  c_arc  – (ARC-DPA only) Asymptotic efficiency at high energies from MD fitting.
           Typical values: Fe=0.286, Ni=0.23, Cu=0.16, Ag=0.257.
  name   – Material name for display only (does not affect calculations).
  Z      – Effective atomic number Z_eff = Σ(c_i·Z_i). Used in Robinson partition
           and CB-DPA model. Default: 49.6 (HfNbZrTiTa).
  A      – Effective atomic mass A_eff in g/mol = Σ(c_i·A_i). Used in Robinson
           partition and CB-DPA model. Default: 118.287 (HfNbZrTiTa).
"""

import argparse
import sys

# ── HfNbZrTiTa material constants ──────────────────────────────────────────
MATERIAL = {
    "name":   "HfNbZrTiTa",
    "Z_eff":  49.6,       # effective atomic number  = Σ(c_i·Z_i)
    "A_eff":  118.287,    # effective atomic mass (g/mol) = Σ(c_i·A_i)
    # Individual elements: Hf(72,178.49) Nb(41,92.91) Zr(40,91.22) Ti(22,47.87) Ta(73,180.95)
    # Z_eff = (72+41+40+22+73)/5 = 49.6
    # A_eff = (178.49+92.91+91.22+47.87+180.95)/5 = 118.287
}


# ═══════════════════════════════════════════════════════════════════════════
#  Partition functions: PKA energy → damage energy
# ═══════════════════════════════════════════════════════════════════════════

def partition_direct(epka):
    """Direct mapping: Ea = E_PKA.

    Physical meaning: assumes that ALL of the PKA's kinetic energy goes
    into atomic collisions (no electronic stopping losses). Valid when
    E_PKA is low (<~5 keV) so that electronic stopping is negligible.
    """
    return epka, {"method": "direct", "f": 1.0}


def partition_robinson(epka, z_eff, a_eff):
    """Robinson analytic fit to the Lindhard energy partition theory.

    Physical meaning: a fast-moving ion loses energy through two channels —
    (1) elastic nuclear collisions with lattice atoms (damage-producing), and
    (2) inelastic electronic excitation/ionization (non-damage-producing).
    The Robinson fit estimates what FRACTION of E_PKA survives as damage
    energy Ea available for creating atomic displacements.

    Parameters
    ----------
    epka  : float  – PKA kinetic energy (eV)
    z_eff : float  – effective atomic number of the target material
    a_eff : float  – effective atomic mass of the target material (g/mol)

    Returns
    -------
    Ea    : float  – damage energy (eV), the fraction of E_PKA deposited into
                     atomic collisions after subtracting electronic losses
    info  : dict   – intermediate quantities for diagnostic output

    Reference: Robinson, M.T. (1970). ORNL-4552.

    Key intermediate quantities
    ---------------------------
    ε (reduced energy) – dimensionless energy scale that characterises
                         whether nuclear or electronic stopping dominates.
                         ε = E_PKA / E_Lindhard, where E_Lindhard ~ Z^(7/3).
                         Small ε → nuclear stopping dominates.
    k – dimensionless parameter that scales electronic stopping strength;
        k ∝ Z^(2/3)·A^(-1/2).
    g(ε) – Robinson's universal correction function:
           g(ε) = ε + 0.40244·ε^(3/4) + 3.4008·ε^(1/6)
    Ea = E_PKA / (1 + k·g(ε))
    """
    if epka <= 0:
        return 0.0, {"method": "robinson", "error": "E_PKA <= 0"}

    # Lindhard characteristic energy (eV) for self-ion irradiation
    e_lindhard = 30.74 * (z_eff ** (7.0 / 3.0))

    # Reduced energy ε (dimensionless)
    epsilon = epka / e_lindhard

    # k parameter
    k_val = 0.1337 * (z_eff ** (2.0 / 3.0)) * (a_eff ** (-0.5))

    # Robinson's universal correction function g(ε)
    g_eps = (epsilon
             + 0.40244 * (epsilon ** 0.75)
             + 3.4008 * (epsilon ** (1.0 / 6.0)))

    # Damage energy (eV)
    ea = epka / (1.0 + k_val * g_eps)
    f = ea / epka  # partition fraction

    return ea, {
        "method":  "robinson",
        "f":       f,
        "epsilon": epsilon,
        "k":       k_val,
        "g_eps":   g_eps,
        "E_L":     e_lindhard,
    }


# ═══════════════════════════════════════════════════════════════════════════
#  DPA models: damage energy → number of displaced atoms per cascade
# ═══════════════════════════════════════════════════════════════════════════

def ed_high_regime(ea, ed_val):
    """Check whether Ea falls in the high-E regime (>2Ed/0.8).

    Returns True if the piecewise function uses Nd = 0.8·Ea/(2Ed)·ξ(Ea).
    Returns False if Nd = 1 (single displacement region) or Nd = 0 (below Ed).
    """
    return ea > 2.0 * ed_val / 0.8


def nd_nrt(ea, ed_val):
    """NRT-DPA: number of atomic displacements per cascade.

    Physical basis: The Norgett-Robinson-Torrens model (1975) modifies the
    original Kinchin-Pease formula with a displacement efficiency g = 0.8,
    obtained from binary collision approximation (BCA) simulations.

    Piecewise definition:
      Nd = 0                          for Ea < Ed
      Nd = 1                          for Ed < Ea < 2Ed/0.8
      Nd = 0.8·Ea / (2·Ed)           for Ea > 2Ed/0.8

    The factor 0.8 accounts for the fact that realistic interatomic scattering
    is NOT hard-sphere-like, reducing the number of displaced atoms by ~20%.

    Reference: Norgett, M.P., Robinson, M.T., Torrens, I.M. (1975).
               Nucl. Eng. Des., 33, 50-54.
    """
    if ea < ed_val:
        return 0
    if ea < 2.0 * ed_val / 0.8:
        return 1
    return 0.8 * ea / (2.0 * ed_val)


def nd_cb(ea, ed_val, z_eff, a_eff):
    """CB-DPA: Chen-Bernard improved model, number of displacements per cascade.

    Physical basis: The NRT model overestimates DPA by ~3x because it ignores
    in-cascade recombination — many Frenkel pairs created in the dense collision
    spike recombine before they can survive as stable defects.

    CB-DPA introduces an athermal recombination efficiency ξ(Ea):
      ξ(Ea) = (2Ed/0.8) / (2Ed/0.8 + β·Ea) + β
      β = Z / (1.5·A)

    where:
      β — asymptotic efficiency at high energies, when sub-cascades form.
          Physically, β ~ 1/3 because at high PKA energies the cascade
          breaks into independent sub-cascades, each with ~1/3 surviving
          fraction relative to NRT.
          1.5 is an empirically determined universal scaling constant.
      The first term describes the decay of efficiency from 1 to β as
      energy increases (from single-cascade to sub-cascade regime).
      The second term β sets the high-energy asymptote.

    Piecewise definition (same as NRT, but with ξ≠1):
      Nd = 0                             for Ea < Ed
      Nd = 1                             for Ed < Ea < 2Ed/0.8
      Nd = (0.8·Ea)/(2·Ed) · ξ(Ea)      for Ea > 2Ed/0.8

    Reference: Chen, S., Bernard, D., Tommasi, J., De Saint Jean, C. (2020).
               EPJ Web of Conferences, 239, 08003.
    """
    if ea < ed_val:
        return 0, {"beta": None, "xi": None}
    if ea < 2.0 * ed_val / 0.8:
        return 1, {"beta": None, "xi": None}

    beta = z_eff / (1.5 * a_eff)
    denom_0 = 2.0 * ed_val / 0.8       # = 2Ed/0.8
    xi = denom_0 / (denom_0 + beta * ea) + beta
    nd = xi * (0.8 * ea / (2.0 * ed_val))

    return nd, {"beta": beta, "xi": xi}


# ═══════════════════════════════════════════════════════════════════════════
#  Output formatting
# ═══════════════════════════════════════════════════════════════════════════

def nd_arcdpa(ea, ed_val, b_arc, c_arc):
    """ARC-DPA (Nordlund et al., 2018): athermal recombination corrected.

    Physical basis: Same as CB-DPA — accounts for in-cascade recombination
    that NRT ignores. Uses a power-law efficiency function fitted to MD
    simulation data for each isotope.

    Efficiency function:
      ξ(Ea) = (1 - c)·(0.8·Ea / 2Ed)^b + c

    where b (typically -0.5 to -1.0) controls the decay rate of efficiency
    from 1 to c, and c (typically ~0.2-0.3) is the asymptotic efficiency at
    high energies. Continuity requires ξ(2Ed/0.8) = 1, which is satisfied
    since (0.8·(2Ed/0.8)/2Ed)^b = 1^b = 1, giving ξ = (1-c)·1 + c = 1.

    Known parameters (Chen et al., 2020 Table 1):
      Fe: b=-0.568, c=0.286    Ni: b=-1.01, c=0.23
      Cu: b=-0.68, c=0.16      Ag: b=-1.06, c=0.257

    Piecewise definition (same structure as NRT/CB-DPA):
      Nd = 0                             for Ea < Ed
      Nd = 1                             for Ed < Ea < 2Ed/0.8
      Nd = (0.8·Ea)/(2·Ed) · ξ(Ea)      for Ea > 2Ed/0.8

    Reference: Nordlund, K., Zinkle, S.J., Sand, A.E., et al. (2018).
               Nat. Commun., 9, 1084.
    """
    if ea < ed_val:
        return 0, {"xi": None, "b_arc": b_arc, "c_arc": c_arc}
    if ea < 2.0 * ed_val / 0.8:
        return 1, {"xi": None, "b_arc": b_arc, "c_arc": c_arc}
    xi = (1.0 - c_arc) * (0.8 * ea / (2.0 * ed_val)) ** b_arc + c_arc
    nd = xi * 0.8 * ea / (2.0 * ed_val)
    return nd, {"xi": xi, "b_arc": b_arc, "c_arc": c_arc}


def print_result(args, ea, part_info, nd_val, ext_info=None):
    """Print formatted DPA calculation result."""
    sep = "=" * 70

    if args.model == "ARC":
        model_label = "ARC-DPA (Nordlund et al. 2018)"
    elif args.model == "CB":
        model_label = "CB-DPA (Chen-Bernard 2020)"
    else:
        model_label = "NRT (Norgett-Robinson-Torrens 1975)"
    part_label  = "Robinson (Lindhard partition)" if args.partition == "robinson" else "Direct (Ea = E_PKA)"

    print(sep)
    print(f"  DPA Calculator — {args.name}")
    print(f"  Model: {model_label}")
    print(f"  Partition function: {part_label}")
    print(f"  Ed = {args.Ed:.0f} eV  |  N_atoms = {args.N_atoms}")
    print(sep)

    # PKA → damage energy
    if args.partition == "robinson":
        print(f"  PKA energy:        {args.E_PKA:,.0f} eV")
        print(f"  Damage energy Ea:  {ea:,.0f} eV  "
              f"(partition fraction f = {part_info['f']:.3f})")
        print(f"  Robinson params:   ε = {part_info['epsilon']:.4f}"
              f"  k = {part_info['k']:.4f}"
              f"  g(ε) = {part_info['g_eps']:.3f}"
              f"  E_L = {part_info['E_L']:,.0f} eV")
    else:
        print(f"  PKA energy = Damage energy Ea = {args.E_PKA:,.0f} eV")

    # DPA computation
    threshold = 2.0 * args.Ed / 0.8
    print(f"{'-'*70}")
    regime = "high-energy" if ea > threshold else ("single-displacement" if ea > args.Ed else "below-Ed")
    print(f"  Ed threshold check:  {ea:,.0f} eV  vs.  {threshold:.0f} eV (2Ed/0.8)")
    print(f"  Regime: {regime}")

    if args.model == "CB" and ext_info and ext_info.get("beta") is not None:
        print(f"  CB-DPA:  β = {ext_info['beta']:.4f}  ξ(Ea) = {ext_info['xi']:.4f}")
    if args.model == "ARC" and ext_info and ext_info.get("xi") is not None:
        print(f"  ARC-DPA:  b_arc = {ext_info['b_arc']:.4f}  c_arc = {ext_info['c_arc']:.4f}  ξ(Ea) = {ext_info['xi']:.4f}")

    print(f"{'-'*70}")
    print(f"  Displacements per cascade (Nd):  {nd_val:,.1f}")
    dpa_per = nd_val / args.N_atoms
    print(f"  DPA per cascade:                 {dpa_per:.6f}")
    if nd_val > 0:
        n_cascades = int(args.target / dpa_per) + 1
        print(f"  Cascades for {args.target:.3f} dpa:  {n_cascades}")
    else:
        print(f"  Nd = 0 — no displacements (Ea < Ed)")

    # References
    print(sep)
    print("  References:")
    print("    NRT model:")
    print("      Norgett, M.P., Robinson, M.T., & Torrens, I.M. (1975).")
    print("      Nucl. Eng. Des., 33, 50-54.")
    if args.model in ("ARC", "CB"):
        print("    ARC-DPA model:")
        print("      Nordlund, K., Zinkle, S.J., Sand, A.E., et al. (2018).")
        print("      Nat. Commun., 9, 1084.")
    if args.model == "CB":
        print("    CB-DPA model:")
        print("      Chen, S., Bernard, D., Tommasi, J., & De Saint Jean, C. (2020).")
        print("      EPJ Web of Conferences, 239, 08003.")
    print("    Robinson partition function:")
    print("      Robinson, M.T. (1970). ORNL-4552.")
    print("    Lindhard stopping theory:")
    print("      Lindhard, J., Scharff, M., & Schiott, H.E. (1963).")
    print("      Kgl. Dan. Vid. Selsk. Mat.-Fys. Medd., 33, 1-42.")
    print(sep)


# ═══════════════════════════════════════════════════════════════════════════
#  Interactive mode
# ═══════════════════════════════════════════════════════════════════════════

def interactive():
    """Prompt the user for all inputs."""
    print("=== DPA Calculator — Interactive Mode ===")
    print(f"Default material: {MATERIAL['name']} (Z_eff={MATERIAL['Z_eff']}, A_eff={MATERIAL['A_eff']})")
    print()

    # ----- Material name -----
    name = input(f"  Material name [{MATERIAL['name']}]: ").strip() or MATERIAL['name']

    # ----- Z_eff -----
    print("Z_eff: Effective atomic number = Σ(c_i·Z_i).")
    print("       Used in Robinson partition and CB-DPA.")
    z_val = float(input(f"  Z_eff [{MATERIAL['Z_eff']}]: ") or str(MATERIAL["Z_eff"]))

    # ----- A_eff -----
    print("A_eff: Effective atomic mass in g/mol = Σ(c_i·A_i).")
    print("       Used in Robinson partition and CB-DPA.")
    a_val = float(input(f"  A_eff [{MATERIAL['A_eff']}]: ") or str(MATERIAL["A_eff"]))

    # ----- PKA energy -----
    print("E_PKA: Primary Knock-on Atom kinetic energy (eV).")
    print("       Energy transferred from incident particle to a lattice atom.")
    epka = float(input("  E_PKA [20000]: ") or "20000")

    # ----- Ed -----
    print("Ed: Threshold displacement energy (eV).")
    print("    Minimum energy to permanently displace a lattice atom.")
    ed_val = float(input("  Ed [62]: ") or "62")

    # ----- Model -----
    print("Model: DPA calculation model.")
    print("  NRT = Norgett-Robinson-Torrens (1975, standard)")
    print("  ARC = Nordlund et al. (2018, MD-fitted recombination correction)")
    print("  CB  = Chen-Bernard (2020, improved, no fitting parameters)")
    model = input("  Model [NRT]: ").strip().upper() or "NRT"

    b_arc = c_arc = None
    if model == "ARC":
        print("b_arc: ARC-DPA power-law exponent from MD fitting.")
        print("       Typical: Fe=-0.568, Ni=-1.01, Cu=-0.68, Ag=-1.06")
        b_arc = float(input("  b_arc: ") or "0")
        print("c_arc: ARC-DPA asymptotic efficiency from MD fitting.")
        print("       Typical: Fe=0.286, Ni=0.23, Cu=0.16, Ag=0.257")
        c_arc = float(input("  c_arc: ") or "0")

    # ----- Partition -----
    print("Partition: how to compute damage energy Ea from E_PKA.")
    print("  direct   = Ea = E_PKA (no electronic loss correction)")
    print("  robinson = Robinson fit to Lindhard theory")
    partition = input("  Partition [direct]: ").strip().lower() or "direct"

    # ----- N_atoms -----
    print("N_atoms: Number of atoms in the simulation cell.")
    n_atoms = int(input("  N_atoms [16000]: ") or "16000")

    # ----- Target -----
    print("Target: Desired DPA dose level.")
    print("        How many cascades are needed to reach this value?")
    target = float(input("  Target dpa [0.5]: ") or "0.5")

    return argparse.Namespace(
        E_PKA=epka, Ed=ed_val, model=model, partition=partition,
        N_atoms=n_atoms, target=target, b_arc=b_arc, c_arc=c_arc,
        name=name, Z=z_val, A=a_val,
        compare=False, test=False,
    )


# ═══════════════════════════════════════════════════════════════════════════
#  Compare all model-partition combinations
# ═══════════════════════════════════════════════════════════════════════════

def compare_all(args):
    """Run all 6 model-partition combinations and print a comparison table."""
    validate_args(args)
    models = [
        ("NRT",     None, None),
        ("ARC",     args.b_arc, args.c_arc),
        ("CB",      None, None),
    ]
    partitions = ["direct", "robinson"]

    # Header
    sep = "=" * 90
    print(sep)
    print(f"  DPA Comparison — {args.name}  |  "
          f"E_PKA = {args.E_PKA:,.0f} eV  |  Ed = {args.Ed:.0f} eV  |  "
          f"N_atoms = {args.N_atoms}")
    print(sep)
    print(f"  {'Model':<8s}  {'Partition':<8s}  {'Ea (eV)':>9s}  "
          f"{'Nd':>8s}  {'DPA/cascade':>12s}  {'Cascades(' + str(args.target) + 'dpa)':>16s}")
    print("  " + "-" * 86)

    for model_name, b_val, c_val in models:
        # Skip ARC if b_arc/c_arc not provided
        arc_skip = False
        if model_name == "ARC" and (b_val is None or c_val is None):
            arc_skip = True

        for part_name in partitions:
            # ── Partition ──
            if part_name == "robinson":
                ea, _ = partition_robinson(args.E_PKA, args.Z, args.A)
            else:
                ea, _ = partition_direct(args.E_PKA)

            # ── DPA model ──
            if arc_skip:
                nd_txt = "(skip)"
                dpa_txt = "(skip)"
                casc_txt = "(b_arc,c_arc needed)"
            else:
                if model_name == "NRT":
                    nd_val = nd_nrt(ea, args.Ed)
                elif model_name == "ARC":
                    nd_val, _ = nd_arcdpa(ea, args.Ed, b_val, c_val)
                else:  # CB
                    nd_val, _ = nd_cb(ea, args.Ed, args.Z, args.A)

                dpa_per = nd_val / args.N_atoms
                cascades = int(args.target / dpa_per) + 1 if nd_val > 0 else 0

                nd_txt = f"{nd_val:,.1f}"
                dpa_txt = f"{dpa_per:.6f}"
                casc_txt = str(cascades)

            md = model_name.rjust(5)
            print(f"  {md:<8s}  {part_name:<8s}  {ea:>9,.0f}  "
                  f"{nd_txt:>8s}  {dpa_txt:>12s}  {casc_txt:>16s}")

    print(sep)
    print("  References:")
    print("    NRT: Norgett, Robinson & Torrens (1975). Nucl. Eng. Des., 33, 50-54.")
    print("    ARC: Nordlund, Zinkle, Sand et al. (2018). Nat. Commun., 9, 1084.")
    print("    CB:  Chen, Bernard, Tommasi, De Saint Jean (2020). EPJ Web Conf., 239, 08003.")
    print("    Robinson: Robinson, M.T. (1970). ORNL-4552.")
    print(sep)


# ═══════════════════════════════════════════════════════════════════════════
#  Input validation
# ═══════════════════════════════════════════════════════════════════════════

def validate_args(args):
    """Check input validity; warn on physically unusual values.

    Returns True if all inputs are valid (program may proceed).
    Exits with an error message if a hard constraint is violated.
    """
    errors = []
    warnings = []

    # ── Hard constraints (program cannot proceed) ─────────────────────
    if args.E_PKA <= 0:
        errors.append(f"E_PKA must be > 0, got {args.E_PKA}")
    if args.Ed <= 0:
        errors.append(f"Ed must be > 0, got {args.Ed}")
    if args.N_atoms <= 0:
        errors.append(f"N_atoms must be > 0, got {args.N_atoms}")
    if args.target <= 0:
        errors.append(f"target dpa must be > 0, got {args.target}")

    if errors:
        print("ERROR — invalid input(s):")
        for e in errors:
            print(f"  {e}")
        sys.exit(1)

    # ── Soft warnings (physically unusual but not impossible) ────────
    if args.E_PKA < 10:
        warnings.append(f"E_PKA = {args.E_PKA:.1f} eV is below typical Ed; "
                         "cascade may not produce any displacements.")
    if args.E_PKA > 5e5:
        warnings.append(f"E_PKA = {args.E_PKA:,.0f} eV exceeds typical PKA "
                         "energy range; Robinson partition may be inaccurate.")
    if args.Ed < 10 or args.Ed > 200:
        warnings.append(f"Ed = {args.Ed:.0f} eV is outside typical range "
                         "(10–200 eV) for most materials.")

    if warnings:
        print("WARNING — physically unusual value(s):")
        for w in warnings:
            print(f"  {w}")
        print()

    return True


# ═══════════════════════════════════════════════════════════════════════════
#  Self-verification tests
# ═══════════════════════════════════════════════════════════════════════════

def run_tests():
    """Run a suite of tests against known reference values.

    Returns True if all tests pass, False otherwise.
    Tests call the mathematical functions directly (no CLI dependency).
    """
    passed = 0
    failed = 0

    def check(desc, actual, expected, tol=1e-6):
        nonlocal passed, failed
        ok = abs(actual - expected) < tol
        if ok:
            passed += 1
            print(f"  PASS  {desc}: {actual:.4f}")
        else:
            failed += 1
            print(f"  FAIL  {desc}: got {actual:.4f}, expected {expected:.4f}")
        return ok

    tol_nd = 0.01   # tolerance for Nd (float rounding, ~0.01 eV sensitivity)

    print("=" * 60)
    print("  Self-verification tests")
    print("=" * 60)

    # ── NRT ──────────────────────────────────────────────────────────────
    # High-energy regime: Nd = 0.8 * 2000 / (2 * 62) = 1600 / 124 = 12.9032...
    check("NRT  Ea=2000 eV Ed=62            ",
          nd_nrt(2000, 62), 12.9032, tol_nd)

    # Below Ed → zero
    check("NRT  Ea=50  eV Ed=62 (0 displace)",
          nd_nrt(50, 62), 0.0, tol_nd)

    # Single displacement: 62 < 100 < 155 (2Ed/0.8 = 155)
    check("NRT  Ea=100 eV Ed=62 (1 displace)",
          nd_nrt(100, 62), 1.0, tol_nd)

    # CB-DPA ────────────────────────────────────────────────────────────
    # β = 49.6 / (1.5 * 118.287) ≈ 0.279546
    # ξ = 155/(155 + 0.279546*20000) + 0.279546
    #   = 155/(155 + 5590.92) + 0.279546
    #   = 155/5745.92 + 0.279546 = 0.02698 + 0.279546 = 0.30653
    # Nd = 0.30653 * 129.032 = 39.55
    nd_cb_val, cb_info = nd_cb(20000, 62, 49.6, 118.287)
    check("CB   β constant (should be ~0.2795) ",
          cb_info["beta"], 0.279546, 1e-5)
    check("CB   ξ(20 keV) (should be ~0.3065)  ",
          cb_info["xi"], 0.3065, 1e-4)
    check("CB   Nd(20 keV) Ed=62               ",
          nd_cb_val, 39.55, tol_nd)

    # CB below Ed
    nd_cb2, _ = nd_cb(50, 62, 49.6, 118.287)
    check("CB   Ea=50  eV Ed=62 (0 displace)   ",
          nd_cb2, 0.0, tol_nd)

    # ── ARC-DPA (Fe parameters: b=-0.568, c=0.286, Ed=40) ───────────────
    # ξ = (1-0.286) * (0.8*20000/80)^(-0.568) + 0.286
    #   = 0.714 * 200.0^(-0.568) + 0.286
    # 200^(-0.568) = exp(-0.568 * ln(200)) = exp(-0.568*5.298) = exp(-3.009) = 0.04934
    # ξ = 0.714*0.04934 + 0.286 = 0.0352 + 0.286 = 0.3212
    # Nd = 0.3212 * (0.8*20000/80) = 0.3212 * 200 = 64.2
    nd_arc_val, arc_info = nd_arcdpa(20000, 40, -0.568, 0.286)
    check("ARC  Fe ξ(20 keV) (should be ~0.321) ",
          arc_info["xi"], 0.3212, 1e-3)
    check("ARC  Fe Nd(20 keV) Ed=40            ",
          nd_arc_val, 64.2, 0.2)

    # ARC below Ed
    nd_arc2, _ = nd_arcdpa(30, 40, -0.568, 0.286)
    check("ARC  Ea=30  eV Ed=40 (0 displace)   ",
          nd_arc2, 0.0, tol_nd)

    # ── Partition ────────────────────────────────────────────────────────
    ea_rob, rob_info = partition_robinson(20000, 49.6, 118.287)
    check("Rob  f(20 keV) < 1 (electronic loss) ",
          float(rob_info["f"] < 1.0), 1.0, 1e-9)
    check("Rob  Ea(20 keV) < E_PKA              ",
          float(ea_rob < 20000), 1.0, 1e-9)
    ea_dir, dir_info = partition_direct(20000)
    check("Dir  f = 1.0, Ea = E_PKA            ",
          ea_dir, 20000.0, tol_nd)

    # ── Piecewise thresholds ─────────────────────────────────────────────
    threshold = 2.0 * 62 / 0.8
    check("Thr  2Ed/0.8 = 155 eV (Ed=62)       ",
          threshold, 155.0, 1e-9)

    # ── Physical ordering: NRT > ARC > CB ─────────────────────────────────
    nd_nrt20 = nd_nrt(20000, 62)
    nd_arc20, _ = nd_arcdpa(20000, 62, -0.568, 0.286)
    nd_cb20, _ = nd_cb(20000, 62, 49.6, 118.287)
    check("Rank Nd(NRT) > Nd(ARC) > Nd(CB)     ",
          float(nd_nrt20 > nd_arc20 > nd_cb20), 1.0, 1e-9)

    # ── Summary ──────────────────────────────────────────────────────────
    total = passed + failed
    print("-" * 60)
    print(f"  Results: {passed} passed, {failed} failed, {total} total")
    print("=" * 60)
    return failed == 0


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="DPA Calculator for HfNbZrTiTa refractory high-entropy alloy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python3 dpa_calculator.py                              # defaults: NRT, direct, 20 keV
  python3 dpa_calculator.py -m CB -p robinson -e 20000   # CB-DPA + Robinson partition
  python3 dpa_calculator.py -c                           # compare all 6 combos
  python3 dpa_calculator.py -c --b_arc -0.568 --c_arc 0.286  # compare incl. ARC
  python3 dpa_calculator.py -i                           # interactive mode
""",
    )

    parser.add_argument(
        "-e", "--E_PKA", type=float, default=20000,
        help="PKA kinetic energy in eV (default: 20000). "
             "This is the energy transferred from the incident particle "
             "(neutron/ion) to a lattice atom, setting it in motion.")
    parser.add_argument(
        "-d", "--Ed", type=float, default=62,
        help="Threshold displacement energy in eV (default: 62). "
             "Minimum kinetic energy for a lattice atom to be permanently "
             "displaced and create a stable Frenkel pair.")
    parser.add_argument(
        "-m", "--model", choices=["NRT", "ARC", "CB"], default="NRT",
        help="DPA model: NRT (Norgett 1975, standard), "
             "ARC (Nordlund 2018, MD-fitted recombination correction), or "
             "CB (Chen-Bernard 2020, no fitting parameters).")
    parser.add_argument(
        "--b_arc", type=float, default=None,
        help="ARC-DPA exponent parameter b_arc (required when --model ARC). "
             "Typical: Fe=-0.568, Ni=-1.01, Cu=-0.68, Ag=-1.06.")
    parser.add_argument(
        "--c_arc", type=float, default=None,
        help="ARC-DPA asymptotic parameter c_arc (required when --model ARC). "
             "Typical: Fe=0.286, Ni=0.23, Cu=0.16, Ag=0.257.")
    parser.add_argument(
        "-p", "--partition", choices=["direct", "robinson"], default="direct",
        help="Partition function: direct (Ea = E_PKA) or "
             "robinson (Lindhard energy partition with electronic loss correction).")
    parser.add_argument(
        "-n", "--N_atoms", type=int, default=16000,
        help="Number of atoms in the simulation cell (default: 16000).")
    parser.add_argument(
        "-t", "--target", type=float, default=0.5,
        help="Target DPA dose level (default: 0.5). "
             "Used to compute the number of cascades required.")
    parser.add_argument(
        "--name", default=MATERIAL["name"],
        help="Material name for display (default: HfNbZrTiTa). "
             "Does not affect calculations — only the output header.")
    parser.add_argument(
        "--Z", type=float, default=MATERIAL["Z_eff"],
        help="Effective atomic number Z_eff (default: 49.6 for HfNbZrTiTa). "
             "Used in Robinson partition and CB-DPA model.")
    parser.add_argument(
        "--A", type=float, default=MATERIAL["A_eff"],
        help="Effective atomic mass A_eff in g/mol (default: 118.287 for HfNbZrTiTa). "
             "Used in Robinson partition and CB-DPA model.")
    parser.add_argument(
        "-i", "--interactive", action="store_true",
        help="Run in interactive mode (prompt for each input).")
    parser.add_argument(
        "-c", "--compare", action="store_true",
        help="Run all 6 model-partition combinations and print a comparison table. "
             "Overrides --model and --partition. For ARC-DPA, use --b_arc --c_arc.")
    parser.add_argument(
        "--test", action="store_true",
        help="Run self-verification tests against known reference values.")

    args = parser.parse_args()
    if args.interactive:
        ict = interactive()
        ict.compare = args.compare
        ict.test = args.test
        args = ict

    validate_args(args)

    # ── Compare mode: run all combinations ───────────────────────────────────
    if args.compare:
        compare_all(args)
        return

    # ── Self-test mode ───────────────────────────────────────────────────────
    if args.test:
        success = run_tests()
        sys.exit(0 if success else 1)

    # ── Step 1: compute damage energy Ea ───────────────────────────────────
    if args.partition == "robinson":
        ea, part_info = partition_robinson(args.E_PKA, args.Z, args.A)
    else:
        ea, part_info = partition_direct(args.E_PKA)

    # ── Step 2: compute displacements ──────────────────────────────────────
    ext_info = None
    if args.model == "ARC":
        if args.b_arc is None or args.c_arc is None:
            print("ERROR: --model ARC requires --b_arc and --c_arc parameters.")
            print("  Example: --b_arc -0.568 --c_arc 0.286")
            sys.exit(1)
        nd_val, ext_info = nd_arcdpa(ea, args.Ed, args.b_arc, args.c_arc)
    elif args.model == "CB":
        nd_val, ext_info = nd_cb(ea, args.Ed, args.Z, args.A)
    else:
        nd_val = nd_nrt(ea, args.Ed)

    # ── Step 3: print ──────────────────────────────────────────────────────
    print_result(args, ea, part_info, nd_val, ext_info)


if __name__ == "__main__":
    main()
