#!/bin/bash
#==============================================================================
# run_Ed.sh - Master script for threshold displacement energy (Ed) calculation
#
# Method: Byggmästar et al. (2024) with binary search energy sweep
# Models: MA (meta-atom), RSS (random solid solution), LCO (short-range ordered)
# Features: multi-config, per-element Ed, 4 defect methods, parallel execution
#
# Usage:  bash run_Ed.sh
#         Edit parameters below and alloy.conf before running.
#==============================================================================

##=================== User Parameters (edit as needed) ========================##
NCORE=16
NJOB=4
LMP=~/lmp-2023                    # LAMMPS executable path
MPIRUN="mpirun"                   # mpirun executable path
CONFIG_SOURCE="auto"               # "auto" = generate random alloy; "custom" = read data files
NCONFIG=5
NDIR_PER_CONFIG=100
SIMTIME=6.0
EMIN=10
EMAX=180
ESTEP=1
TEMP=40
DEFECT_METHOD=ovito
VERBOSE=0                         # 0=quiet (completions + periodic progress), 1=full binary search detail
##============================================================================##


##=================== Setup ===================================================##
set -u
set -m  # enable job control: background jobs get own process group (Ctrl+C -> foreground only)
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR" || { echo "ERROR: cannot cd to $WORKDIR"; exit 1; }
mkdir -p ../results logs
DUMPDIR="dump_files"
mkdir -p ${DUMPDIR}

##-- Cleanup trap: Ctrl+C kills all background jobs and mpirun processes
cleanup() {
    echo ""
    echo ">>> Interrupted. Stopping all jobs..."
    killall -KILL lmp-2023 mpirun orted 2>/dev/null
    sleep 1
    for pid in "${ALL_PIDS[@]:-}"; do
        kill -KILL $pid 2>/dev/null
    done
    killall -KILL lmp-2023 mpirun orted 2>/dev/null
    echo ">>> All processes stopped."
    exit 1
}
trap cleanup SIGINT SIGTERM

## Source environment configuration (overrides LMP and MPIRUN defaults)
if [ -f ../env.conf ]; then
    source ../env.conf
fi

##-- Generate RSS substitution include file (if CONFIG_SOURCE=auto and N_ELEMENTS>1)
generate_rss_include() {
    local n=$1
    local seed=$2
    local out="rss_substitute.lmp"
    > "$out"   # truncate/create empty
    if [ "$n" -le 1 ]; then
        return  # single element, no substitution needed
    fi
    local denom=$n
    local i=2
    while [ "$i" -le "$n" ]; do
        local s=$((seed + i - 2))
        echo "variable f_${i} equal 1/${denom}" >> "$out"
        echo "set type 1 type/fraction ${i} \${f_${i}} ${s}" >> "$out"
        i=$((i + 1))
        denom=$((denom - 1))
    done
}
export -f generate_rss_include

## Source alloy configuration
if [ -f alloy.conf ]; then
    source alloy.conf
    ELEMENTS="${ELEMENT_NAMES}"
    echo "  Alloy config loaded: ${ALLOY_NAME}"
else
    echo "WARNING: alloy.conf not found, using hardcoded defaults."
    ALLOY_NAME="HfNbZrTiTa"
    N_ELEMENTS=5
    POTENTIAL_FILE="HfNbZrTiTa.eam.alloy.0.zbl"
    ELEMENTS="Hf Nb Zr Ti Ta"
    LATTICE_CONSTANT=3.40
    LCO_DATA_PREFIX="data.mc"
    CONFIG_SOURCE="auto"
    CRYSTAL_STRUCTURE="bcc"
fi

## Use potential and elements from alloy config
POTENTIAL="${POTENTIAL_FILE}"
ELEMENTS_PAIR="${ELEMENTS}"

TOTAL_DIRS=$((NCONFIG * NDIR_PER_CONFIG))
NUMPROC=$((NJOB * NCORE))

## Element type -> name mapping for display
declare -a ELEM_NAME=()
read -ra ELEM_NAME <<< "${ELEMENTS}"

echo "================================================"
echo "  Ed Calculation - Byggmästar (2024) method"
echo "  Alloy: ${ALLOY_NAME}  |  Configs: ${NCONFIG}  |  Dir/Config: ${NDIR_PER_CONFIG}"
echo "  Total directions: ${TOTAL_DIRS}"
echo "  NJOB=${NJOB}  NCORE=${NCORE}  |  Total cores: ${NUMPROC}"
echo "  Defect method: ${DEFECT_METHOD}"
echo "  Sim time: ${SIMTIME} ps  |  E range: ${EMIN}-${EMAX} eV"
echo "================================================"


##=================== Helper Functions ========================================##

##-- Element type number -> symbol
elem_name() {
    local t=$1
    if [ "${t}" -ge 1 ] && [ "${t}" -le "${#ELEM_NAME[@]}" ] 2>/dev/null; then
        echo "${ELEM_NAME[$((t-1))]}"
    else
        echo "t${t}"
    fi
}

##-- Run a single recoil simulation and check for defects
run_recoil() {
    local energy=$1 dx=$2 dy=$3 dz=$4 idx=$5 seed=$6 config=$7
    local logfile="logs/log.c${config}.d${idx}.e${energy}"
    local outfile="logs/out.c${config}.d${idx}.e${energy}"

    ${MPIRUN} -n ${NCORE} ${LMP} \
        -log ${logfile} \
        -var EPKA ${energy} \
        -var DX ${dx} -var DY ${dy} -var DZ ${dz} \
        -var TEMP ${TEMP} \
        -var SIMTIME ${SIMTIME} \
        -var SEED ${seed} \
        -var DIR_IDX ${idx} \
        -var CONFIG ${config} \
        -var DUMPDIR ${DUMPDIR} \
        -var POTENTIAL ${POTENTIAL} \
        -var ELEMENTS "${ELEMENTS_PAIR}" \
        -var CRYSTAL_STRUCTURE ${CRYSTAL_STRUCTURE} \
        < in.ed.recoil.lmp > ${outfile} 2>&1

    local ndefects result
    case "${DEFECT_METHOD}" in
        ovito)
            result=$(ovitos check_defects.py \
                "${DUMPDIR}/dump.recoil.c${config}.d${idx}.e${energy}.final.gz" \
                "${DUMPDIR}/dump.recoil.c${config}.d${idx}.e${energy}.init.gz" 2>/dev/null)
            ;;
        cna)
            ndefects=$(grep "CNA_DEFECTS:" ${outfile} | tail -1 | awk '{print $2}')
            if [ -z "${ndefects}" ]; then result="ERROR"
            elif [ "${ndefects}" -gt 0 ]; then result="DEFECT"
            else result="NO_DEFECT"; fi
            ;;
        ptm)
            ndefects=$(grep "PTM_DEFECTS:" ${outfile} | tail -1 | awk '{print $2}')
            if [ -z "${ndefects}" ]; then result="ERROR"
            elif [ "${ndefects}" -gt 0 ]; then result="DEFECT"
            else result="NO_DEFECT"; fi
            ;;
        *)
            ndefects=$(grep "NDEFECTS:" ${outfile} | tail -1 | awk '{print $2}')
            if [ -z "${ndefects}" ]; then result="ERROR"
            elif [ "${ndefects}" -gt 0 ]; then result="DEFECT"
            else result="NO_DEFECT"; fi
            ;;
    esac

    echo "${result}"
}

##-- Binary search for Ed in a given direction
binary_search_ed() {
    local dx=$1 dy=$2 dz=$3 idx=$4 seed=$5 config=$6
    local e_low=${EMIN} e_high=${EMAX} e_mid result

    while [ $((e_high - e_low)) -gt ${ESTEP} ]; do
        e_mid=$(( ((e_low + e_high) / 2) / ESTEP * ESTEP ))
        [ ${e_mid} -le ${e_low} ]  && e_mid=$((e_low + ESTEP))
        [ ${e_mid} -ge ${e_high} ] && e_mid=$((e_high - ESTEP))

        result=$(run_recoil ${e_mid} "${dx}" "${dy}" "${dz}" ${idx} ${seed} ${config})

        case "${result}" in
            DEFECT)    e_high=${e_mid} ;;
            NO_DEFECT) e_low=${e_mid} ;;
            *)         echo "ERROR" >&2; return 1 ;;
        esac

        [ "${VERBOSE}" = "1" ] && echo "  c${config}:d${idx} E=${e_mid} eV -> ${result}" >&2
    done

    echo ${e_high}
}

export -f run_recoil binary_search_ed

##-- Poll and remove finished PIDs from a list
clean_pids() {
    local name=$1
    local -n __arr=$2
    local new=()
    for p in "${__arr[@]}"; do
        kill -0 $p 2>/dev/null && new+=($p)
    done
    __arr=("${new[@]}")
}

##-- Count active PIDs
count_active() {
    local -n __arr=$1
    local c=0
    for p in "${__arr[@]}"; do
        kill -0 $p 2>/dev/null && c=$((c + 1))
    done
    echo $c
}


##=================== Main Config Loop ========================================##

GLOBAL_START_TIME=$(date +%s)
ALL_PIDS=()

for CONFIG in $(seq 1 ${NCONFIG}); do
    CONFIG_START_TIME=$(date +%s)
    echo ""
    echo "========================================================================"
    echo "  CONFIG ${CONFIG}/${NCONFIG}  |  Alloy: ${ALLOY_NAME}  |  Source: ${CONFIG_SOURCE}"
    echo "========================================================================"

    CONFIG_RESULTS="../results/config_${CONFIG}"
    mkdir -p ${CONFIG_RESULTS}

    ##--- Step 1: Equilibrate this config
    echo ""
    echo "--- Config ${CONFIG}: Step 1 - Equilibration ---"

    EQUIL_DATA="data.equilibrated.config_${CONFIG}"
    EQUIL_MODEL_TAG=".equilibrated_model.config_${CONFIG}"

    if [ -f "${EQUIL_DATA}" ] && [ -f "${EQUIL_MODEL_TAG}" ]; then
        prev_source=$(cat ${EQUIL_MODEL_TAG})
        if [ "${prev_source}" != "${CONFIG_SOURCE}" ]; then
            echo "WARNING: ${EQUIL_DATA} was for '${prev_source}', re-equilibrating..."
            rm -f ${EQUIL_DATA} ${EQUIL_MODEL_TAG}
        fi
    fi

    if [ ! -f "${EQUIL_DATA}" ]; then
        echo "  Equilibrating ${ALLOY_NAME} config ${CONFIG} at ${TEMP} K ..."
        LCO_DATA=""
        RSS_SEED=$((10000 + CONFIG * 300))
        if [ "${CONFIG_SOURCE}" == "custom" ]; then
            LCO_DATA="${LCO_DATA_PREFIX}.${CONFIG}"
        fi
        if [ "${CONFIG_SOURCE}" != "custom" ]; then
            generate_rss_include ${N_ELEMENTS} ${RSS_SEED}
        fi
        EQUIL_LOG="logs/log.equilibrate.c${CONFIG}"
        EQUIL_OUT="logs/out.equilibrate.c${CONFIG}"

        ${MPIRUN} -n ${NCORE} ${LMP} \
            -log ${EQUIL_LOG} \
            -var CONFIG_SOURCE ${CONFIG_SOURCE} \
            -var N_ELEMENTS ${N_ELEMENTS} \
            -var TEMP ${TEMP} \
            -var CONFIG ${CONFIG} \
            -var DATA_OUT ${EQUIL_DATA} \
            -var LCO_DATA "${LCO_DATA}" \
            -var POTENTIAL ${POTENTIAL} \
            -var ELEMENTS "${ELEMENTS_PAIR}" \
            -var LATTICE ${LATTICE_CONSTANT} \
            -var CRYSTAL_STRUCTURE ${CRYSTAL_STRUCTURE} \
            -var ALLOY_NAME ${ALLOY_NAME} \
            < in.ed.equilibrate.lmp > ${EQUIL_OUT} 2>&1

        if [ ! -f "${EQUIL_DATA}" ]; then
            echo "ERROR: Equilibration failed for config ${CONFIG}. Check ${EQUIL_OUT}"
            exit 1
        fi
        echo "${CONFIG_SOURCE}" > ${EQUIL_MODEL_TAG}
        echo "  Equilibration complete -> ${EQUIL_DATA}"
    else
        echo "  ${EQUIL_DATA} found, skipping."
    fi

    ##--- Step 2: Generate directions
    echo ""
    echo "--- Config ${CONFIG}: Step 2 - Directions ---"

    DIRFILE="directions.config_${CONFIG}.txt"
    if [ ! -f "${DIRFILE}" ]; then
        python3 generate_directions.py ${NDIR_PER_CONFIG} $((42 + CONFIG * 100)) > ${DIRFILE}
        echo "  Directions generated -> ${DIRFILE}"
    else
        echo "  ${DIRFILE} found, skipping."
        local_ndir=$(wc -l < ${DIRFILE})
        echo "  (${local_ndir} directions available)"
    fi

    ##--- Step 3: Parallel energy sweep
    echo ""
    echo "--- Config ${CONFIG}: Step 3 - Ed calculation (${NDIR_PER_CONFIG} dir, NJOB=${NJOB}) ---"
    echo ""

    ## Collect pending directions
    declare -a PENDING=()
    while read -r idx dx dy dz; do
        [ -z "${idx:-}" ] && continue
        [ -f "${CONFIG_RESULTS}/Ed_direction_${idx}.txt" ] && continue
        PENDING+=("${idx}|${dx}|${dy}|${dz}")
    done < ${DIRFILE}

    TOTAL_PENDING=${#PENDING[@]}
    echo "  Pending: ${TOTAL_PENDING} directions"
    if [ ${TOTAL_PENDING} -eq 0 ]; then
        echo "  All directions already complete."
        CONFIG_ELAPSED=0
        continue
    fi

    LAUNCHED=0
    CONFIG_PIDS=()
    LAST_REPORTED=0

    while [ ${LAUNCHED} -lt ${TOTAL_PENDING} ]; do
        ## Wait for a slot
        while true; do
            clean_pids CONFIG_PIDS CONFIG_PIDS
            [ ${#CONFIG_PIDS[@]} -lt ${NJOB} ] && break
            sleep 0.3
        done

        ## Launch next pending direction
        IFS='|' read idx dx dy dz <<< "${PENDING[$LAUNCHED]}"
        LAUNCHED=$((LAUNCHED + 1))
        seed=$((1000 + CONFIG * 10000 + idx * 7))

        (
            t0=$(date +%s)
            ed=$(binary_search_ed "${dx}" "${dy}" "${dz}" "${idx}" "${seed}" "${CONFIG}")
            elapsed=$(( $(date +%s) - t0 ))
            pka_type=$(grep -h "PKA_TYPE:" logs/out.c${CONFIG}.d${idx}.e* 2>/dev/null | tail -1 | awk '{print $2}')
            pka_type=${pka_type:-0}
            pka_elem=$(elem_name ${pka_type})
            if [ -f "${CONFIG_RESULTS}/Ed_direction_${idx}.txt" ]; then
                exit 0
            fi
            if [ "${ed}" = "ERROR" ] || [ -z "${ed}" ]; then
                echo "ERROR ${pka_type}" > ${CONFIG_RESULTS}/Ed_direction_${idx}.txt
                echo "  [c${CONFIG}:d${idx}]  ERROR  (${elapsed}s)" >&2
            else
                echo "${ed} ${pka_type}" > ${CONFIG_RESULTS}/Ed_direction_${idx}.txt
                echo "  [c${CONFIG}:d${idx}]  Ed=${ed} eV  PKA=${pka_elem}  (${elapsed}s)" >&2
            fi
        ) &
        CONFIG_PIDS+=($!)
        ALL_PIDS+=($!)

        ## Progress: print only when done count crosses a 10-boundary or first/last
        DONE_COUNT=$(ls ${CONFIG_RESULTS}/Ed_direction_*.txt 2>/dev/null | wc -l)
        if [ ${DONE_COUNT} -gt ${LAST_REPORTED} ] && \
           ( [ $((DONE_COUNT % 10)) -eq 0 ] || [ ${DONE_COUNT} -eq ${TOTAL_PENDING} ] || [ ${DONE_COUNT} -le 2 ] ); then
            LAST_REPORTED=${DONE_COUNT}
            ELAPSED=$(( $(date +%s) - GLOBAL_START_TIME ))
            ELAPSED_FMT=$(printf '%02d:%02d:%02d' $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60)))
            if [ ${DONE_COUNT} -gt 0 ]; then
                AVG_TIME=$(( ELAPSED / DONE_COUNT ))
                REMAINING=$(( AVG_TIME * (TOTAL_PENDING - DONE_COUNT) ))
                REMAINING_FMT=$(printf '%02d:%02d:%02d' $((REMAINING/3600)) $(((REMAINING%3600)/60)) $((REMAINING%60)))
                echo "  Config ${CONFIG}: ${DONE_COUNT}/${TOTAL_PENDING} ($((100 * DONE_COUNT / TOTAL_PENDING))%) done  |  ${#CONFIG_PIDS[@]}/${NJOB} running  |  ETA ${REMAINING_FMT}"
            else
                echo "  Config ${CONFIG}: ${DONE_COUNT}/${TOTAL_PENDING}  |  ${#CONFIG_PIDS[@]}/${NJOB} running  |  Elapsed ${ELAPSED_FMT}"
            fi
        fi
    done

    ## Wait for all jobs in this config to finish
    for pid in "${CONFIG_PIDS[@]}"; do
        wait $pid 2>/dev/null || true
    done
    DONE_COUNT=$(ls ${CONFIG_RESULTS}/Ed_direction_*.txt 2>/dev/null | wc -l)
    echo "  Config ${CONFIG}: ${DONE_COUNT}/${TOTAL_PENDING} (100%) done  |  0/${NJOB} running"

    CONFIG_ELAPSED=$(( $(date +%s) - CONFIG_START_TIME ))
    CONFIG_FMT=$(printf '%02d:%02d:%02d' $((CONFIG_ELAPSED/3600)) $(((CONFIG_ELAPSED%3600)/60)) $((CONFIG_ELAPSED%60)))
    echo "  Config ${CONFIG} completed in ${CONFIG_FMT}"

done

##=================== Step 4: Collect results =================================##
echo ""
echo "========================================================================"
echo "  Collecting results across all configurations..."
echo "========================================================================"
echo ""

python3 collect_results.py ../results/ "${ELEMENTS}" > ../Ed_summary.txt
cat ../Ed_summary.txt

echo ""
echo "================================================"
echo "  Ed calculation complete."
echo "  Results: ../Ed_summary.txt"
echo "  Per-config: ../results/config_*/Ed_direction_*.txt"
echo "================================================"
