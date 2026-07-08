#!/usr/bin/env ovitos
"""Check for Frenkel pairs using OVITO Wigner-Seitz analysis.

Usage: ovitos check_defects.py <final_dump> <reference_dump>
       (or run via ./check_defects.py after chmod +x)
Output: DEFECT or NO_DEFECT (printed to stdout)
"""
import sys
from ovito.io import import_file
from ovito.modifiers import WignerSeitzAnalysisModifier
from ovito.pipeline import FileSource

final_file = sys.argv[1]
ref_file = sys.argv[2]

try:
    pipeline = import_file(final_file)
    ws = WignerSeitzAnalysisModifier()
    ws.reference = FileSource()
    ws.reference.load(ref_file)
    pipeline.modifiers.append(ws)
    data = pipeline.compute()
    ni = data.attributes['WignerSeitz.interstitial_count']
    nv = data.attributes['WignerSeitz.vacancy_count']
    if ni > 0 or nv > 0:
        print("DEFECT")
    else:
        print("NO_DEFECT")
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)
