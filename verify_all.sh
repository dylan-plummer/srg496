#!/usr/bin/env bash
# Quick verification of srg(496,231,102,112) and its complement (about 2 minutes).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
echo "== graph files"
(cd graphs && sha256sum -c SHA256SUMS)
echo "== original constructor and standard-library audit (seed/)"
WORK=$(mktemp -d)
cp seed/construct.py seed/audit.py data/starter.tsv "$WORK"/
(cd "$WORK" && $PY construct.py && $PY audit.py)
cmp "$WORK"/srg_496_231_102_112.g6 graphs/srg_496_231_102_112.g6
cmp "$WORK"/srg_496_264_144_136.g6 graphs/srg_496_264_144_136.g6
cmp "$WORK"/steiner_2_61_5.txt data/steiner_2_61_5.txt
rm -rf "$WORK"
echo "== independent rebuild from the written description"
$PY audit/independent_rebuild.py
echo "== Steiner half-frame"
$PY audit/steiner_subframe.py
echo "== three random family members (srg check + nauty)"
$PY family/family.py 3 11
echo "== fiber recovery on the exported graph"
$PY family/fibre_planes.py graphs/srg_496_231_102_112.g6
echo "== template scan"
$PY template/scan_template.py > /dev/null && echo "scan_template OK"
echo "== three designs, random members in all rows, recovery under switching, the lower bound"
(cd extensions && $PY designs_audit.py)
(cd extensions && sha256sum -c SHA256SUMS)
echo "ALL CHECKS PASSED"
# Longer runs:
#   python seed/symmetry.py             # |Aut| = 4 (nauty and bliss)
#   python audit/twograph.py            # two-graph identification (nauty on the double covers)
#   python family/family.py 50 2026     # the sample in family/family_sample.out
#   python family/fibre_planes_family.py 6 7
#   (cd extensions && python seeds_audit.py)        # order-84 seeds, all exported graphs, Burnside (several minutes)
#   (cd extensions && python starters/recover.py && python starters/plane_count.py && python starters/blind_audit.py)
#   companion reproduction of the starter census and the catalogue: see extensions/CATALOGUE.md
