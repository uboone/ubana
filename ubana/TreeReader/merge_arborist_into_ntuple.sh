#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  ./merge_arborist_into_ntuple.sh <arborist_out.root> <source_ntuple.root> <output.root>

Example:
  ./merge_arborist_into_ntuple.sh arborist_20251219T195436_out.root random_ncpi0_ntuple_file.root random_ncpi0_with_arborist.root

What it does:
  - Copies <source_ntuple.root> -> <output.root> (preserves all original directories/objects)
  - Adds/overwrites a top-level TTree named:
      spline_weights
    with branches:
      run (int), subrun (int), event (int), entry (long64)
      one branch per mcweight key (each is vector<double>)
    filled from Arborist's arborist/eventweight_tree branch: mcweight.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -ne 3 ]]; then
  usage
  exit 2
fi

ARB_FILE="$1"
SRC_FILE="$2"
OUT_FILE="$3"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Look for the macro in the script dir first (dev environment),
# then fall back to the installed source location (via UBANA_DIR).
if [[ -f "${SCRIPT_DIR}/merge_arborist_into_ntuple.C" ]]; then
  MACRO="${SCRIPT_DIR}/merge_arborist_into_ntuple.C"
elif [[ -n "${UBANA_DIR:-}" && -f "${UBANA_DIR}/source/ubana/TreeReader/merge_arborist_into_ntuple.C" ]]; then
  MACRO="${UBANA_DIR}/source/ubana/TreeReader/merge_arborist_into_ntuple.C"
else
  echo "ERROR: Cannot find merge_arborist_into_ntuple.C" >&2
  echo "Looked in: ${SCRIPT_DIR}/" >&2
  [[ -n "${UBANA_DIR:-}" ]] && echo "       and: ${UBANA_DIR}/source/ubana/TreeReader/" >&2
  exit 3
fi

root -l -b <<EOF
.L ${MACRO}
merge_arborist_into_ntuple("${ARB_FILE}","${SRC_FILE}","${OUT_FILE}");
.q
EOF


