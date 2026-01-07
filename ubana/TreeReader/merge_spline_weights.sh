#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  (jobsub endscript mode, no args)
    ./merge_spline_weights.sh

  (explicit mode)
    ./merge_spline_weights.sh <treereader_artroot.root> <source_ntuple.root> <output_ntuple.root> [run_arborist.fcl]

What it does (single endscript):
  1) Runs Arborist on the TreeReader/EventWeight artroot file to produce:
       arborist_*_out.root
     (by default in the current directory)
  2) Runs merge_arborist_into_ntuple.sh to copy <source_ntuple.root> -> <output_ntuple.root>
     and add a top-level TTree named "spline_weights".

Notes:
  - You must run this in an environment where `lar` and `root` are available and
    `run_arborist.fcl` is resolvable (either via FHICL_FILE_PATH or by providing
    the optional 4th argument with an explicit path).
  - In jobsub endscript mode, the script auto-detects inputs in the current working
    directory:
      - treereader_artroot: newest file matching treereader*.root
      - source_ntuple: a .root file containing the 'singlephotonana' directory, excluding treereader*/arborist*
      - output_ntuple: ${source_ntuple%.root}_with_splines.root unless SPLINE_MERGE_OUTPUT is set
    You can override auto-detection with environment variables:
      - TREEREADER_ARTROOT
      - SOURCE_NTUPLE
      - SPLINE_MERGE_OUTPUT
      - ARBORIST_FCL
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MERGE_SCRIPT="${SCRIPT_DIR}/merge_arborist_into_ntuple.sh"

if [[ ! -x "${MERGE_SCRIPT}" ]]; then
  echo "ERROR: merge script not found/executable: ${MERGE_SCRIPT}" >&2
  exit 3
fi

detect_ntuple_file() {
  # Return first file that looks like the original ntuple:
  # has TDirectory "singlephotonana" and isn't an arborist/treereader output.
  for f in *.root; do
    [[ -f "$f" ]] || continue
    [[ "$f" == treereader* ]] && continue
    [[ "$f" == arborist_*_out.root ]] && continue
    # Check for "singlephotonana" directory.
    if root -l -b -q <<EOF >/dev/null 2>&1
TFile tf("$f");
if (tf.IsZombie()) { gSystem->Exit(1); }
auto* d = tf.GetDirectory("singlephotonana");
gSystem->Exit(d ? 0 : 2);
EOF
    then
      echo "$f"
      return 0
    fi
  done
  return 1
}

detect_treereader_file() {
  # Prefer newest treereader*.root in cwd.
  local f
  f="$(ls -1t treereader*.root 2>/dev/null | head -n 1 || true)"
  [[ -n "$f" ]] && echo "$f"
}

if [[ $# -eq 0 ]]; then
  TREEREADER_ARTROOT="${TREEREADER_ARTROOT:-$(detect_treereader_file)}"
  SOURCE_NTUPLE="${SOURCE_NTUPLE:-$(detect_ntuple_file || true)}"
  ARBORIST_FCL="${ARBORIST_FCL:-run_arborist.fcl}"
  if [[ -z "${SOURCE_NTUPLE}" ]]; then
    echo "ERROR: could not auto-detect source ntuple in $(pwd)" >&2
    echo "Set SOURCE_NTUPLE=/path/to/ntuple.root (or run in explicit mode)." >&2
    exit 2
  fi
  OUTPUT_NTUPLE="${SPLINE_MERGE_OUTPUT:-${SOURCE_NTUPLE%.root}_with_splines.root}"
elif [[ $# -eq 3 || $# -eq 4 ]]; then
  TREEREADER_ARTROOT="$1"
  SOURCE_NTUPLE="$2"
  OUTPUT_NTUPLE="$3"
  ARBORIST_FCL="${4:-run_arborist.fcl}"
else
  usage
  exit 2
fi

if [[ -z "${TREEREADER_ARTROOT}" || ! -f "${TREEREADER_ARTROOT}" ]]; then
  echo "ERROR: treereader artroot file not found: ${TREEREADER_ARTROOT:-<empty>}" >&2
  echo "Tip: set TREEREADER_ARTROOT=/path/to/treereader*.root (or pass args in explicit mode)." >&2
  exit 4
fi
if [[ ! -f "${SOURCE_NTUPLE}" ]]; then
  echo "ERROR: source ntuple file not found: ${SOURCE_NTUPLE}" >&2
  exit 5
fi

echo "Treereader artroot: ${TREEREADER_ARTROOT}"
echo "Source ntuple:      ${SOURCE_NTUPLE}"
echo "Output ntuple:      ${OUTPUT_NTUPLE}"
echo "Arborist fcl:       ${ARBORIST_FCL}"

echo "Running Arborist in: $(pwd)"

# Arborist writes to a local relative filename via TFileService:
#   arborist_%tc_out.root
lar -c "${ARBORIST_FCL}" "${TREEREADER_ARTROOT}"

ARB_OUT="$(ls -1t arborist_*_out.root 2>/dev/null | head -n 1 || true)"
if [[ -z "${ARB_OUT}" ]]; then
  echo "ERROR: Arborist did not produce arborist_*_out.root in $(pwd)" >&2
  echo "Tip: check that ${ARBORIST_FCL} sets TFileService.fileName to arborist_%tc_out.root" >&2
  exit 6
fi

ARB_OUT_ABS="$(pwd)/${ARB_OUT}"
echo "Arborist output:   ${ARB_OUT_ABS}"
echo "Merging into output ntuple: ${OUTPUT_NTUPLE}"

"${MERGE_SCRIPT}" "${ARB_OUT_ABS}" "${SOURCE_NTUPLE}" "${OUTPUT_NTUPLE}"

echo "Done."


