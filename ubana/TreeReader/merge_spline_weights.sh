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
      - OUTDIR (pnfs directory to copy output to; auto-detected if not set)
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

fetch_consumed_file() {
  # If consumed_files.list exists, try to fetch the input file via ifdh.
  local consumed_list="consumed_files.list"
  if [[ ! -f "${consumed_list}" ]]; then
    echo "DEBUG: No consumed_files.list found" >&2
    return 1
  fi
  
  local sam_filename
  sam_filename="$(head -n 1 "${consumed_list}" | tr -d '[:space:]')"
  if [[ -z "${sam_filename}" ]]; then
    echo "DEBUG: consumed_files.list is empty" >&2
    return 1
  fi
  
  echo "DEBUG: Fetching consumed file from SAM: ${sam_filename}" >&2
  
  # Try to get the file location and copy it locally
  if command -v ifdh >/dev/null 2>&1; then
    local sam_uri
    sam_uri="$(samweb locate-file "${sam_filename}" 2>/dev/null | head -n 1 || true)"
    if [[ -n "${sam_uri}" ]]; then
      echo "DEBUG: SAM location: ${sam_uri}" >&2
      # Convert enstore:/pnfs/... or similar to just the path
      local file_path="${sam_uri#enstore:}"
      file_path="${file_path#dcache:}"
      file_path="${file_path%%\(*}"  # Remove any trailing (...)
      
      if [[ -n "${file_path}" ]]; then
        echo "DEBUG: Copying ${file_path}/${sam_filename} to $(pwd)/" >&2
        if ifdh cp "${file_path}/${sam_filename}" "./${sam_filename}" >&2 2>&1; then
          echo "DEBUG: Successfully fetched ${sam_filename}" >&2
          echo "${sam_filename}"
          return 0
        else
          echo "DEBUG: ifdh cp failed, trying direct path" >&2
        fi
      fi
    fi
    
    # Fallback: try ifdh fetchInput if available
    echo "DEBUG: Trying ifdh fetchInput for ${sam_filename}" >&2
    if ifdh fetchInput "${sam_filename}" >&2 2>&1; then
      echo "DEBUG: fetchInput succeeded" >&2
      echo "${sam_filename}"
      return 0
    fi
  fi
  
  echo "DEBUG: Could not fetch consumed file" >&2
  return 1
}

detect_ntuple_file() {
  # Return first file that looks like the original ntuple:
  # has TDirectory "singlephotonana" and isn't an arborist/treereader output.
  echo "DEBUG: Looking for ntuple files in $(pwd):" >&2
  echo "DEBUG: All files in directory:" >&2
  ls -la >&2 || true
  echo "DEBUG: ROOT files only:" >&2
  ls -la *.root >&2 2>&1 || echo "DEBUG: No .root files found" >&2
  for f in *.root; do
    [[ -f "$f" ]] || continue
    [[ "$f" == treereader* ]] && continue
    [[ "$f" == arborist_*_out.root ]] && continue
    echo "DEBUG: Checking $f for singlephotonana directory..." >&2
    # Check for "singlephotonana" directory.
    if root -l -b -q <<EOF >/dev/null 2>&1
TFile tf("$f");
if (tf.IsZombie()) { gSystem->Exit(1); }
auto* d = tf.GetDirectory("singlephotonana");
gSystem->Exit(d ? 0 : 2);
EOF
    then
      echo "DEBUG: Found ntuple: $f" >&2
      echo "$f"
      return 0
    fi
  done
  
  # If no local ntuple found, try to fetch from SAM
  echo "DEBUG: No local ntuple found, attempting to fetch from SAM..." >&2
  fetch_consumed_file
}

detect_treereader_file() {
  # Prefer newest treereader*.root in cwd.
  local f
  f="$(ls -1t treereader*.root 2>/dev/null | head -n 1 || true)"
  [[ -n "$f" ]] && echo "$f"
}

if [[ $# -eq 0 ]]; then
  echo "DEBUG: Running in auto-detect mode" >&2
  echo "DEBUG: Current directory: $(pwd)" >&2
  echo "DEBUG: Parent directory contents:" >&2
  ls -la .. >&2 || true
  
  TREEREADER_ARTROOT="${TREEREADER_ARTROOT:-$(detect_treereader_file)}"
  # Try to detect ntuple: first check SOURCE_NTUPLE env var, then auto-detect
  if [[ -z "${SOURCE_NTUPLE:-}" ]]; then
    SOURCE_NTUPLE="$(detect_ntuple_file || true)"
  fi
  ARBORIST_FCL="${ARBORIST_FCL:-run_arborist.fcl}"
  if [[ -z "${SOURCE_NTUPLE}" ]]; then
    echo "ERROR: could not auto-detect source ntuple in $(pwd)" >&2
    echo "Set SOURCE_NTUPLE=/path/to/ntuple.root (or run in explicit mode)." >&2
    echo "DEBUG: Environment variables that might help:" >&2
    env | grep -iE 'input|file|sam|source' >&2 || true
    echo "DEBUG: Checking for root files anywhere under working_dir:" >&2
    find "$(dirname "$(pwd)")" -name "*.root" -type f >&2 2>/dev/null || true
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

# Verify the output file was created
if [[ ! -f "${OUTPUT_NTUPLE}" ]]; then
  echo "ERROR: Output ntuple was not created: ${OUTPUT_NTUPLE}" >&2
  exit 7
fi

echo "Output file created: ${OUTPUT_NTUPLE} ($(stat -c%s "${OUTPUT_NTUPLE}") bytes)"

# Copy output to pnfs output directory
# project.py passes --outdir to the job wrapper, we need to find it
OUTPUT_DIR="${OUTDIR:-}"

echo "DEBUG: Looking for output directory..." >&2

# Try various methods to find the output directory
if [[ -z "${OUTPUT_DIR}" ]]; then
  # Check for JOBSUB_SCRIPT_OUTPUT_DIR (set by some job wrappers)
  OUTPUT_DIR="${JOBSUB_SCRIPT_OUTPUT_DIR:-}"
  [[ -n "${OUTPUT_DIR}" ]] && echo "DEBUG: Found OUTDIR via JOBSUB_SCRIPT_OUTPUT_DIR" >&2
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  # Try to get from parent process command line (the job wrapper)
  # This works because endscripts run as children of the main job script
  PARENT_PID="$$"
  while [[ "${PARENT_PID}" != "1" ]]; do
    PARENT_PID="$(ps -o ppid= -p "${PARENT_PID}" 2>/dev/null | tr -d ' ')" || break
    [[ -z "${PARENT_PID}" ]] && break
    if [[ -f "/proc/${PARENT_PID}/cmdline" ]]; then
      CMDLINE="$(tr '\0' ' ' < /proc/${PARENT_PID}/cmdline 2>/dev/null)"
      if [[ "${CMDLINE}" == *"--outdir"* ]]; then
        OUTPUT_DIR="$(echo "${CMDLINE}" | grep -oP '(?<=--outdir\s)\S+' || true)"
        [[ -n "${OUTPUT_DIR}" ]] && echo "DEBUG: Found OUTDIR from parent process cmdline" >&2 && break
      fi
    fi
  done
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  # Try to extract from job script arguments - look for --outdir in command*.txt
  for cmdfile in command*.txt commandStage*.txt; do
    if [[ -f "$cmdfile" ]]; then
      OUTPUT_DIR="$(grep -oP '(?<=--outdir\s)\S+' "$cmdfile" 2>/dev/null | head -1 || true)"
      [[ -n "${OUTPUT_DIR}" ]] && echo "DEBUG: Found OUTDIR in $cmdfile" >&2 && break
    fi
  done
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  # Try to find it in env.txt if it exists
  if [[ -f "env.txt" ]]; then
    OUTPUT_DIR="$(grep -oP '(?<=OUTDIR=)\S+' env.txt 2>/dev/null || true)"
    [[ -n "${OUTPUT_DIR}" ]] && echo "DEBUG: Found OUTDIR in env.txt" >&2
  fi
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  # Look for the job wrapper script in parent directories
  for wrapper in ../*.sh ../../*.sh; do
    if [[ -f "$wrapper" ]]; then
      # Check if this script contains --outdir handling
      OUTPUT_DIR="$(grep -oP '(?<=outdir=")[^"]+' "$wrapper" 2>/dev/null | head -1 || true)"
      if [[ -z "${OUTPUT_DIR}" ]]; then
        OUTPUT_DIR="$(grep -oP "(?<=outdir=')[^']+" "$wrapper" 2>/dev/null | head -1 || true)"
      fi
      [[ -n "${OUTPUT_DIR}" ]] && echo "DEBUG: Found OUTDIR in wrapper script $wrapper" >&2 && break
    fi
  done
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  # Try to infer from the current working directory path pattern
  # Working dir is often like: /pnfs/.../project_name/work
  # Output dir would be: /pnfs/.../project_name/reco
  WORK_DIR="$(pwd)"
  if [[ "${WORK_DIR}" == *"/work"* ]]; then
    INFERRED_OUT="${WORK_DIR/\/work//reco}"
    # Make sure the path looks like pnfs
    if [[ "${INFERRED_OUT}" == /pnfs/* ]] || [[ "${INFERRED_OUT}" == /srv/* ]]; then
      echo "DEBUG: Inferring OUTDIR from work dir pattern: ${INFERRED_OUT}" >&2
      OUTPUT_DIR="${INFERRED_OUT}"
    fi
  fi
fi

if [[ -n "${OUTPUT_DIR}" ]]; then
  echo "Copying output to: ${OUTPUT_DIR}/"
  # Create the output directory if needed (for user-level pnfs)
  ifdh mkdir_p "${OUTPUT_DIR}" >&2 2>/dev/null || true
  
  if ifdh cp "${OUTPUT_NTUPLE}" "${OUTPUT_DIR}/${OUTPUT_NTUPLE}" >&2; then
    echo "Successfully copied ${OUTPUT_NTUPLE} to ${OUTPUT_DIR}/"
  else
    echo "WARNING: Failed to copy output to ${OUTPUT_DIR}/, trying alternative method..." >&2
    # Try with full path
    if ifdh cp "$(pwd)/${OUTPUT_NTUPLE}" "${OUTPUT_DIR}/${OUTPUT_NTUPLE}" >&2; then
      echo "Successfully copied with full path"
    else
      echo "ERROR: Failed to copy output file" >&2
    fi
  fi
else
  echo "WARNING: Could not determine output directory, file remains in $(pwd): ${OUTPUT_NTUPLE}" >&2
  echo "DEBUG: Listing environment variables:" >&2
  env | grep -iE 'out|dir|pnfs|scratch' >&2 | head -20 || true
  echo "DEBUG: Listing files in working directory:" >&2
  ls -la *.txt 2>/dev/null >&2 || true
fi

echo "Done."


