#! /bin/bash
#------------------------------------------------------------------
#
# reco2 <initsource> for ExtGen GiBUU overlay samples (BNB): wraps the
# combined-reco-tree fcl and sets the Wire-Cell analysis parameters,
# as def_filetype_kdar_overlay.sh did for the HepMC GiBUU samples:
#   FileType "prodgibuu_GB25_<beam>_nu_overlay_<run period>"
#   (a label in the trees; GB25 = the ExtGen GiBUU tune GB25_05a_00_000).
# For NuMI samples set beam=numi and flag_numi=true below.
#
# The run period comes from the run number of the stage's input file.
#
#------------------------------------------------------------------

beam="bnb"
flag_numi="false"
flag_SaveLeeWeights="false"
BucketTimeSigma='1.3'

if [ x$FCL = x ]; then
  echo "Variable FCL not defined."
  exit 1
fi
if [ ! -f $FCL ]; then
  echo "Fcl file $FCL does not exist."
  exit 1
fi

mv $FCL mix_wrapper.fcl

# run number of the file the next stage reads
next_stage_input=`ls -t1 *.root | egrep -v 'celltree|hist|larlite|larcv|Supplemental|TGraphs' | artroot_filter.py | head -n1`
run_number=`echo $next_stage_input | cut -d '-' -f3`
[[ "$run_number" =~ ^[0-9]+$ ]] && run_number=$((10#$run_number)) || run_number=0
echo "$next_stage_input: run $run_number"

FT_STREAM="unknown"
if   [ "$run_number" -ge  3420 ] && [ "$run_number" -le  8316 ]; 
then 
FT_STREAM="run1"
BucketTimeSigma="2.8"
elif [ "$run_number" -ge  8317 ] && [ "$run_number" -le 13696 ]; then FT_STREAM="run2"
elif [ "$run_number" -ge 13697 ] && [ "$run_number" -le 18960 ]; then FT_STREAM="run3"
elif [ "$run_number" -ge 18961 ] && [ "$run_number" -le 19752 ]; then FT_STREAM="run4a"
elif [ "$run_number" -ge 19753 ] && [ "$run_number" -le 21285 ]; then FT_STREAM="run4b"
elif [ "$run_number" -ge 21286 ] && [ "$run_number" -le 22269 ]; then FT_STREAM="run4c"
elif [ "$run_number" -ge 22270 ] && [ "$run_number" -le 24319 ]; then FT_STREAM="run4d"
elif [ "$run_number" -ge 24320 ] && [ "$run_number" -le 25768 ]; then FT_STREAM="run5"
fi
echo "FileType stream $FT_STREAM"

cat <<EOF2 > $FCL
#include "mix_wrapper.fcl"

physics.analyzers.wcpselection.IsNuMI:               ${flag_numi}
physics.analyzers.wcpweights.IsNuMI:                 ${flag_numi}

physics.analyzers.wcpselection.ssmBDT:               ${flag_numi}

physics.analyzers.wcpselection.get_redk2nu_time:     ${flag_numi}

physics.analyzers.wcpselection.SaveLeeWeights:       ${flag_SaveLeeWeights}
physics.analyzers.wcpweights.SaveLeeWeights:         ${flag_SaveLeeWeights}

physics.analyzers.wcpselection.get_reboone_time:     true
physics.analyzers.wcpselection.TimeBetweenBuckets: 18.936
physics.analyzers.wcpselection.BucketTimeSigma: ${BucketTimeSigma}
physics.analyzers.wcpselection.NBucketsPerBatch: 84
physics.analyzers.wcpselection.NFilledBucketsPerBatch: 81
physics.analyzers.wcpselection.BatchIntensities: [1]

physics.analyzers.wcpselection.MC:              true
physics.analyzers.wcpselection.SaveWeights:     false
physics.analyzers.wcpselection.POT_inputTag:    "generator"

physics.analyzers.wcpselection.FileType: "prodgibuu_GB25_${beam}_nu_overlay_${FT_STREAM}"

physics.ana: [wcpselection]
physics.end_paths: [ stream1, ana ]

services.TFileService.fileName: 	"reco_stage_2_hist.root"
microboone_tfile_metadata.JSONFileName: ["reco_stage_2_hist.root.json"]

EOF2

if ! lar --debug-config=/dev/stdout -c $FCL | grep -q IFDH:; then
  cat <<EOF2 >> $FCL
services.IFDH:
{
}

EOF2
fi
