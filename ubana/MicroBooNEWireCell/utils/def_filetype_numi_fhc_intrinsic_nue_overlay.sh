#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: This script is intended to update the configuration
#          of the fhicl using the initilization source hook. It
#          makes a wrapper fcl file that overrides certain fcl
#          parameters used by module WireCellAnaTree and
#          WireCellEventWeightTree, specifically:
#        
#          physics.analyzers.wcpselection.FileType: "data_extbnb_run3_G1"
#          physics.analyzers.wcpweights.FileType: "data_extbnb_run3_G1"
#
# Created: Wenqiang Gu, 14-Jan-2021
#         
# Updated: Liang Liu,   03-Jan-2025
#          Unified this script for all runs
#
#------------------------------------------------------------------

# Make sure batch environment variables needed by this script are defined.

if [ x$FCL = x ]; then
  echo "Variable FCL not defined."
  exit 1
fi

# Make sure fcl file $FCL exists.

if [ ! -f $FCL ]; then
  echo "Fcl file $FCL does not exist."
  exit 1
fi

# Rename the existing fcl file $FCL to something else.

mv $FCL mix_wrapper.fcl

# Generate wrapper according to run_number
# get the file that will be used as input for next stage
next_stage_input=`ls -t1 *.root | egrep -v 'celltree|hist|larlite|larcv|Supplemental|TGraphs' | artroot_filter.py | head -n1`
echo $next_stage_input
# get the run number
run_number=`echo $next_stage_input | cut -d '-' -f3`
echo $run_number

FT_STREAM="run1"

STANDARD_OVERLAY="DataOverlayNoTPC"
if [ "$run_number" -ge "3420"  ] && [  "8316" -ge "$run_number"  ];    # in the run1 run number interval
then
        echo "run run1 fhicl"
        FT_STREAM="run1"
STANDARD_OVERLAY="DataOverlayOptical"
elif [ "$run_number" -ge "0008317"  ] && [  "0011048" -ge "$run_number"  ];   # in the run2 run number interval
then
        echo "run run2a fhicl"
        FT_STREAM="run2"
STANDARD_OVERLAY="DataOverlayOptical"
elif [ "$run_number" -ge "0011049"  ] && [  "0013696" -ge "$run_number"  ];   # in the run2 run number interval
then
        echo "run run2b fhicl"
        FT_STREAM="run2"
elif [ "$run_number" -ge "0013697"  ] && [  "0014116" -ge "$run_number"  ];    # in the run3 run number interval
then
        echo "run run3a fhicl"
        FT_STREAM="run3"
elif [ "$run_number" -ge "14117"  ] && [  "0018960" -ge "$run_number"  ];    # in the run3 run number interval
then
        echo "run run3b fhicl"
        FT_STREAM="run3"
elif [ "$run_number" -ge "0018961"  ] && [  "0019752" -ge "$run_number"  ];  # in the run4 run number interval, run4a has four epochs, we may need to split it
then
        echo "run run4a fhicl"
        FT_STREAM="run4a"
elif [ "$run_number" -ge "0019753"  ] && [  "0021285" -ge "$run_number"  ];  # in the run4 run number interval, run4b has four epochs, we may need to split it
then
        echo "run run4b fhicl"
        FT_STREAM="run4b"
elif [ "$run_number" -ge "0021286"  ] && [  "0022269" -ge "$run_number"  ];  # in the run4 run number interval, run4c has four epochs, we may need to split it
then
        echo "run run4c fhicl"
        FT_STREAM="run4c"
elif [ "$run_number" -ge "0022270"  ] && [  "0024319" -ge "$run_number"  ];  # in the run4 run number interval, run4d has four epochs, we may need to split it
then
        echo "run run4d fhicl"
        FT_STREAM="run4d"
elif [ "$run_number" -ge "0024320"  ] && [  "0025768" -ge "$run_number"  ];  # in the run5 run number interval
then
        echo "run run5 fhicl"
        FT_STREAM="run5"
fi

sample_type="intrinsic_nue_overlay"
flag_SaveLeeWeights="false"
flag_numi="true"
beam="numi"
horncur='_fhc'

cat <<EOF > $FCL
#include "mix_wrapper.fcl"

physics.filters.nuselection.AnalysisTools.default.NuMISWTriggerProcName: "${STANDARD_OVERLAY}"
physics.analyzers.wcpselection.IsNuMI:               ${flag_numi}
physics.analyzers.wcpweights.IsNuMI:                 ${flag_numi}

physics.analyzers.wcpselection.ssmBDT:               ${flag_numi}

physics.analyzers.wcpselection.get_redk2nu_time:     ${flag_numi}

physics.analyzers.wcpselection.SaveLeeWeights:       ${flag_SaveLeeWeights}
physics.analyzers.wcpweights.SaveLeeWeights:         ${flag_SaveLeeWeights}


physics.analyzers.wcpselection.FileType: "prodgenie_${beam}${horncur}_${sample_type}_${FT_STREAM}"
physics.analyzers.wcpweights.FileType: "prodgenie_${beam}${horncur}_${sample_type}_${FT_STREAM}"

physics.end_paths: [ stream1, ana ]

services.TFileService.fileName: 	"reco_stage_2_hist.root"
microboone_tfile_metadata.JSONFileName: ["reco_stage_2_hist.root.json"]

EOF

# Make sure IFDH service is configured in fcl file.

if ! lar --debug-config=/dev/stdout -c $FCL | grep -q IFDH:; then
  cat <<EOF >> $FCL
services.IFDH:
{
}

EOF
fi

