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
next_stage_input=`ls -t1 *.root | egrep -v 'celltree|hist|larlite|larcv|Supplemental|TGraphs' | head -n1`
echo $next_stage_input
# get the run number
run_number=`echo $next_stage_input | cut -d '-' -f3`
echo $run_number

FT_STREAM="run1"
BucketTimeSigma='1.3'
if [ "$run_number" -ge "3420"  ] && [  "8316" -ge "$run_number"  ];    # in the run1 run number interval
then
        echo "run run1 fhicl"
        FT_STREAM="run1"
	BucketTimeSigma="2.5"
elif [ "$run_number" -ge "0008317"  ] && [  "0011048" -ge "$run_number"  ];   # in the run2 run number interval
then
        echo "run run2a fhicl"
        FT_STREAM="run2"
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



sample_type='intrinsic_nue'
flag_SaveLeeWeights="true"
flag_numi="false"
flag_redk2nu="false"
beam="bnb"
horncur=""



detvar_type="cv"
if ls *reg4_SCE*.fcl 1> /dev/null 2>&1; then
	detvar_type="SCE"
elif ls *reg4_recombination*.fcl 1> /dev/null 2>&1; then
        detvar_type="recombination"
elif ls *wiremod_ScaleX*.fcl 1> /dev/null 2>&1; then
        detvar_type="WMX"
elif ls *wiremod_ScaleYZ*.fcl 1> /dev/null 2>&1; then
        detvar_type="WMYZ"
elif ls *wiremod_ScaleAngleYZ*.fcl 1> /dev/null 2>&1; then
        detvar_type="WMThetaYZ"
elif ls *wiremod_ScaleAngleXZ*.fcl 1> /dev/null 2>&1; then
        detvar_type="WMThetaXZ"
elif ls *reg4_LY_reyliegh*.fcl 1> /dev/null 2>&1; then
        detvar_type="LY_reyliegh"
elif ls *reg4_LY_suppression75*.fcl 1> /dev/null 2>&1; then
        detvar_type="LY_down"
elif ls *reg4_LY_attenuation8m*.fcl 1> /dev/null 2>&1; then
        detvar_type="LY_attenuation"
else
	flag_redk2nu="true"
	echo "WARNING: Unable to establish DetVar type."
        echo "Setting this to tbe the CV sample."
        echo "Please check fhicls to ensure this is correct."
fi


cat <<EOF > $FCL
#include "mix_wrapper.fcl"

physics.analyzers.wcpselection.IsNuMI:               ${flag_numi}
physics.analyzers.wcpweights.IsNuMI:                 ${flag_numi}

physics.analyzers.wcpselection.ssmBDT:               ${flag_numi}

physics.analyzers.wcpselection.get_reboone_time:     ${flag_reboone}
physics.analyzers.wcpselection.TimeBetweenBuckets: 18.831
physics.analyzers.wcpselection.BucketTimeSigma: ${BucketTimeSigma}
physics.analyzers.wcpselection.NBucketsPerBatch: 84
physics.analyzers.wcpselection.NFilledBucketsPerBatch: 81
physics.analyzers.wcpselection.BatchIntensities: [1]

physics.analyzers.wcpselection.SaveLeeWeights:       ${flag_SaveLeeWeights}

physics.analyzers.wcpselection.FileType: "prodgenie_${beam}${horncur}_${sample_type}_detvar_${detvar_type}_${FT_STREAM}"
physics.analyzers.wcpweights.FileType: "prodgenie_${beam}${horncur}_${sample_type}_detvar_${detvar_type}_${FT_STREAM}"

physics.ana: [wcpselection]
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

