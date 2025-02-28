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
isRun3="false"
ShiftOffset="0"
if [ "$run_number" -ge "3420"  ] && [  "8316" -ge "$run_number"  ];    # in the run1 run number interval
then
        echo "run run1 fhicl"
        FT_STREAM="run1"
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
	isRun3="true"
        ShiftOffset="387.8"
elif [ "$run_number" -ge "14117"  ] && [  "0018960" -ge "$run_number"  ];    # in the run3 run number interval
then
        echo "run run3b fhicl"
        FT_STREAM="run3"
	isRun3="true"
        ShiftOffset="387.8"
elif [ "$run_number" -ge "0018961"  ] && [  "0019752" -ge "$run_number"  ];  # in the run4 run number interval, run4a has four epochs, we may need to split it
then
        echo "run run4a fhicl"
        FT_STREAM="run4a"
	ShiftOffset="118.3"
elif [ "$run_number" -ge "0019753"  ] && [  "0021285" -ge "$run_number"  ];  # in the run4 run number interval, run4b has four epochs, we may need to split it
then
        echo "run run4b fhicl"
        FT_STREAM="run4b"
	ShiftOffset="118.3"
elif [ "$run_number" -ge "0021286"  ] && [  "0022269" -ge "$run_number"  ];  # in the run4 run number interval, run4c has four epochs, we may need to split it
then
        echo "run run4c fhicl"
        FT_STREAM="run4c"
	ShiftOffset="118.3"
elif [ "$run_number" -ge "0022270"  ] && [  "0024319" -ge "$run_number"  ];  # in the run4 run number interval, run4d has four epochs, we may need to split it
then
        echo "run run4d fhicl"
        FT_STREAM="run4d"
	ShiftOffset="118.3"
elif [ "$run_number" -ge "0024320"  ] && [  "0025768" -ge "$run_number"  ];  # in the run5 run number interval
then
        echo "run run5 fhicl"
        FT_STREAM="run5"
	ShiftOffset="118.3"
fi

ccnd1_a="0.529594"
ccnd1_b="7.13804"
ccnd2_a="0.068752"
ccnd2_b="2.32023"
ccnd3_a="0.4697"
ccnd3_b="0.004233"
ccnd3_c="0.000001006"
ccnd3_d="-0.195"
ccnd4_a="0"
ccnd4_b="0"

flag_numi="false"
beam="bnb"
if ls Physics*numi*.root 1> /dev/null 2>&1 || ls Beam*numi*.root 1> /dev/null 2>&1 ; then
	flag_numi="true"
        beam="numi"
        flag_numi="true"
        ccnd1_a="0.4343"
        ccnd1_b="6.2884"
        ccnd2_a="0.0637"
        ccnd2_b="1.489"
        ccnd3_a="0"
        ccnd3_b="0"
        ccnd3_c="0"
        ccnd3_d="0"
        ccnd4_a="0.0125"
        ccnd4_b="2.3152"
        isRun3="false"
        ShiftOffset="0"
fi

flag_ext="data"
if ls Physics*ext*.root 1> /dev/null 2>&1 || ls Beam*ext*.root 1> /dev/null 2>&1 ; then
	flag_ext="ext"
fi


cat <<EOF > $FCL
#include "mix_wrapper.fcl"

physics.analyzers.wcpselection.ssmBDT:               ${flag_numi}

physics.analyzers.wcpselection.ccnd1_a: ${ccnd1_a} 
physics.analyzers.wcpselection.ccnd1_b: ${ccnd1_b}
physics.analyzers.wcpselection.ccnd2_a: ${ccnd2_a} 
physics.analyzers.wcpselection.ccnd2_b: ${ccnd2_b} 
physics.analyzers.wcpselection.ccnd3_a: ${ccnd3_a}
physics.analyzers.wcpselection.ccnd3_b: ${ccnd3_b}
physics.analyzers.wcpselection.ccnd3_c: ${ccnd3_c}
physics.analyzers.wcpselection.ccnd3_d: ${ccnd3_d}
physics.analyzers.wcpselection.ccnd4_a: ${ccnd4_a} 
physics.analyzers.wcpselection.ccnd4_b: ${ccnd4_b} 

physics.analyzers.wcpselection.ShiftOffset: ${ShiftOffset}
physics.analyzers.wcpselection.isRun3: ${isRun3}

physics.analyzers.wcpselection.FileType: "${flag_ext}_${beam}_${FT_STREAM}"
physics.analyzers.wcpweights.FileType: "${flag_ext}_${beam}_${FT_STREAM}"

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

