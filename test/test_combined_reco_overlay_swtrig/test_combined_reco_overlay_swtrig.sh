#! /bin/bash

# Loop over all combined reco overlay ntuple fcls.

ok=0   # 0 = success

while read fcl
do
  fclname=`basename $fcl`
  echo
  echo "Testing fcl file $fclname"

  # Parse this fcl file.

  fclout=`basename ${fcl}`.out
  larout=`basename ${fcl}`.lar.out
  larerr=`basename ${fcl}`.lar.err
  lar -c $fcl --debug-config $fclout > $larout 2> $larerr

  # Any non-zero exit status exit immediately (fail).

  stat=$?
  if [ $stat -ne 0 ]; then
    echo "Error parsing ${fcl}."
    exit $stat
  fi

  # Calculate expected process name.

  pname=''
  if [[ $fcl == *run1* ]]; then
    if [[ $fcl == *numi* ]]; then
      echo "Run 1 numi"
      trigsim=standard_overlay_optical_numi_uboone_updated.fcl
    else
      echo "Run 1 bnb"
      trigsim=standard_overlay_optical_uboone.fcl
    fi
  elif [[ $fcl == *run3* ]]; then
    if [[ $fcl == *numi* ]]; then
      echo "Run 3 numi"
      trigsim=standard_overlay_notpc_numi_uboone_updated.fcl
    else
      echo "Run 3 bnb"
      trigsim=standard_overlay_notpc_uboone.fcl
    fi
  elif [[ $fcl == *run4* ]]; then
    if [[ $fcl == *numi* ]]; then
      echo "Run 4 numi"
      trigsim=standard_overlay_notpc_numi_uboone_updated.fcl
    else
      echo "Run 4 bnb"
      trigsim=standard_overlay_notpc_uboone.fcl
    fi
  else
    echo "Unknown epoch."
    exit 1
  fi
  echo "Trigsim fcl file = $trigsim"
  pname=`fhicl-dump $trigsim | grep process_name | awk '{print $2}' | tr -d '"'`
  echo "Expected swtrig process name = $pname"

  # Check swtrig process name.

  pname2=`grep NuMISWTriggerProcName $fclout | awk '{print $2}' | tr -d '"' | sort -u`
  echo "Actual swtrig process name   = $pname2"
  if [ "$pname" != "$pname2" ]; then
    echo "Process name mismatch."
    ok=1
  else
    echo "Process name OK."
  fi

done < <( find $MRB_BUILDDIR/ubana/job -name run_combinedrecotree_run\*overlay\*.fcl -print )


# Done

echo
if [ $ok -eq 0 ]; then
  echo "Overall result: pass"
else
  echo "Overall result: fail"
fi

exit $ok
