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
      pname=DataOverlayOpticalNuMI
    else
      echo "Run 1 bnb"
      pname=DataOverlayOptical
    fi
  elif [[ $fcl == *run3* ]]; then
    if [[ $fcl == *numi* ]]; then
      echo "Run 3 numi"
      pname=DataOverlayNoTPCNuMI
    else
      echo "Run 3 bnb"
      pname=DataOverlayNoTPC
    fi
  elif [[ $fcl == *run4* ]]; then
    if [[ $fcl == *numi* ]]; then
      echo "Run 4 numi"
      pname=DataOverlayNoTPCNuMI
    else
      echo "Run 4 bnb"
      pname=DataOverlayNoTPC
    fi
  else
    echo "Unknown epoch."
    exit 1
  fi
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
