#! /bin/bash

TOP_DIR_NAME="plots"

for DIR0 in "pidvalidcalo" "pidvalidcaloSCE" "pidvalidcali" "pidvalidcaliSCE" "pidvalidcaliSCEtracksFromShower"; do
  for DIR1 in "all_tracks" "muon_candidates" "non_muon_candidates"; do
    for DIR2 in "all" "contained" "uncontained"; do
      FULL_DIR="$TOP_DIR_NAME"/"$DIR0"/"$DIR1"/"$DIR2"/
      echo "making directory" $FULL_DIR
      mkdir -p $FULL_DIR
    done
  done
done
