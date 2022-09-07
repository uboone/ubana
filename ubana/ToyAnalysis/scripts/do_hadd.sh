#!/bin/bash

# Purpose: Combine lots of ntuples into a single file
# Useage: source do_hadd.sh 
# Set value of treename first, include the sam definition.

treename=analysisOutput_
filelist=filesana.list

#Prepare list of files
onelinefilelist=$(cat ${filelist} | grep 'OutputTrees.root' | tr \\n ' ')

echo
echo $onelinefilelist
echo

hadd -f ${treename} ${onelinefilelist}
