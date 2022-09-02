#!/bin/bash

# Purpose: Setup sam processes and create wrapper script
# Inputs: 
# The name of the fhicl you want to make the wrapper for
# The sam definition you want to use as input

fhiclname=run_ToyAnalysis.fcl
def=cthorpe_make_hyperon_events_reco2_numi_fhc_run1_pt1_hyperon_reco2_reco2_first_5_files

#----------Don't edit below this line-----------------------

echo "Setting up sam for fhicl $fhiclname and definition $def"
prjname=${USER}_${def}_`date +%Y%m%d_%H%M%S`
prjurl=`samweb start-project --defname=$def $prjname`
export TMPDIR=/uboone/data/users/$USER/temp
mkdir -p $TMPDIR
cpid=`samweb start-process --appfamily=art --appname=lar --appversion=v03_01_00 $prjurl`
make_sam_wrapper.sh $fhiclname $prjurl $cpid > wrapper.fcl
echo services.FileCatalogMetadata.processID: \"$cpid\" >> wrapper.fcl 
