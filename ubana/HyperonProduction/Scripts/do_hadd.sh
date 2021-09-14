#!/bin/bash

treename=analysisOutputFHC_NuWro_Overlay_Hyperon_cthorpe_rerun_reco2_alloutcomes_NuWro_run1_fhc_Lambda_rerun_reco2_reco2.root
filelist=/uboone/data/users/cthorpe/gridjobs/hyperon/v08_00_00_40/hyperon_tree/process_dataset_cthorpe_rerun_reco2_alloutcomes_NuWro_run1_fhc_Lambda_rerun_reco2_reco2/filesana.list

# Prepare list of files

onelinefilelist=$(cat ${filelist} | tr \\n ' ')

echo
echo $onelinefilelist
echo

hadd -f ${treename} ${onelinefilelist}

