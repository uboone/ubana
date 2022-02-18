#!/bin/bash

while read file; do
  lar -c standard_overlay_optical_numi_uboone_updated.fcl -s ${file} 
done < /uboone/app/users/ohanabr/muDAR/gridsub_Mar01_21/muDAR_reco2_list.txt