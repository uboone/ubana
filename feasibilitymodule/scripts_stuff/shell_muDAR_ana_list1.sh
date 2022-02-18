#!/bin/bash

while read file; do
  lar -c ../run_feasibility.fcl -s ${file} 
done < /uboone/app/users/ohanabr/muDAR/ubana_v08_00_00_52/srcs/ubana/feasibilitymodule/scripts_stuff/reco2_opt_updated_list1.txt