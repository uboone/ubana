#!/bin/bash

# Purpose: Prestage a dataset
# Inputs: $samdef = the definition you want to prestage 
# Useage: ./prestage_dataset.sh $samdef

nohup samweb prestage-dataset --defname=${1} >& prestage_${1}.log &
