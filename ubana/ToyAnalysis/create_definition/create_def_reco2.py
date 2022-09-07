#!/usr/bin/env python
import os
import subprocess

project="make_events_BNB_TEST"
stage="reco2"
stage_name=stage
usr="cthorpe"

common="file_type mc and ub_project.version prod_v08_00_00_61 and ub_project.name "

defname    = ["reco2"]
fileformat = ["artroot"]

files = [0 for i in range(len(defname))]
for i in range(len(defname)):
    thisdefname = usr + "_" + project+"_"+stage_name+"_"+ defname[i]
    cmd1 = "samweb delete-definition " + thisdefname
    os.system(cmd1)
    cmd = "samweb create-definition " + thisdefname
    cmd = cmd + " '" + common + project + " and file_format " +fileformat[i] + " and ub_project.stage " + stage
    cmd = cmd + " with availability physical'"
    print ""
    print "Creating definition with the command:"
    print cmd
    os.system(cmd)
    count_cmd = "samweb count-definition-files " + thisdefname
    files[i] = subprocess.check_output(count_cmd,shell=True).strip()

print ""
print "Your Definitions:" 
for i in range(len(defname)):
    thisdefname = "cthorpe_" + project+ "_" +stage_name+"_" + defname[i] + " has " + files[i] + " files"
    print thisdefname 

print ""
print "Done!"
