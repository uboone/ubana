#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: This script is intended to update the configuration
#          of the fhicl using the initilization source hook. It
#          makes a wrapper fcl file that overrides certain fcl
#          parameters.
#
# Created: Wenqiang Gu, 14-Jan-2021
#
#------------------------------------------------------------------

# Make sure batch environment variables needed by this script are defined.

if [ x$FCL = x ]; then
  echo "Variable FCL not defined."
  exit 1
fi

# Make sure fcl file $FCL exists.

if [ ! -f $FCL ]; then
  echo "Fcl file $FCL does not exist."
  exit 1
fi

# Rename the existing fcl file $FCL to something else.

mv $FCL wrapper_update_wcpf_port.fcl

flag_numi="false"
if ls Physics*numi*.root 1> /dev/null 2>&1 || ls Beam*numi*.root 1> /dev/null 2>&1 ; then
        flag_numi="true"
        echo Running NuMI fhcl
elif ls *numi*.fcl 1> /dev/null 2>&1; then
        flag_numi="true"
        echo Running NuMI fhcl
else
	echo Running BNB fhcl
fi



# Generate wrapper.

cat <<EOF > $FCL
#include "wrapper_update_wcpf_port.fcl"

physics.producers.wirecellPF.ssmBDT:                      ${flag_numi}
physics.producers.wirecellPF.PFInput_prefix:              "WCPwork/nue"

EOF

# Make sure IFDH service is configured in fcl file.

if ! lar --debug-config=/dev/stdout -c $FCL | grep -q IFDH:; then
  cat <<EOF >> $FCL
services.IFDH:
{
}

EOF
fi
