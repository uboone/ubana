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

mv $FCL wrapper_update_slimmed_port_data.fcl

# Generate wrapper.

cat <<EOF > $FCL
#include "wrapper_update_slimmed_port_data.fcl"

physics.producers.nuselMetrics.PortInput:                 "./WCPwork/merge.root"
physics.producers.portedFlash.PortInput:                  "./WCPwork/merge.root"
physics.producers.portedSpacePointsThreshold.PortInput:   "./WCPwork/merge.root"
physics.producers.portedThresholdhit.PortInput:           "./WCPwork/merge.root"

EOF

# Make sure IFDH service is configured in fcl file.

if ! lar --debug-config=/dev/stdout -c $FCL | grep -q IFDH:; then
  cat <<EOF >> $FCL
services.IFDH:
{
}

EOF
fi
