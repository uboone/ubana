#!/bin/bash

sed -i 's/reco2.*.root/reco2.root/g' $1
while read p; do 

  NAME=`basename $p`
  echo $NAME >> ${1}out

done < $1


