offset=( 0.0 -1.0 -2.0 -3.0 )
sigma=( 0.5 1.0 2.0 )
chanceOfFlip=( 0.0 0.25 0.5 0.75 1.0 )

for i in ${!offset[@]}; do
  for j in ${!sigma[@]}; do
    for k in ${!chanceOfFlip[@]}; do 
      echo "${offset[i]} ${sigma[j]} ${chanceOfFlip[k]}"
      ./pidaStudies prodgenie_bnb_nu_uboone_mcc8.6_reco2.list ${offset[i]} ${sigma[j]} ${chanceOfFlip[k]}
    done
  done
done
