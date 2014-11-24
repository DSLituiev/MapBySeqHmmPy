#!/bin/bash
species=( 'Arabidopsis_lyrata' \
'Aethionema_arabicum' \
'Capsella_rubella' \
'Brapa_197' \
'Brassica_oleracea' \
'Leavenworthia_alabamica' \
'Sisymbrium_irio' \
'Solanum_lycopersicum' \
'Vitis_vinifera' \
)

for ((i = 0; i < ${#species[@]}; i++))
do    
     python3 ./run_loc_blast_sql.py $1 -o "${species[i]}"

    echo "              blasting against " ${species[i]} ;
    echo "======================================================";

done

echo "finished"

