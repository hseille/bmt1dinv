#! /bin/bash

project=$1
params='inversionParameters.txt'



echo "\n 
############################################
##   Trans-dimensional MCMC 1D Inversion  ##
##        of Magnetotelluric Data         ##
##                                        ##
##       Gerhard Visser / Hoël Seillé     ##
##            CSIRO DEI FSP 2020          ##
############################################"

#read transdMT inversion parameters file
echo ' '
echo '\nInversion parameters:'
while IFS='=' read -r param val
do
   i=$((i+1))
   eval "P$i"=$val
   echo "   $param = $val"
done < "$params"
echo ' '


i=0
while IFS='=' read -r param val
do
   i=$((i+1))
   eval "P$i"=$val
done < "$params"




cd ../projects/$project/transdMT

input_dir=../csv
output_dir=outfolder


#Run the trans-d
for file in $input_dir/*.csv; do #loop1
	siteid=$(basename "$file" .csv)
	mkdir -p ${output_dir}/${siteid}
	bash multithreads.sh $params $siteid $input_dir $output_dir $P1 $P2 $P3 $P4 $P5 $P6 $P7 $P8 $P9 $P10
done #end loop1

echo "goodbye!
"

