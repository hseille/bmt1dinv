#! /bin/bash

project=$1
params='inversionParameters.txt'



echo " # Trans-dimensional MCMC 1D Inversion of Magnetotelluric Data #" 


cd ../projects/$project/transdMT/


#read transdMT inversion parameters file
echo ' '
echo 'Inversion parameters:'
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


input_dir=../projects/$project/csv
output_dir=../projects/$project/transdMT


cd ../../../scripts

# create output folders
mkdir $output_dir/outfolder/samps/ $output_dir/outfolder/plots/ $output_dir/outfolder/logs/ $output_dir/outfolder/csv/ 


#Run the trans-d
for file in $input_dir/*.csv; do 
	siteid=$(basename "$file" .csv)
	mkdir -p ${output_dir}/outfolder/${siteid}

	bash multithreads.sh $params $siteid $input_dir $output_dir $P1 $P2 $P3 $P4 $P5 $P6 $P7 $P8 $P9 $P10

done 
echo "goodbye!"

