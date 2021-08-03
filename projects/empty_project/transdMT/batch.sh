#! /bin/bash

params=$1

#read transdMT parameters file
echo 'Inversion parameters:'
while IFS='=' read -r param val
do
   echo "   $param = $val"
done < "$params"
echo ' '


input_dir=../csv
output_dir=outfolder


#Run the trans-d
for file in $input_dir/*.csv; do #loop1
	siteid=$(basename "$file" .csv)
	mkdir -p ${output_dir}/${siteid}
	bash multithreads.sh $params $siteid $input_dir $output_dir
done #end loop1



