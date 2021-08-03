#! /bin/bash

params=$1
siteid=$2
input_dir=$3
output_dir=$4

#read transdMT parameters file
i=0
while IFS='=' read -r param val
do
   i=$((i+1))
   eval "P$i"=$val
done < "$params"

nChainsPar=$P1
nChainsTot=$P2
runc=$P3
ipfile=$P4
gamma=$P5
samples=$P6
prior=$P7
min=$P8
max=$P9

#echo "$gamma $ipfile $runc ${input_dir}/${siteid}.csv ${output_dir}/${i} -s $RANDOM -p $samples $prior $min $max"

maxChains=$((nChainsTot/nChainsPar))
#maxChains=1

echo "running 1D transdMT for site "$siteid

# nChainsPar chains at a time (nChainsPar nodes, 1 thread by node)

output_dir=outfolder/${siteid}

for i in `seq 1 $maxChains`
do
	echo "	running $(($i*$nChainsPar))/$nChainsTot chains  ..."
	
	# start 8 chains in parallel	
	for i in `seq 1 $nChainsPar`
	do
		../../../src/transdMT5 $gamma $ipfile $runc ${input_dir}/${siteid}.csv ${output_dir}/${i} -s $RANDOM -p $samples $prior $min $max &
	done

	# wait for all chains to finish
	for job in `jobs -p`
	do
	  wait $job
	done	

	# merge generated files
	rm -f ${output_dir}/sampModels
	for i in `seq 1 $nChainsPar`
	do
		cat ${output_dir}/${i}sampModels >> ${output_dir}/${siteid}CsampModels
	done

	rm -f ${output_dir}/sampResps
	for i in `seq 1 $nChainsPar`
	do
		cat ${output_dir}/${i}sampResps >> ${output_dir}/${siteid}CsampResps
	done

	# clear temporary files
	for i in `seq 1 $nChainsPar`
	do
	rm -f ${output_dir}/${i}sampModels ${output_dir}/${i}sampResps 	
	done
done

cp ${input_dir}/${siteid}.csv outfolder/
cp ${output_dir}/${siteid}CsampModels outfolder/
cp ${output_dir}/${siteid}CsampResps outfolder/

rm -r ${output_dir}


echo "inversion done."
echo "output files saved"
echo "--------------------------------"

