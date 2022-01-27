#! /bin/bash

params=$1
siteid=$2
input_dir=$3
output_dir=$4

nChainsPar=$5
nChainsTot=$6
runc=$7
ipfile=$8
gamma=$9
samples=${10}
prior=${11}
min=${12}
max=${13}
log=${14}

#echo "$gamma $ipfile $runc ${input_dir}/${siteid}.csv ${output_dir}/${i} -s $RANDOM -p $samples $prior $min $max"

maxChains=$((nChainsTot/nChainsPar))
#maxChains=1

echo "--------------------------------"
echo "running 1D transdMT for site "$siteid

# nChainsPar chains at a time (nChainsPar nodes, 1 thread by node)

output_dir=outfolder/${siteid}

for i in `seq 1 $maxChains`
do
	echo "	running $(($i*$nChainsPar))/$nChainsTot chains  ..."
	
	# start 8 chains in parallel	
	for i in `seq 1 $nChainsPar`
	do
		if [ "$log" = true ] ; then
			../../../src/transdMT7 $gamma $ipfile $runc ${input_dir}/${siteid}.csv ${output_dir}/${i} -s $RANDOM -p $samples $prior $min $max -l&

		else
	    		../../../src/transdMT7 $gamma $ipfile $runc ${input_dir}/${siteid}.csv ${output_dir}/${i} -s $RANDOM -p $samples $prior $min $max &
		fi
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
		cat ${output_dir}/${i}sampModels >> ${output_dir}/${siteid}_sampModels
	done

	rm -f ${output_dir}/sampResps
	for i in `seq 1 $nChainsPar`
	do
		cat ${output_dir}/${i}sampResps >> ${output_dir}/${siteid}_sampResps
	done

	if [ "$log" = true ] ; then
		rm -f ${output_dir}/log.csv
		for i in `seq 1 $nChainsPar`
	  	do
	    		cat ${output_dir}/${i}log.csv >> ${output_dir}/${siteid}_log.csv
	  	done
	fi

	#clear temporary files
	for i in `seq 1 $nChainsPar`
	do
		rm -f ${output_dir}/${i}sampModels ${output_dir}/${i}sampResps ${output_dir}/${i}log.csv 
	done
done

cp ${input_dir}/${siteid}.csv outfolder/
cp ${output_dir}/${siteid}_sampModels outfolder/
cp ${output_dir}/${siteid}_sampResps outfolder/
cp ${output_dir}/${siteid}_log.csv outfolder/

rm -r ${output_dir}


echo "inversion done"
echo "output files saved"
echo " "
