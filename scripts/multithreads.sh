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


maxChains=$((nChainsTot/nChainsPar))
#maxChains=1

echo "--------------------------------"
echo "running 1D transdMT for site "$siteid

# nChainsPar chains at a time (nChainsPar nodes, 1 thread by node)

temp=$output_dir/outfolder/${siteid}

for i in `seq 1 $maxChains`
do
	echo "	running $(($i*$nChainsPar))/$nChainsTot chains  ..."
	
	# start 8 chains in parallel	
	for i in `seq 1 $nChainsPar`
	do
		if [ "$log" = true ] ; then
			../src/transdMT12 $gamma $output_dir/$ipfile $runc ${input_dir}/${siteid}.csv ${temp}/${i} -s $RANDOM -p $samples $prior $min $max -l 1000&

		else
	    		../src/transdMT12 $gamma $ipfile $runc ${input_dir}/${siteid}.csv ${temp}/${i} -s $RANDOM -p $samples $prior $min $max &
		fi
	done


	# wait for all chains to finish
	for job in `jobs -p`
	do
		wait $job
	done	

	# merge generated files
	rm -f ${temp}/sampModels
	for i in `seq 1 $nChainsPar`
	do
		cat ${temp}/${i}sampModels >> ${temp}/${siteid}_sampModels
	done

	rm -f ${temp}/sampResps
	for i in `seq 1 $nChainsPar`
	do
		cat ${temp}/${i}sampResps >> ${temp}/${siteid}_sampResps
	done

	if [ "$log" = true ] ; then
		rm -f ${temp}/log.csv
		for i in `seq 1 $nChainsPar`
	  	do
	    		cat ${temp}/${i}log.csv >> ${temp}/${siteid}_log.csv
	  	done
	fi

	#clear temporary files
	for i in `seq 1 $nChainsPar`
	do
		rm -f ${temp}/${i}sampModels ${temp}/${i}sampResps ${temp}/${i}log.csv 
	done
done

cp ${input_dir}/${siteid}.csv $output_dir/outfolder/csv/
cp ${temp}/${siteid}_sampModels $output_dir/outfolder/samps/
cp ${temp}/${siteid}_sampResps $output_dir/outfolder/samps/
cp ${temp}/${siteid}_log.csv  $output_dir/outfolder/logs/

rm -r ${temp}


echo "inversion done"
echo "output files saved"
echo " "
