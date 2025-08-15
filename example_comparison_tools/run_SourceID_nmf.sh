#!/bin/bash

K="10"

meta_data="/home/joungmin/mst/metadata_feast/metadata_K_"$K".txt"
res_dir="/home/joungmin/mst/SourceID-NMF/results/"

for nmajor in "2" "5"
do
	res_dir_second=$res_dir"K_"$K"_nmajor_"$nmajor"/"
	mkdir $res_dir_second

	for i in {1..1000}; do 
		sim_data="/home/joungmin/mst/simul_data/K_"$K"_nmajor_"$nmajor"/simul_"$i".txt"
		res_data=$res_dir_second"sim_"$i".txt"
		test_timefile=$res_dir"run_time_K_"$K"_nmajor_"$nmajor".txt"

		python SourceID-NMF.py -i $sim_data -n $meta_data -o $res_data
	done
done
