#!/bin/bash

K="10"
metadata="/home/joungmin/mst/metadata_st/metadata_K_"$K".txt"

for nmajor in "2" "5"
do
	res_dir="/home/joungmin/mst/sourcetracker2/results/K_"$K"_nmajor_"$nmajor"/"
	mkdir $res_dir
	test_timefile=$res_dir"run_time_K_"$K"_nmajor_"$nmajor".txt"

	for i in {1..1000}; do
		sim_data="/home/joungmin/mst/simul_data/K_"$K"_nmajor_"$nmajor"/modified/simul_"$i".txt"
		st_res_dir=$res_dir"simul_"$i"/"
		sourcetracker2 gibbs -i $sim_data -m $metadata -o $st_res_dir
	done
done
