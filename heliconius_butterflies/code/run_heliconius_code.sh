#!/usr/bin/env bash

## How to run butterfly analysis in "Leveraging shared ancestral variation to detect local introgression" Lesly Lopez Fang, Diego Ortega Del-Vecchyo, Emily Jane McTavish and Emilia Huerta-Sanchez.


virtualenv -p python3 venv-butterflies

source venv-butterflies/bin/activate
pip install -r requirements.txt

mkdir data
mkdir results
mkdir results/graphs

git clone https://github.com/LeslyLopezFang/Dplus/

#Download data set from Martin, Simon H. et al. (2013), Data from: Genome-wide evidence for speciation with gene flow in Heliconius butterflies, Dryad, Dataset, https://doi.org/10.5061/dryad.dk712
wget -O data/set31.Zupdated.union.geno.part_1_of_2.gz https://datadryad.org/stash/downloads/file_stream/63229
wget -O data/set31.Zupdated.union.geno.part_2_of_2.gz https://datadryad.org/stash/downloads/file_stream/63228

cat data/set31.Zupdated.union.geno.part_1_of_2.gz data/set31.Zupdated.union.geno.part_2_of_2.gz > data/heliconius_geno.csv.gz

#Run statistics analysis on heliconius genome data
cd code

infile="../data/heliconius_geno.csv.gz"
outfile="../results/helico_pop-pi_stats.csv"

population_string="pop1[ag108,ag572,ag112,ag569];pop2[am216,am160,am48,am293];pop3[tiP86,tiP313,tiP84,tiP57];outgroup[hec273,eth67,ser202,par371]"

python3 heliconius_butterfly_sliding_windows.py -i ${infile} -o ${outfile} --d_statistic --dplus_statistic --fd_statistic --fdm_statistic --df_statistic -ws 5000 -p ${population_string}

#Graph the results
results_infile="../results/helico_pop-pi_stats.csv"
outfile_path="../results/graphs"

python3 graph_butterfly_figures.py -i ${results_infile} -o ${graphs_outfile_path}



