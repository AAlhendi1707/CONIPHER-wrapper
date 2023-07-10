#!/bin/bash

################################################################################## Input parameters
###################################################################################################
case_id=$1
scriptDir="path/to/src/"
Workdir="/path/to/CONIPHER/${case_id}/"

# prepare
#SNV is maf single sample each
#CNV is ASCAT raw segments
SNVDir="path/to/SNV/combined/Splitted_bySample/"
CNVDir="path/to/CNV/ASCAT/Combined/"

# CONIPHER
inputDir="${Workdir}/input/"
outDir="${Workdir}/results/"


############################################################### Running clustering and treebuilding
###################################################################################################


clusteringDir=${outDir}"/Clustering/"
treeDir=${outDir}"/TreeBuilding/"

mkdir -p ${clusteringDir}
mkdir -p ${treeDir}

## prepare forcecalling annovar table to conipher tsv file

## prepare forcecalling annovar table to conipher tsv file
Rscript ${scriptDir}cross_sectional_annovar_to_conipher.R $case_id $SNVDir $CNVDir $inputDir


# run clustering
Rscript ${scriptDir}run_clustering.R \
--case_id ${case_id} \
--script_dir ${scriptDir} \
--input_tsv ${inputDir}${case_id}".input.tsv" \
--working_dir ${clusteringDir} \
--nProcs 20

# run treebuild
Rscript ${scriptDir}run_treebuilding.R \
--input_tsv ${clusteringDir}${case_id}".SCoutput.CLEAN.tsv" \
--out_dir ${treeDir} \
--script_dir ${scriptDir} \
--prefix ${case_id}

module load R
# plot clonemap and write the shannon index
Rscript ${scriptDir}clonemap.diversity.R $treeDir $case_id
