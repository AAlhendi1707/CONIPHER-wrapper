#!/bin/bash

################################################################################## Input parameters
###################################################################################################
case_id=$1
scriptDir=`pwd`"/src/"
Workdir=""
sampleinfo="${Workdir}/SampleInfo.txt"
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
Rscript ${scriptDir}force_calling_annovar_to_conipher.R $case_id $sampleinfo $inputDir


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
