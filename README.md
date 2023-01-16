# CONIPHER 

## CONIPHER clustering and tree building wrapper

This is a README detailing how to run both mutation clustering and phylogenetic tree building using CONIPHER (COrrecting Noise In PHylogenetic Evaluation and Reconstruction). NOTE: our R package for CONIPHER tree building is available for download [here](https://github.com/McGranahanLab/CONIPHER).

--- 
### Setup

Clone the github repo using the following command from your terminal and enter the directory:
```
git clone https://github.com/McGranahanLab/CONIPHER-wrapper/
cd CONIPHER-wrapper
```

#### Create CONIPHER conda environment
To be able to run CONIPHER clustering and tree building with one wrapper script, follow the steps below. 

To create the conda environment to successfully install and run CONIPHER clustering and tree building, please manually build the conda environment using the instructions below.

On your terminal, ensure you are located in the `CONIPHER-wrapper` directory.

1. Create the conda environment with the following libaries installed
```
conda create -n conipher -c bioconda -c conda-forge pyclone r-base=3.6.1 r-essentials r-tidyverse r-cowplot r-ggpubr r-fst r-biocmanager r-devtools r-seqminer
```

2. Once this has been run, activate the conda environment and start R

```
conda activate conipher
R
```

3. Subsequently install the below list of packages from an R session. NOTE: please do not update related packages during installation when prompted to do so. 

```
# Packages required for CONIPHER clustering

install.packages("mclust")
BiocManager::install("GenomicRanges")
BiocManager::install("Rsamtools")
install.packages("gplots")
install.packages("gdata")
install.packages("future")
install.packages("optparse")
install.packages("bootstrap")
BiocManager::install("copynumber")
devtools::install_version("sequenza", version = "2.1.2")
install.packages("coin")
install.packages("wordcloud")


# CONIPHER treebuilding R package
devtools::install_github("McGranahanLab/CONIPHER")
```

4. Once all of these have been installed quit R and deactivate the conda environment

```
q()
conda deactivate
```

--- 

### Quickstart
#### Running CONIPHER clustering + tree building

We provide a wrapper bash script to run CONIPHER clustering and tree building end-to-end. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then run the following command:

```
sh 0_runningClusteringTreeBuilding.sh
```


#### Running CONIPHER tree building

We additionally provide a wrapper script to run CONIPHER tree building by itself. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then run the following command:

```
sh 0_runningTreeBuilding.sh
```

#### Input parameters

option(c("-t", "--input_tsv"), type="character", default=NULL, 
              help="A mutation table in long format (mutations x tumour regions)", metavar="character")

**This is essentail** to run CONIPHER. Columns:
SAMPLE	CHR	POS	REF	ALT	REF_COUNT	VAR_COUNT	DEPTH	COPY_NUMBER_A	COPY_NUMBER_B	PLOIDY	ACF

  
option(c("--input_seg_tsv"), type="character", default=NULL, 
              help="A copy number segment .tsv used for plotting only", metavar="character")
**This optional**  as `seg.mat.copy` can replace it if input_tsv get the copy number info. Columns:
COPY_NUMBER_A	COPY_NUMBER_B	PLOIDY	ACF

option(c("--subclonal_copy_correction"), type="logical", default="TRUE", 
              help="should subclonal copy number correction be used", metavar="character")
**This is not required** as no use for it in the script, as only mode  to work without it “FALSE mode”. By defult not needed as seg.mat.copy will automatically convert into phylo seg.mat.phylo in FALSE mode. 
  
option(c("--only_truncal_subclonal_copy_correction"), type="logical", default="TRUE", 
help="should only truncal subclonal copy number correction be used",metavar="character")
**This is required** and it is TRUE by default, It is a sort of fix for CNA correction when it didn't make the mutation clonal

option(c("--nProcs"), type = "integer", default = 1,
              help="Number of cores allocated to run script in parallel", metavar="character")
**This can be increased** to speed up the pyclone run

##### Rest of options are good
option(c("--clean_clusters"), type="logical", default=TRUE, 
              help="should clusters be cleaned and merged?", metavar="character")
              
option(c("--clonal_cutOff"), type="numeric", default=0.9, 
              help="lower threshold CCF to be considered clonal", metavar="character")
              
option(c("--propClonal_threshold"), type="numeric", default=0.25, 
              help="proportion of cluster that needs to be considered clonal to merge", metavar="character")
              
option(c("--fix_absentCCFs"), type="logical", default=TRUE, 
              help="should CCF of absent mutations be set to zero?", metavar="character")
              
option(c("--driver_filter"), type="character", default="1A,1,2A", 
              help="what filter to use for drivers", metavar="character")
              
option(c("--burn_in"), type="integer", default=1000, 
              help="burn-in for DP clustering", metavar="character")

