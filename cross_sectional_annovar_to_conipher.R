args <- commandArgs(TRUE)
PatientID <- as.character(args[1])
SNV.dir <-as.character(args[2])
CNV.dir <-as.character(args[3])
Output.dir <-as.character(args[4])

## library required
library(dplyr);library(data.table);library(tidyr)


## library required

## Define some functions
## function to phrase chr of dataframe name into actual dataframe
f <- function(expression) {
  return(eval(parse(text=expression)))
}

## function to phrase chr of dataframe name into actual dataframe

input.file <- paste0(SNV.dir, "/", PatientID, ".deduplicated.combined.mutect2.varscan.pass.hg19_multianno.maf")
dfsnv <- fread(input.file)

## Protein_change, Variant_Classification Variant_Classification
dfsnv <- dfsnv %>% select(Tumor_Sample_Barcode, Hugo_Symbol,Chromosome,Start_Position,End_Position, Reference_Allele,Tumor_Seq_Allele2,t_ref_count, t_alt_count, aaChange, Variant_Classification, Variant_Type)%>%
mutate(Start_Position = as.integer(Start_Position), End_Position = as.integer(End_Position))

## filter for autosmal chrs
dfsnv <- dfsnv %>% filter(Chromosome %in% c(1:22))

## make a snvid
dfsnv$mutation_id <- paste0(dfsnv$Chromosome,":", dfsnv$Start_Position, ":", dfsnv$Reference_Allele)
dfsnv$Chromosome <- as.character(dfsnv$Chromosome)


cnv.file <- paste0(CNV.dir, "/", PatientID, ".segments.txt")
dfcnv <- fread(cnv.file)

dfcnv$CASE_ID <- unique(dfsnv$Tumor_Sample_Barcode)
dfcnv$SAMPLE <- unique(dfsnv$Tumor_Sample_Barcode)
dfcnv$SampleID <- unique(dfsnv$Tumor_Sample_Barcode)

# merge
dfcnv <- dfcnv %>% select(CASE_ID, SAMPLE, SampleID, chr, startpos, endpos, Ploidy, ACF, nMajor, nMinor, nAraw, nBraw)

# rename columns
colnames(dfcnv)[1:6] <- c("CASE_ID", "SAMPLE", "SampleID", "Chromosome","Start_Position","End_Position")
# filter for only Autosomal chromosomes
dfcnv <- dfcnv %>% filter(Chromosome %in% c(1:22))
dfcnv$Chromosome <- as.character(dfcnv$Chromosome)

#Adding CNV minor and major to snv
# combine SNV + CNV and rbind them

cnvi <- as.data.table(dfcnv)
snvi <- as.data.table(dfsnv)
setkey(cnvi, Chromosome, Start_Position, End_Position)
dfxall <- foverlaps(snvi, cnvi)


dfxall <- dfxall %>% filter(nMajor !=0, !is.na(nMajor))

## reformat some columns
dfxall <- dfxall %>% mutate(DEPTH = as.integer(t_ref_count) + as.integer(t_alt_count)) %>%
select(CASE_ID, SAMPLE, Hugo_Symbol,Chromosome,i.Start_Position,Reference_Allele,Tumor_Seq_Allele2,t_ref_count, t_alt_count, DEPTH, aaChange, Variant_Classification, Variant_Type, Ploidy, ACF, nMajor, nMinor, nAraw, nBraw)

# remame some columns
dfxall <- dfxall %>% 
rename(CHR =Chromosome, POS=i.Start_Position, REF =Reference_Allele, ALT = Tumor_Seq_Allele2, REF_COUNT = t_ref_count, VAR_COUNT = t_alt_count,
Protein_change = aaChange, COPY_NUMBER_A = nMajor, COPY_NUMBER_B= nMinor, PLOIDY = Ploidy)


## write the SNV (with local cnv) files  

fwrite(dfxall, paste0(Output.dir, "/", PatientID,".input.tsv"), sep="\t")
