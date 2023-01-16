args <- commandArgs(TRUE)
PatientID <- as.character(args[1])
SampleInfo <-as.character(args[2])
SNV.file <-as.character(args[3])
Work.dir <-as.character(args[4])

## library required
library(dplyr);library(data.table);library(tidyr)

suppressPackageStartupMessages(require(maftools))

## library required

## Define some functions
## function to phrase chr of dataframe name into actual dataframe
f <- function(expression) {
  return(eval(parse(text=expression)))
}


df <- fread(SNV.file)
SampleInfo <- fread(SampleInfo)

# sort as in vcf info columns
SampleInfo <- SampleInfo[order(SampleInfo$Tumor_Sample_Barcode),]

### convert into barcode
Samples <- SampleInfo$SampleID
Germline <- "Germline"

## change column namses for samples
colnumbers <- length(colnames(df))
colnames(df)[125:colnumbers] <- Samples

## drop germline column '71005172'
df <- df %>% select(-Germline)

## melt to genrate barcode
ids <- colnames(df)[1:124]
df2 <- melt(df, id=ids, variable.name = "Barcode", value.name = "VCF_INFO")

filename.barcoded <- paste0(PatientID, ".somatic_oncefiltered.mutect2.pass.hg19_multianno.withBarcode.txt")
fwrite(df2, filename.barcoded, sep="\t")

## convert to maf

var.maf <- annovarToMaf(annovar = filename.barcoded, Center = NULL, refBuild = "hg19", 
tsbCol = "Barcode", table = "refGene", sep = "\t")
maf.filename <- paste0(PatientID, ".somatic_oncefiltered.mutect2.pass.hg19_multianno.withBarcode.maf")
fwrite(var.maf, maf.filename, sep="\t")


## step 1 - processing combined SNV maf
dfsnv <- var.maf

## format columns - this depent on somatic caller, this works for mutect2
dfsnv <- dfsnv %>% separate(VCF_INFO, c("Tumor_GT", "Tumor_AD", "Tumor_AF", "Tumor_DP", "Tumor_F1R2", "Tumor_F2R1", "Tumor_SB"), sep=":")

## Protein_change, Variant_Classification Variant_Classification
dfsnv <- dfsnv %>% separate(Tumor_AD, c("t_ref_count", "t_alt_count"), sep=",") %>% 
select(Tumor_Sample_Barcode, Hugo_Symbol,Chromosome,Start_Position,End_Position, Reference_Allele,Tumor_Seq_Allele2,t_ref_count, t_alt_count, aaChange, Variant_Classification, Variant_Type)%>%
mutate(Start_Position = as.integer(Start_Position), End_Position = as.integer(End_Position))

## filter for autosmal chrs
dfsnv <- dfsnv %>% filter(Chromosome %in% c(1:22))

## make a snvid
dfsnv$mutation_id <- paste0(dfsnv$Chromosome,":", dfsnv$Start_Position, ":", dfsnv$Reference_Allele)
dfsnv$Chromosome <- as.character(dfsnv$Chromosome)
## processing cnv
SampleInfo <- SampleInfo[!SampleInfo$SampleID == "Germline", ]


dfcnv <- data.frame()

for ( x in unique(SampleInfo$Tumor_Sample_Barcode)){
    cnv.file <- paste0(Work.dir, "/", x, ".segments.txt")
    cnv.temp <- fread(cnv.file)
    cnv.temp$sample <- x
    dfcnv <- rbind(dfcnv, cnv.temp)
    }

## with dfcnv with sampleinfo
SampleInfo$Tumor_Sample_Barcode <- as.character(SampleInfo$Tumor_Sample_Barcode)
dfcnv$Barcode <- as.character(dfcnv$sample)
dfcnv$sample <- NULL

dfcnv <- merge(dfcnv, SampleInfo, by.x="Barcode", by.y="Tumor_Sample_Barcode")

dfcnv <- dfcnv %>% select(CASE_ID, SAMPLE, SampleID, chr, startpos, endpos, Ploidy, ACF, nMajor, nMinor, nAraw, nBraw)
# rename columns
colnames(dfcnv)[1:6] <- c("CASE_ID", "SAMPLE", "SampleID", "Chromosome","Start_Position","End_Position")
# filter for only Autosomal chromosomes
dfcnv <- dfcnv %>% filter(Chromosome %in% c(1:22))
dfcnv$Chromosome <- as.character(dfcnv$Chromosome)

SampleNames <- unique(SampleInfo$SampleID)

## step 4 - Split maf and CNV by sample id
for ( x in SampleNames){
    assign (paste0("SNV.",x) , dfsnv[dfsnv$Tumor_Sample_Barcode == x, ] )
    assign (paste0("CNV.",x) , dfcnv[dfcnv$SampleID == x, ] )
}

## step 6 - Adding CNV minor and major to snv
# loop throught sample to combine SNV + CNV and rbind them
dfxall = data.frame()
for ( i in SampleNames){
    cnvi <- f(paste0("CNV.", i))
    snvi <- f(paste0("SNV.", i))
    setkey(cnvi, Chromosome, Start_Position, End_Position)
    df.temp <- foverlaps(snvi, cnvi)
    dfxall <- rbind(dfxall,df.temp)         
}

dfxall <- dfxall %>% filter(nMajor !=0, !is.na(nMajor))

## reformat some columns
dfxall <- dfxall %>% mutate(DEPTH = as.integer(t_ref_count) + as.integer(t_alt_count)) %>%
select(CASE_ID, SAMPLE, Hugo_Symbol,Chromosome,i.Start_Position,Reference_Allele,Tumor_Seq_Allele2,t_ref_count, t_alt_count, DEPTH, aaChange, Variant_Classification, Variant_Type, Ploidy, ACF, nMajor, nMinor, nAraw, nBraw)

# remame some columns
dfxall <- dfxall %>% 
rename(CHR =Chromosome, POS=i.Start_Position, REF =Reference_Allele, ALT = Tumor_Seq_Allele2, REF_COUNT = t_ref_count, VAR_COUNT = t_alt_count,
Protein_change = aaChange, COPY_NUMBER_A = nMajor, COPY_NUMBER_B= nMinor, PLOIDY = Ploidy)


## step 7 - write the SNV (with local cnv) files  

fwrite(dfxall, paste0(Work.dir, "/", PatientID,".input.tsv"), sep="\t")
