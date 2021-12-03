## ---------------------------
##
## Script name: process-protein-groups.R 
##
## Purpose of script: process proteinGroups.txt file from MaxQuant to include in Expression Atlas
##
## Author: Dr. Andrew Jarnuczak
##
## Date Created: 2019-07-01
##
## Copyright (c) Andrew Jarnuczak, 2019
## Email: jarnuczak@ebi.ac.uk
##
## ---------------------------
##
## Notes: Steps performed in the script
## 1. clean up proteinGroups file (remove CONTAMINANTS and REVERSE)
## 2. for MS1 based quantification results (i.e. label free or SILAC) normalise iBAQ intensities to ppb
## 3. perform proteinGroup to Gene ID mapping based on mapping provided by UniProt 
## 4. perform any other necessary clean-up steps (this will depend on the dataset, for example if there are any "technical" assays included in the original dataset that should be removed before data is included in Expression Atlas) 
##   
##
## ---------------------------

library(matrixStats)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)

# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")

library(mygene)

setwd("/Users/ananth/Documents/MSstats_scripts_for_Andrew")

#In case that the library mygene was not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("mygene")

if ( !exists("EXP_TYPE")) warning("Please specify experiment type variable: EXP_TYE")
## Specify experiment type, i.e. what type of quantification is used, MS1-based or MS2-based. 
## As a rule of thumb, label free and SILAC is MS1-based, and iTRAQ and TMT is MS2-based.
 EXP_TYPE <- "MS2-quant"
#EXP_TYPE <- "MS1-quant"



##############
# define FOT normalisation function
# FOT stands for Fraction Of Total. In this normalisation method each protein iBAQ intensity value is scaled to the total amount of signal in a given MS run (column) and transformed to parts per billion (ppb)
fot.normalise <- function(x){
  data.sum <-   apply(x, 2, function(y){sum(y, na.rm=TRUE)})
  # barplot((data.sum), log = "y")
  #
  ##### do ppm normalisation
  x.mat <- as.matrix(x)
  x.mat.ppb <- apply(x.mat, 2, function(i) i/sum(i, na.rm = T) * 1000000000 )
  x.mat.ppb <- as.data.frame(x.mat.ppb)
  colnames(x.mat.ppb) <- paste("ppb.", colnames(x.mat.ppb), sep = "")
  return(x.mat.ppb)
}
##############


##### read the protein expression matrix produced by MaxQuant
tmp  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

## clean up 
tmp <- tmp[ tmp[, "Reverse"] != "+", ]
tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]
## consider those protein groups with more than 1 peptides matched
tmp <- tmp[ tmp[, "Peptides"] > 1, ]

## 
# if experiment is label free use iBAQ quant:
if(EXP_TYPE == "MS1-quant"){
  message("Collecting iBAQ quantification")
  tmp <- tmp[ ,c(2, grep("iBAQ.", colnames(tmp))) ]
}


# if experiment is TMT, use reporter intensities, if experiment iTRAQ use Intensities. Note, MaxQuant might change how it reports iTRAQ and TMT intensities.
if(EXP_TYPE == "MS2-quant"){
  message("Collecting MS2 intensities")
  if( any(grepl("Reporter.intensity.corrected", colnames(tmp))) ){
    
    # code modified by Ananth:
    # Filter to consider only protein groups for which the at least 50 % of the sample replicates 
    # have non-zero intensity values
    
    ##### Read sample replicate information
    ##### Note: User defined annotation file
    sample_replicates <- read.table("sample_replicates.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
    
    for(i in 1:nrow(sample_replicates))
    {
      sample_replicates[i,] <- gsub("^","Reporter.intensity.corrected.", sample_replicates[i,], perl=TRUE)
      sample_replicates[i,] <- gsub(" ",".", sample_replicates[i,], perl=TRUE)
      sample_replicates[i,] <- gsub("-",".", sample_replicates[i,], perl=TRUE)
      
      subdata <- tmp[, c(unname(unlist(sample_replicates[i,-1]))), drop=FALSE]
      tmp$replicate_nonzero_instensity_count <- apply(subdata[1:ncol(subdata)], 1, function(x) length(which(x != 0)) )
      
      tmp <- tmp[tmp$replicate_nonzero_instensity_count >= ncol(subdata)/2,]
      
    }
    ### code modification stops
    tmp <- tmp[ ,c(2, grep("Reporter.intensity.corrected.[0-9].{1,}.", colnames(tmp))), drop=FALSE ]
  } else {
    tmp <- tmp[ ,c(2, grep("Intensity.", colnames(tmp))), drop=FALSE ]
  }
}

# In case of 10 TMT/iTRAQ samples check tmp and remove the coloumn 'Reporter.intensity.corrected.10'
#tmp <- tmp[,c(-2)]



#####
Majority.protein.IDs <- tmp$Majority.protein.IDs
tmp <- tmp[ , -1]
# for iTRAQ and TMT ppb normalization might not be the best method
# For TMT dataset PXD007160 normalisation, see NOTE at the end fo the script

if(EXP_TYPE == "MS1-quant"){
  tmp <- fot.normalise(tmp)  
}
#
tmp <- data.frame( cbind(Majority.protein.IDs, tmp, stringsAsFactors = FALSE) )
##############
tmp[tmp == 0] <- NA
tmp[ tmp == "NaN"] <- NA


##### perform the protein group to gene mapping
# this loop will take some time


data.to.map <- tmp
data.to.map$"ENSG" <- "NA"
data.to.map$"Gene.Symbol" <- "NA"
data.to.map$"unique.gene.count" <- "NA"

for(i in 1:nrow(data.to.map)){
    
    x <- data.frame(strsplit(data.to.map[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
    x_temp <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
    #x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
    #x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
    #x[,1] <- gsub("-\\d+", "", x[,1], perl=TRUE) # remove isoform number
    f = file()
    sink(file=f)
    res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("human")),
             error = function(e) {print(0)})
    # Covid2 taxonomy id: 2697049
    sink()
    close(f)
                 
    if (class(res)=="DFrame" | class(res) == "DataFrame"){
    data.to.map[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    data.to.map[ i, "Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")}
    temp_symb <- data.to.map[i,"Gene.Symbol"]
    data.to.map[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
    
    print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(data.to.map),1)),"%"))
  }



# remove protein groups that have no mapping to an ENSG gene IDs

#For SARS2, the ENSG are not present/mapped, so instead of leaving the field blank, copy Gene.Symbols
data.to.map$ENSG[grep("SARS2", data.to.map$Majority.protein.IDs, perl=TRUE)] <- data.to.map[grep("SARS2",data.to.map$Majority.protein.IDs, perl=TRUE), c("Majority.protein.IDs")] 
# This protein ID could not be mapped to gene symbols, therefore manually added
data.to.map$ENSG[data.to.map$Majority.protein.IDs == "sp|P0DTD2|ORF9B_SARS2"] <- "sp|P0DTD2|ORF9B_SARS2" 
data.to.map$Gene.Symbol[data.to.map$Majority.protein.IDs == "sp|P0DTD2|ORF9B_SARS2"] <- "9b"

#take a backup of data before filtering
data.to.map_before_filtering <- data.to.map

data.to.map <- data.to.map[ data.to.map$ENSG != "" , ]
data.to.map <- data.to.map[ data.to.map$ENSG != "NA" , ]
data.to.map <- data.to.map[ data.to.map[, "unique.gene.count"] == 1, ]

#tmp <- data.to.map
# remove all protein groups that map to multiple ENSG gene IDs (this removes a lot of proteins) - the reasoning to remove these cases is that we cannot establish for sure which gene is contributing the signal to the protein abundance; all genes contribute equally or one gene is a majority? 
data.to.map <- data.to.map[ grep(";", data.to.map$ENSG, invert = TRUE) , ]
# for genes that map to multiple proteins, in order to determine the amount of protein that gene is producing we sum the protein quantification values
#xx.Majority.protein.IDs <- aggregate(data.to.map$Majority.protein.IDs, list(ESNG = data.to.map$ENSG ), function(x) paste0( (x) )  )
#
colnames(data.to.map)
## Select which coloumns to aggregate 
data.to.map <- aggregate(data.to.map[ , 2:(ncol(data.to.map)-3)], list("Gene ID" = data.to.map$ENSG, "Gene.Symbol" = data.to.map$Gene.Symbol), median, na.rm =TRUE)

colnames(data.to.map) <- gsub("Reporter.intensity.corrected.", "", colnames(data.to.map), perl=TRUE)
# data.to.map <- cbind(Majority.protein.IDs = (as.character(xx.Majority.protein.IDs$x)), data.to.map )

#data.to.map <- data.to.map[,-c(13:252)]

if(EXP_TYPE == "MS2-quant"){
  write.table(data.to.map, "proteinGroups_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )} else {
    write.table(data.to.map, "proteinGroups_ppb_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )   
  }

#Write to file the list of protein groups which are mapped to more than 1 gene ID
ambiguous_gene_mapped_protein_groups <- data.to.map_before_filtering[ data.to.map_before_filtering[, "unique.gene.count"] > 1 | data.to.map_before_filtering[, "ENSG"] == "" | is.na(data.to.map_before_filtering$ENSG), ]
ambiguous_gene_mapped_protein_groups <- data.to.map_before_filtering[ data.to.map_before_filtering[, "unique.gene.count"] > 1 ,]
ambiguous_gene_mapped_protein_groups$ENSG <- unlist(lapply(ambiguous_gene_mapped_protein_groups$ENSG,function(x) paste(sort(strsplit(x, split=";")[[1]]), collapse=";")))
ambiguous_gene_mapped_protein_groups$Gene.Symbol <- unlist(lapply(ambiguous_gene_mapped_protein_groups$Gene.Symbol,function(x) paste(sort(strsplit(x, split=";")[[1]]), collapse=";")))
#check for non sample columns that start with ppb.iBAQ.
#ambiguous_gene_mapped_protein_groups <- ambiguous_gene_mapped_protein_groups[-c(2)]

ambiguous_gene_mapped_protein_groups <- ambiguous_gene_mapped_protein_groups[,c(1, ncol(ambiguous_gene_mapped_protein_groups)-2, ncol(ambiguous_gene_mapped_protein_groups)-1, ncol(ambiguous_gene_mapped_protein_groups), 2:(ncol(ambiguous_gene_mapped_protein_groups)-3) )]

write.table(ambiguous_gene_mapped_protein_groups, "ambiguous_gene_mapped_protein_groups.txt", sep = "\t", row.names = FALSE, quote = FALSE )

####### END ##########




# NOTE ###################
# For dataset PXD007160 (TMT) the "proteinGroups_final.txt" was further processed.
# Each batch was internally normalised with mean of the two Global Internal Standards (GIS-1 and GIS-2).
# The internal normalised file is named "proteinGroups_final_GIS_BatchNormalised.txt"

foo1 <- data.to.map
#data.to.map<- data.to.map[,c(-3)]

# NORMALISATION
# Batch1 AnteriorCingulateGyrus : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards 
data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus <- median(c(data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.1),
                                                                (data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.1),
                                                                na.rm=TRUE)
# Normalise channels 2 to 9 of Batch1 using the median values of the AnteriorCingulateGyrus Batch-1 GISs.
data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.1 <- data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.1/data.to.map$Median_Batch1_GIS_AnteriorCingulateGyrus

# Batch2 AnteriorCingulateGyrus : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus <- median(c(data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.2),
                                                                (data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.2),
                                                                na.rm=TRUE)
# Normalise channels 2 to 9 of Batch2 using the median values of the AnteriorCingulateGyrus Batch-2 GISs.
data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.2 <- data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.2/data.to.map$Median_Batch2_GIS_AnteriorCingulateGyrus

# Batch3 AnteriorCingulateGyrus : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus <- median(c(data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.3),
                                                               (data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.3),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch3 using the median values of the AnteriorCingulateGyrus Batch-3 GISs.
data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.3 <- data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.3/data.to.map$Median_Batch3_GIS_AnteriorCingulateGyrus

# Batch4 AnteriorCingulateGyrus : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus <- median(c(data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.4),
                                                               (data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.4),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch4 using the median values of the AnteriorCingulateGyrus Batch-4 GISs.
data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.4 <- data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.4/data.to.map$Median_Batch4_GIS_AnteriorCingulateGyrus

# Batch5 AnteriorCingulateGyrus : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus <- median(c(data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.5),
                                                               (data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.5),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch5 using the median values of the AnteriorCingulateGyrus Batch-5 GISs.
data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.1.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.2.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.3.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.4.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.5.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.6.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.7.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.8.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.9.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus
data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.5 <- data.to.map$Reporter.intensity.corrected.10.anterior.cingulate.gyrus_batch.5/data.to.map$Median_Batch5_GIS_AnteriorCingulateGyrus




# Batch1 FrontalCortex : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards 
data.to.map$Median_Batch1_GIS_FrontalCortex <- median(c(data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.1),
                                                               (data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.1),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch1 using the median values of the FrontalCortex Batch-1 GISs.
data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.1 <- data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.1/data.to.map$Median_Batch1_GIS_FrontalCortex

# Batch2 FrontalCortex : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch2_GIS_FrontalCortex <- median(c(data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.2),
                                                               (data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.2),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch2 using the median values of the FrontalCortex Batch-2 GISs.
data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.2 <- data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.2/data.to.map$Median_Batch2_GIS_FrontalCortex

# Batch3 FrontalCortex : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch3_GIS_FrontalCortex <- median(c(data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.3),
                                                               (data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.3),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch3 using the median values of the FrontalCortex Batch-3 GISs.
data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.3 <- data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.3/data.to.map$Median_Batch3_GIS_FrontalCortex

# Batch4 FrontalCortex : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch4_GIS_FrontalCortex <- median(c(data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.4),
                                                               (data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.4),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch4 using the median values of the FrontalCortex Batch-4 GISs.
data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.4 <- data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.4/data.to.map$Median_Batch4_GIS_FrontalCortex

# Batch5 FrontalCortex : Get median value of GIS-1 (channel1) and GIS-2 (channel10) internal standards
data.to.map$Median_Batch5_GIS_FrontalCortex <- median(c(data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.5),
                                                               (data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.5),
                                                               na.rm=TRUE)
# Normalise channels 2 to 9 of Batch5 using the median values of the FrontalCortex Batch-5 GISs.
data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.1.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.2.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.3.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.4.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.5.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.6.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.7.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.8.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.9.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex
data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.5 <- data.to.map$Reporter.intensity.corrected.10.frontal.cortex_batch.5/data.to.map$Median_Batch5_GIS_FrontalCortex

colnames(data.to.map)
# Remove the median values from dataframe
normalised_data <- data.to.map[,c(-103:-112)]

write.table(normalised_data, "proteinGroups_final_GIS_BatchNormalised.txt", sep = "\t", row.names = FALSE, quote = FALSE )

proteingroups_normalised <- normalised_data 
colnames(proteingroups_normalised) <- gsub("Reporter.intensity.corrected.", "", colnames(proteingroups_normalised))

# RENAMING CHANNELS
#Anterior cingulate gyrus
#Batch-1
colnames(proteingroups_normalised) <- gsub("1.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.anterior.cingulate.gyrus_batch.1", "Anterior cingulate gyrus_Batch-1. GIS-2", colnames(proteingroups_normalised))

#Batch-2
colnames(proteingroups_normalised) <- gsub("1.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.anterior.cingulate.gyrus_batch.2", "Anterior cingulate gyrus_Batch-2. GIS-2", colnames(proteingroups_normalised))

#Batch-3
colnames(proteingroups_normalised) <- gsub("1.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.anterior.cingulate.gyrus_batch.3", "Anterior cingulate gyrus_Batch-3. GIS-2", colnames(proteingroups_normalised))

#Batch-4
colnames(proteingroups_normalised) <- gsub("1.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.anterior.cingulate.gyrus_batch.4", "Anterior cingulate gyrus_Batch-4. GIS-2", colnames(proteingroups_normalised))

#Batch-5
colnames(proteingroups_normalised) <- gsub("1.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.anterior.cingulate.gyrus_batch.5", "Anterior cingulate gyrus_Batch-5. GIS-2", colnames(proteingroups_normalised))


#Frontal Cortex
#Batch-1
colnames(proteingroups_normalised) <- gsub("1.frontal.cortex_batch.1", "Frontal cortex_Batch-1. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.frontal.cortex_batch.1", "Frontal cortex_Batch-1. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.frontal.cortex_batch.1", "Frontal cortex_Batch-1. GIS-2", colnames(proteingroups_normalised))

#Batch-2
colnames(proteingroups_normalised) <- gsub("1.frontal.cortex_batch.2", "Frontal cortex_Batch-2. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.frontal.cortex_batch.2", "Frontal cortex_Batch-2. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.frontal.cortex_batch.2", "Frontal cortex_Batch-2. GIS-2", colnames(proteingroups_normalised))

#Batch-3
colnames(proteingroups_normalised) <- gsub("1.frontal.cortex_batch.3", "Frontal cortex_Batch-3. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.frontal.cortex_batch.3", "Frontal cortex_Batch-3. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.frontal.cortex_batch.3", "Frontal cortex_Batch-3. GIS-2", colnames(proteingroups_normalised))

#Batch-4
colnames(proteingroups_normalised) <- gsub("1.frontal.cortex_batch.4", "Frontal cortex_Batch-4. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.frontal.cortex_batch.4", "Frontal cortex_Batch-4. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.frontal.cortex_batch.4", "Frontal cortex_Batch-4. GIS-2", colnames(proteingroups_normalised))

#Batch-5
colnames(proteingroups_normalised) <- gsub("1.frontal.cortex_batch.5", "Frontal cortex_Batch-5. GIS-1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("2.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("3.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("4.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Control 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("5.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Control 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("6.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("7.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("8.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's & Parkinson's disease 1", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("9.frontal.cortex_batch.5", "Frontal cortex_Batch-5. Alzheimer's & Parkinson's disease 2", colnames(proteingroups_normalised))
colnames(proteingroups_normalised) <- gsub("10.frontal.cortex_batch.5", "Frontal cortex_Batch-5. GIS-2", colnames(proteingroups_normalised))

write.table(proteingroups_normalised, "proteinGroups_final_GIS_BatchNormalised-ChannelsRenamed.txt", sep = "\t", row.names = FALSE, quote = FALSE )


proteingroups_modify <- proteingroups_normalised[ -c(1,2) ]
meltData <- melt(proteingroups_modify)
meltData$Batch <- meltData$variable
meltData$Batch <- gsub("\\.\\s+.*", "", meltData$Batch, perl=TRUE)
meltData$Batch <- gsub(".*_", "", meltData$Batch, perl=TRUE)
meltData$Condition <- meltData$variable
meltData$Condition <- gsub("Batch-\\d+\\.\\s", "", meltData$Condition, perl=TRUE)
meltData$Condition <- gsub(".*_", "", meltData$Condition, perl=TRUE)
meltData$Tissue <- meltData$variable
meltData$Tissue <- gsub("_.*", "", meltData$Tissue, perl=TRUE)


# Plot
ggplot( meltData, aes(x=Condition, y=value)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("PXD007160 (TMT-10plex); MaxQuant output\nAfter Global Internal Standard (GIS) Normalisation (Median) per Batch") +
  xlab("")+
  ylab("Reporter intensity corrected")+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(Tissue ~ Batch)
  