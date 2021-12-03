# MSstats tutorial for post-processing MaxQuant outputs - TMT/iTRAQ

# Webpage  - http://msstats.org/msstats-2/
# Tutorial - TMT/iTRAQ: https://www.bioconductor.org/packages/release/bioc/vignettes/MSstatsTMT/inst/doc/MSstatsTMT.html
# Tutorial - LabelFree: http://msstats.org/wp-content/uploads/2019/11/MSstats_v3.18.1_manual.pdf

# http://msstats.org/msstats-2/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstatsTMT")
library('MSstatsTMT', warn.conflicts = F, quietly = T, verbose = F)



# https://bioconductor.org/packages/release/bioc/manuals/mygene/man/mygene.pdf
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
 ##install.packages("remotes")
 ##remotes::install_github("Bioconductor/GenomeInfoDb")
 ##remotes::install_github("Bioconductor/GenomicFeatures")
#BiocManager::install("mygene")
library(mygene)


library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)

setwd('/Users/ananth/Documents/MSstats_scripts_for_Andrew/')


#### 1. Read proteinGroups to get proteinID information.
proteinGroups_inp  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#### 2. Read evidence.txt file
evidence_file  <- read.table( "evidence.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#### 3. Read annotation file
#### Note: In the user defined the annotation file the RAW file names in coloumn 'Run' should not have .raw extension!
####       the Channel coloumn is case sensitive, and entries must only be in the format 'channel.1', 'channel.2', 'channel.3', ......
annot <- read.table("MSstat_annot_factors.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


##### 4. Read sample replicate information
##### Note: User defined annotation file
sample_replicates <- read.table("sample_replicates.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#####A. First we filter the proteinGroups input file before processing with MSstatsTMT
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "Reverse"] != "+", ]
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "Potential.contaminant"] != "+", ]
## consider those protein groups with more than 1 peptides matched
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "Peptides"] > 1, ]


# # Only consider those protein groups for which at least 50% of the sample replicates have non-zero intensity values

for(i in 1:nrow(sample_replicates))
{
  sample_replicates[i,] <- gsub("^","Reporter.intensity.corrected.", sample_replicates[i,], perl=TRUE)
  sample_replicates[i,] <- gsub(" ",".", sample_replicates[i,], perl=TRUE)
  
  subdata <- proteinGroups_inp[, c(unname(unlist(sample_replicates[i,-1])))]
  proteinGroups_inp$replicate_nonzero_instensity_count <- apply(subdata[1:ncol(subdata)], 1, function(x) length(which(x != 0)) )
  
  proteinGroups_inp <- proteinGroups_inp[proteinGroups_inp$replicate_nonzero_instensity_count >= ncol(subdata)/2,]
  
}


#####B. Add gene name and protein name
#### perform protein group to gene mapping
#### this loop will take some time

proteinGroups_inp$"ENSG" <- "NA"
proteinGroups_inp$"Gene.Symbol" <- "NA"
proteinGroups_inp$"Gene.Name" <- "NA"
proteinGroups_inp$"unique.gene.count" <- "NA"

for(i in 1:nrow(proteinGroups_inp))
{
  x <- data.frame(strsplit(as.character(proteinGroups_inp[ i, "Majority.protein.IDs"]), split = ";"), stringsAsFactors = FALSE)
  
  # x[,1] <- gsub("sp\\||tr\\|", "", x[,1], perl=TRUE)
  # x[,1] <- gsub("\\|.*", "", x[,1], perl=TRUE)
  # x[,1] <- gsub("-\\d+", "", x[,1], perl=TRUE) # remove isoform number
  
  x[,1] <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))
  
  f = file()
  sink(file=f)
  res <- tryCatch(queryMany(x[,1], scopes="uniprot", fields=c("ensembl.gene", "symbol", "name"), species="human"), error = function(e) {print(0)})
  sink()
  close(f)
  
  #print(res)
  
  if (class(res)=="DFrame" | class(res) == "DataFrame")
  {
    proteinGroups_inp[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    proteinGroups_inp[ i, "Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")
    proteinGroups_inp[ i, "Gene.Name"] <- paste( unique(unlist(res$name[!is.na(res$name)])), collapse = ";") 
    #print(proteinGroups_inp[i,"Gene.Symbol"])
    temp_symb <- proteinGroups_inp[i,"Gene.Symbol"]
    proteinGroups_inp[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
  }
  print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(proteinGroups_inp),1)),"%"))
} #end of for


#####C Protein groups that are mapped to a unique gene
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp$ENSG != "" , ]
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp$ENSG != "NA" , ]
proteinGroups_inp <- proteinGroups_inp[ proteinGroups_inp[, "unique.gene.count"] == 1, ]
proteinGroups_inp <- proteinGroups_inp[ grep(";", proteinGroups_inp$ENSG, invert = TRUE) , ]

colnames(proteinGroups_inp)
# Although the first step of postprocessing  MaxQuant output in "process-protein-groups-with-gene-names.R"
# does aggregation, here it is not advised as it changes MaxQuant output format and MSstatsTMT does not recognise
# change in format.
# There will slight change in number of rows before and after aggregation. 
# Normalised intensity values in examples such as TMT PXD007160 may slightly differ some entries.
## Aggregate the columns "Reporter.intensity.corrected". Check their column numbers
# data.to.map <- aggregate(data.to.map[ , 122:221], list("Gene ID" = data.to.map$ENSG, "Gene.Name" = data.to.map$Gene.Name), sum, na.rm =TRUE)



# Copy dataframe here, normalise in case of TMT using the "process-protein-groups-with-gene-names.R"
# script and use the normalised proteinGroups_inp as input to MSstatsTMT

#### Start MSstatsTMT
#### I . Convert input files into MSstat format and do initial filtering
input.mq <- MaxQtoMSstatsTMTFormat(evidence_file, proteinGroups_inp, annot,
                                   rmPSM_withMissing_withinRun = TRUE)


#### II Protein summarization.
quant.msstats <- proteinSummarization(input.mq,
                                      method="msstats",
                                      global_norm = TRUE,
                                      reference_norm = TRUE,
                                      remove_norm_channel = TRUE,
                                      remove_empty_channel = TRUE,
                                      MBimpute = FALSE
                                      #imputation should be set to false
                                      )


#### III Test for all the possible pairs of conditions
#### For pairwise comparisons among all conditions: 
# test.pairwise <- groupComparisonTMT(quant.msstats)

#### If comparisons have to be made between certain conditions then specify
#### Check how many conditions are there to be compared. 
levels(quant.msstats$Condition)
# [1] "EarlyOnset_Alzheimer's", "EarlyOnset_Non_demented", "LateOnset_Alzheimer's", "LateOnset_Non_demented"

#### The numerator is set to 1 and denominator is set to -1
#### All other places in the matrix related to other conditions are set to 0

# comparison between "EarlyOnset_Alzheimer's" and "EarlyOnset_Non_demented"
#### Here for example "Control" is set to denominator and hence -1
comparison_1 <- matrix(c(1,-1,0,0),nrow=1)
# comparison between "LateOnset_Alzheimer's" and "LateOnset_Non_demented"
comparison_2 <- matrix(c(0,0,1,-1),nrow=1)



#### Set the column names. (Has to be in the same order as levels(quant.msstats$Condition) )
colnames(comparison_1) <- levels(quant.msstats$Condition)
colnames(comparison_2) <- levels(quant.msstats$Condition)

#### Set the row name for the matrix
row.names(comparison_1)<-"EarlyOnset_Alzheimer's/EarlyOnset_Non_demented"
row.names(comparison_2)<-"LateOnset_Alzheimer's/LateOnset_Non_demented"


test.comparison_1 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_1)
test.comparison_2 <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison_2)


# The output (test.comparison) has multiple entries of same proteins because of their association with multiple protein groups.
# To filter out duplicate entries and to retain only one protein group, we retain that entry which has all the proteins that are also the majority proteins
# ie., those proteins that have > 50% of the peptides belonging to that protein group.
deduped_proteingroups_1 <- merge(x=test.comparison_1, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                               by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                               all.x=FALSE, all.y=FALSE) 

deduped_proteingroups_2 <- merge(x=test.comparison_2, y=data.frame(proteinGroups_inp[,c('Majority.protein.IDs','ENSG','Gene.Symbol','Gene.Name')]),
                                 by.x=c('Protein'), by.y=c('Majority.protein.IDs'),
                                 all.x=FALSE, all.y=FALSE) 



group1 <- deduped_proteingroups_1[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group1) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "EarlyOnset_Alzheimer's vs. EarlyOnset_Non_demented.foldChange", 
                      "EarlyOnset_Alzheimer's vs. EarlyOnset_Non_demented.pValue")

group2 <- deduped_proteingroups_2[,c("Protein", "ENSG", "Gene.Symbol", "Gene.Name", "log2FC", "adj.pvalue")]
colnames(group2) <- c("Protein", "Gene ID", "Gene Name", "Gene description", 
                      "LateOnset_Alzheimer's vs. LateOnset_Non_demented.foldChange",
                      "LateOnset_Alzheimer's vs. LateOnset_Non_demented.pValue")


combined_foldchange_file <- merge(x=group1, y=group2,
                                  by.x=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  by.y=c("Protein", "Gene ID", "Gene Name", "Gene description"),
                                  all.x=FALSE, all.y=FALSE) 


write.table(combined_foldchange_file, "FoldChange_groupcomparison.txt", sep = "\t", row.names = FALSE, quote = FALSE )


##### Additional processing
## This section looks at postprocessd output from 
# 1) FOT normalisation and 
# 2) MS-StatsTMT postprocessing after calculating log2 foldchange
# to keep same gene ids common between files.
# This setp is done as it is desired to submit to Expression Atlas in this format.

# NOTE: Postprocessing for FOT normalisation should be done and the output file be kept ready
# before running this step.

FOT_normalised_file <- read.table( "proteinGroups_final.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
log2_FC_file <- read.table( "FoldChange_groupcomparison.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
# Remove NA Pvalues
log2_FC_file <- log2_FC_file[complete.cases(log2_FC_file),]

FOT_genes <- FOT_normalised_file[,c("Gene.ID", "Gene.Symbol")]
log2_FC_genes <- log2_FC_file[,c("Gene.ID", "Gene.Name")]

# get common genes
merged_ids <- merge(x=FOT_genes, y=log2_FC_genes,
                    by.x=c("Gene.ID", "Gene.Symbol"), by.y=c("Gene.ID", "Gene.Name"),
                    all.x=FALSE, all.y=FALSE)

# filter these genes from both files
# the resulting file will have common genes between them.
FOT_normalised_file <- merge(x=merged_ids, y=FOT_normalised_file,
                             by.x=c("Gene.ID", "Gene.Symbol"),
                             by.y=c("Gene.ID", "Gene.Symbol"),
                             all.x=FALSE, all.y=FALSE)

log2_FC_file <- merge(x=merged_ids, y=log2_FC_file,
                             by.x=c("Gene.ID", "Gene.Symbol"),
                             by.y=c("Gene.ID", "Gene.Name"),
                             all.x=FALSE, all.y=FALSE)


write.table(FOT_normalised_file, "ExpressionAtlas_proteinGroups_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )
write.table(log2_FC_file, "ExpressionAtlas_FoldChange_groupcomparison.txt", sep = "\t", row.names = FALSE, quote = FALSE )
