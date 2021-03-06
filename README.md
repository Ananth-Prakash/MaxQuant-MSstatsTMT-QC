# MaxQuant-MSstatsTMT-QC
For TMT/iTRAQ experiments

The directory includes scripts for post-processing raw output from MaxQuant, calculating foldchanges using MSstatsTMT and Quality Assessment script for generating pdf document.

These files should be run in this sequential order.

1_process-protein-groups-with-gene-names.R

2_MSstats_MaxQuant-TMT-iTRAQ-pipeline.R

3_MaxQuant_Summary_ExpressionAtlas_Differential_TMT-iTRAQ_Rmarkdown.Rmd

The directory inlcudes all input and output files of PXD012203 (iTRAQ) for running these scripts.

To note: This QC script is specific to this dataset, which is different from the QC script for baseline experiments, where one script can be used for all datasets. The reasons are i) while some TMT and iTRAQ experiments include channels for Global Internal Standard (GIS) which are used for normalisation, some do not, and therefore scripts need to be modified in those cases, ii) also factors that needs comparing to compute fold-changes are different fir different experiments and they need to be modified within the code.

To generate the IDF, make use of Annotare https://www.ebi.ac.uk/fg/annotare/login/

