###########################################################################
#                                                                         #
# In this R script PSI values in BRCA samples are calculated by first     #
#   fetching a RangedSummarizedExperiment object via recount3, followed   #
#   by filtering for samples and junctions of interest. At the end PSI    #
#   values are calculated as PSI = IJcounts / (IJcounts + SJcounts).      #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(recount3)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

tcgaDir <- paste0(projectDir, "/8_TCGA_analysis")

#############################################################
#          Fetch RangedSummarizedExperiment object          #
#############################################################

# Load the sample meta information created in Calculate_expression_BRCA.R
metaInformation <- readRDS(paste0(tcgaDir,"/rds_files/BRCA_metaInformation.rds"))

humanProjects <- available_projects()
projectBRCA <- dplyr::filter(humanProjects, project == "BRCA" & file_source == "tcga")

# This call creates the RangedSummarizedExperiment, which has
#   4,703,960 (junctions) and 1.256 columns (samples)
rse <- create_rse(projectBRCA, type = "jxn", jxn_format="UNIQUE")

###############################
#          Filtering          #     
###############################

# Subset to samples that survived the filtering in Calculate_expression_BRCA.R
rse <- rse[,rse$tcga.tcga_barcode %in% metaInformation$tcga.tcga_barcode]

# Load the junctions of interest 
jxGRanges <- readRDS(paste0(tcgaDir,"/rds_files/jxGRanges.rds"))

# recount junctions start at the first and end at the last
#   intron position. Need to adjust the junctions of interest.
jxGRanges <- jxGRanges-1

# reduce the TCGA junctions to the ones of interest by accessing their
#   ranges (rowRanges() returns a GRanges object). 3.021 of 3.089 junctions
#   are present in the RangedSummarizedExperiment
rse <- rse[rowRanges(rse) %in% jxGRanges,]

##################################
#          Compute PSIs          #     
##################################

# Split junctions GRanges objects of junction pairs by eventIDs
jxGRangesList <- split(jxGRanges, jxGRanges$eventID)

# Extract the junction counts
junctionCounts <- assays(rse)$counts %>% as.matrix

# For each entry (CE or ALE event) calculate the PSI per sample
#   Note that for the HNRNPH1 NMD isoform the first and second
#   inclusion junction are used, while for all other events
#   only the first one is used.
PSIs <- lapply(1:length(jxGRangesList), function(i){
    
    # Grep the GRanges object with the two junctions
    gr <- jxGRangesList[[i]]
    
    # Determine the matches in the RangedSummarizedExperiment object
    rseMatches <- match(gr, rowRanges(rse))
    
    # If one or both junctions have no match return an empty data.frame
    if(any(is.na(rseMatches))){
        return(data.frame())
    }
    
    # Extract the counts of the two junctions
    tmpCounts <- junctionCounts[rseMatches,]
    
    # Subset the RangedSummarizedExperiment object
    tmpRse <- rse[rseMatches,]
    
    # Determine samples (columns) with enough counts 
    keep <- colSums(tmpCounts) >= 25
    
    # If there are less than 10 samples with enough counts return an empty
    #   data.frame
    if(sum(keep) < 10){
        return(data.frame())
    }
    
    # Remove the samples that do not have enough counts
    tmpCounts <- tmpCounts[,keep]
    tmpRse <- tmpRse[, keep]
    
    # For the HNRNPH1 NMD isoform both inclusion junctions are used
    if(gr$eventType[1] != "HNRNPH1nmdIsoform"){
        PSIs <- prop.table(tmpCounts,2)[which(gr$jxType == "IJ"),]
        tcgaPatientBarcodes <- tmpRse$tcga.gdc_cases.submitter_id
        tcgaSampleBarcodes <- tmpRse$tcga.tcga_barcode
        sampleTypes <- tmpRse$tcga.gdc_cases.samples.sample_type 
        return(
            data.frame(eventID=gr$eventID[1],
                       hillCat=gr$hillCat[1],
                       geneName=gr$geneName[1],
                       geneID=gr$geneID[1],
                       eventType=gr$eventType[1],
                       PSI=PSIs,
                       tcgaPatientBarcode=tcgaPatientBarcodes,
                       tcgaSampleBarcode=tcgaSampleBarcodes,
                       sampleType = sampleTypes)
        )
    } else {
        PSIs <- tmpCounts %>%
            # col[1] is the IJ1, col[2] is the IJ2 and col[3] the SJ
            apply(., 2, function(col){
                (col[1] + col[2]) / (col[1] + col[2] + col[3] *2)
                })
        tcgaPatientBarcodes <- tmpRse$tcga.gdc_cases.submitter_id
        tcgaSampleBarcodes <- tmpRse$tcga.tcga_barcode
        sampleTypes <- tmpRse$tcga.gdc_cases.samples.sample_type 
        return(
            data.frame(geneName=gr$geneName[1],
                       geneID=gr$geneID[1],
                       eventType=gr$eventType[1],
                       PSI=PSIs,
                       tcgaPatientBarcode=tcgaPatientBarcodes,
                       tcgaSampleBarcode=tcgaSampleBarcodes,
                       sampleType = sampleTypes)
        )
    }
    
}) %>% bind_rows()

rownames(PSIs) <- c()

# Save
saveRDS(PSIs, paste0(tcgaDir,"/rds_files/BRCA_PSI.rds"))
    