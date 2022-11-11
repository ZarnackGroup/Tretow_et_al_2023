###########################################################################
#                                                                         #
# In this R script expression in BRCA samples is calculated by first      #
#   fetching a RangedSummarizedExperiment object via recount3, followed   #
#   by calculation of read counts, normalization,vst-transformation and   #
#   batch-correction. This script runs quite long.                        #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(recount3)
library(DESeq2)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

tcgaDir <- paste0(projectDir, "/8_TCGA_analysis")

#############################################################
#          Fetch RangedSummarizedExperiment object          #
#############################################################

humanProjects <- available_projects()

projectBRCA <- dplyr::filter(humanProjects, project == "BRCA" &
                                 file_source == "tcga")

# This call creates the RangedSummarizedExperiment, which has
#   63.856 rows (Genes) and 1.256 columns (samples)
rse <- create_rse(projectBRCA)

##########################################
#          Filtering of samples          #     
##########################################

# Keep samples of type Primary Tumor and Solid Tissue Normal obtained
#   from female participants

sampleType <- rse$tcga.gdc_cases.samples.sample_type
gender <- rse$tcga.gdc_cases.demographic.gender

keep <- which(sampleType %in% c("Primary Tumor", "Solid Tissue Normal") &
                  gender == "female")

rse <- rse[, keep]

# Yet there are some cases where I have multiple samples from the same
#   participant. In the case of one Primary Tumor and one Solid Tissue Normal
#   sample everything is fine. In rare cases (e.g. TCGA-A7-A0DB) ther
#   are for example 3 Primary Tumor and 1 Solid Tissue Normal samples.
#   These cases are removed as these cases are rare and hard to handle.

submitterID <- rse$tcga.gdc_cases.submitter_id

# These are participants with more than one sample and potential candidates
#   for removal
removalCandidates <- which(submitterID %>% table > 1) %>% names

# For all participiants with more than one sample I check if there is exactly
#   1 Primary Tumor and 1 Solid Tissue Normal sample. If this is the case
#   the participant is removed from the removal list
removalCandidates <- rse %>%
    colData %>%
    as.data.frame %>%
    dplyr::filter(., tcga.gdc_cases.submitter_id %in% removalCandidates) %>%
    # Count for each participant the number of samples per sample type
    dplyr::count(tcga.gdc_cases.submitter_id,
                 tcga.gdc_cases.samples.sample_type) %>%
    # Group counts by patient
    group_by(tcga.gdc_cases.submitter_id) %>%
    # Here for each participant it is define if he should be kept on the removal
    #   list (TRUE) or not (FALSE)
    summarize(
        tcga.gdc_cases.submitter_id = tcga.gdc_cases.submitter_id[1],
        remove=case_when(
            # more than two sample types  
            n() != 2 ~ TRUE,
            # one of the two sample types has more than one sample
            (n[1] != 1 | n[2] != 1) ~ TRUE,
            # at least one of the sample types is not Primary Tumor or Solid Tissue Normal
            any(!tcga.gdc_cases.samples.sample_type %in%
                    c("Primary Tumor", "Solid Tissue Normal")) ~ TRUE,
            TRUE ~ FALSE)) %>%
    dplyr::filter(remove) %>% 
    pull(tcga.gdc_cases.submitter_id)
    
#remove all samples of the removal participants 
rse <- rse[,!(submitterID %in% removalCandidates)]

#############################################################
#          Extract read counts and prepare colData          #
#############################################################

readCounts <- compute_read_counts(rse) 

# Set TCGA barcode as column name
colnames(readCounts) <- rse$tcga.tcga_barcode

colData <- data.frame(row.names = rse$tcga.tcga_barcode,
                      sampleType = rse$tcga.gdc_cases.samples.sample_type,
                      batchID=factor(rse$tcga.cgc_case_batch_number))

# Rename sample types
colData$sampleType <- colData$sampleType %>%
    gsub("Primary Tumor", "tumor", .) %>%
    gsub("Solid Tissue Normal", "normal", .)

############################################################################
#          Normalization, vst-transformation and batch correction          #
############################################################################

# Create a dds object
dds <- DESeqDataSetFromMatrix(countData = readCounts, colData = colData,
                              design = ~ batchID + sampleType)

# vst-transformation
vsd <- vst(dds, blind=FALSE)

# Batch correction  
vsdBatchCorrected <- vsd
mat <- assay(vsdBatchCorrected)
mm <- model.matrix(~sampleType, colData(vsdBatchCorrected))
mat <- limma::removeBatchEffect(mat, batch=vsdBatchCorrected$batchID, design=mm)
assay(vsdBatchCorrected) <- mat

# Save 
metaInformation <- colData(rse) %>% as.data.frame
rownames(metaInformation) <- c()

saveRDS(vsd, paste0(tcgaDir,"/rds_files/BRCA_vsd.rds"))
saveRDS(vsdBatchCorrected, paste0(tcgaDir,"/rds_files/BRCA_vsdBatchCorrected.rds"))
saveRDS(metaInformation, paste0(tcgaDir,"/rds_files/BRCA_metaInformation.rds"))
