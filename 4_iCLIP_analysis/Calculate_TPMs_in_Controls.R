###########################################################################
#                                                                         #
# In this R script for the two control experiments (knockdown and         #
#   overexpression) TPM values are computed, which will be used in        #
#   upcoming analyses to define expressed genes.                          #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(DESeq2)
library(GenomicFeatures)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

datarDir <- paste0(projectDir, "/data")
iClipDir <- paste0(projectDir, "/4_iCLIP_analysis")

###############################
#          Load data          #
###############################

# Extract unnormalized counts for knockdown experiments
load(paste0(datarDir, "/DESeq2/DESeq2_knockdown.RData"))
KD_counts <- counts(dds, normalized=FALSE)
rm(conts, dds, gtf, pairwise.dds, res, rld, add_factors)


# Extract unnormalized counts for knockdown experiments
load(paste0(datarDir, "/DESeq2/DESeq2_overexpression.RData"))
OE_counts <- counts(dds, normalized=FALSE)
rm(conts, dds, pairwise.dds, res, rld, add_factors)

# Extract the control replicates
KD_counts <- KD_counts[, c("ctrl_rep1", "ctrl_rep2", "ctrl_rep3")]
OE_counts <- OE_counts[, c("Ctrl_rep1", "Ctrl_rep2", "Ctrl_rep3")]

# Merge knockdwon and overexpression
counts <- cbind(KD_counts, OE_counts)

##################################
#          Computations          #
##################################

# The gtf object was loaded as part of the the DESeq2 output of the
#   overexpression and is the same as the one loaded from the knockdown DESeq2
#   output
txdb <- makeTxDbFromGRanges(gtf)

# Extract all exons of a gen
exonsByGene <- exonsBy(txdb, by = "gene")

# reduce overlapping exons to a single region
reducedExonsByGene <- reduce(exonsByGene)

# Approximate the gene length as the sum of the its exons lengths
geneLengths <- sum(width(reducedExonsByGene))

# Re-order the vector of gene lengths to match the order in the counts
geneLengths <- geneLengths[match(rownames(counts), names(geneLengths))]

# First step1 in TPM calculation
controlTPMs <- counts / geneLengths

# Second step in TPM calculation
controlTPMs <- t(t(controlTPMs) * 1e6 / colSums(controlTPMs)) %>% as.data.frame

# Create a final data.frame with the average TPM of knockdown and overexpression
#   TPMs
controlTPMs <- controlTPMs %>%
    tibble::rownames_to_column() %>%
    dplyr::rename(geneID=rowname) %>%
    rowwise() %>%
    mutate(KDmerge=mean(c(ctrl_rep1, ctrl_rep2, ctrl_rep3))) %>% 
    mutate(OEmerge=mean(c(Ctrl_rep1, Ctrl_rep2, Ctrl_rep3))) %>% 
    mutate(Allmerge=mean(c(ctrl_rep1, ctrl_rep2, ctrl_rep3,
                           Ctrl_rep1, Ctrl_rep2, Ctrl_rep3))) %>%
    ungroup

# The compute TPMs are stored in a RDS-File
saveRDS(controlTPMs, paste0(iClipDir,"/rds_files/controlTPMs.rds"))
