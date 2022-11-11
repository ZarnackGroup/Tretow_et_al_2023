###########################################################################
#                                                                         #
# In this R script the adjusted Rand Index for the separation of          #
#   normal from tumor samples and for the separation for the              #
#   Basal subtype from the other subtypes should be calculated            #
#   for the regulated CE and ALE events as well as for sampled            #
#   non-regulated CE events. This script runs for multiple hours.                                               #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(missMDA)
library(matrixStats)
library(TCGAbiolinks)
library(mclust)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

tcgaDir <- paste0(projectDir, "/8_TCGA_analysis")

##################################
#          Prepare data          #
##################################

PSIs <- readRDS(paste0(tcgaDir,"/rds_files/BRCA_PSI.rds"))

regulatedPSIs <- PSIs %>% dplyr::filter(hillCat %in% c("Coop-Enh", "Coop-Rep"))
nonRegulatedPSIs <- PSIs %>% dplyr::filter(eventType  == "nonRegulatedCE")

#############################################################################
#          Remove events that are rarely quantified or have low dynamic     #
#############################################################################

keep1 <- regulatedPSIs %>%
    dplyr::count(eventID) %>%
    dplyr::filter(n >= 100) %>%
    pull(eventID)

keep2 <- regulatedPSIs %>%
    dplyr::count(eventID, PSI) %>%
    dplyr::count(eventID) %>%
    dplyr::filter(n >= 10) %>%
    pull(eventID)

keep <- intersect(keep1, keep2)

regulatedPSIs <- regulatedPSIs %>% dplyr::filter(eventID %in% keep)


keep1 <- nonRegulatedPSIs %>%
    dplyr::count(eventID) %>%
    dplyr::filter(n >= 100) %>%
    pull(eventID)

keep2 <- nonRegulatedPSIs %>%
    dplyr::count(eventID, PSI) %>%
    dplyr::count(eventID) %>%
    dplyr::filter(n >= 10) %>%
    pull(eventID)

keep <- intersect(keep1, keep2)

nonRegulatedPSIs <- nonRegulatedPSIs %>% dplyr::filter(eventID %in% keep)

########################################################################
#          Calculate Rand Indices for normal sample and Basal          # 
#               subtype separation                                     #
########################################################################


computeAdjustedRandIndices <- function(PSIs, setLabel){
    
    mat <- PSIs %>%
        dplyr::select(eventID, tcgaSampleBarcode, PSI) %>%
        arrange(eventID, tcgaSampleBarcode) %>%
        pivot_wider(names_from=tcgaSampleBarcode, values_from = PSI) %>% 
        column_to_rownames("eventID") %>%
        as.matrix
    
    mat <- imputePCA(mat %>% t, ncp=2)
    mat <- mat$completeObs %>% t
    
    # Z-score transformation
    matZscore <- (mat  - (mat %>% rowMeans(., na.rm=T))) / rowSds(mat, na.rm=T)
    
    sampleType <- PSIs$sampleType[
        match(colnames(matZscore), PSIs$tcgaSampleBarcode)]
    
    # PAM50 subtype via PanCancerAtlas_subtypes() of the TCGAbiolinks package
    PAM50subtype <- PanCancerAtlas_subtypes() %>%
        dplyr::filter(cancer.type=="BRCA")
    PAM50subtype <- PAM50subtype$Subtype_mRNA[
        match(colnames(matZscore),PAM50subtype$pan.samplesID)]
    
    set.seed(1)
    
    colClusters <- kmeans(matZscore %>% t, centers=5)
    
    
    df <- data.frame(cluster=colClusters$cluster,
                     sampleType=sampleType,
                     subType=PAM50subtype)
    
    #Identify the "normal" cluster - the one with the highest number of
    #   normal samples
    normalCluster <- df %>% 
        dplyr::count(cluster, sampleType) %>%
        dplyr::filter(sampleType == "Solid Tissue Normal") %>%
        arrange(desc(n)) %>% dplyr::slice(1) %>% pull(cluster)
    
    #Identify the "Basal" cluster - the one with the highest number of
    #   Basal samples
    basalCluster <- df %>% 
        dplyr::count(cluster, subType) %>%
        dplyr::filter(subType == "Basal") %>%
        arrange(desc(n)) %>% dplyr::slice(1) %>% pull(cluster)
    
    # Rename the clusters as either being the normal (C_Normal) or another 
    #   cluster (C_Other)
    df <- df %>%
        mutate(cluster_new=ifelse(cluster == normalCluster, "C_Normal", "C_Other"))
    
    # Calculate the first adjusted Rand index
    normalAdjustedRandIndex <- adjustedRandIndex(df$sampleType, df$cluster_new)
    
    # Rename the clusters as either being the Basal (C_Basal) or another 
    #   cluster (C_other). Similarly renme the subtypes as being either
    #   Basal or another subtype (other)
    df <- df %>%
        mutate(cluster_new=ifelse(cluster == basalCluster, "C_Basal", "C_Other")) %>%
        mutate(subType_new=ifelse(subType == "Basal", "Basal", "Other"))
    
    # Calculate the second adjusted Rand index
    basalAdjustedRandIndex <- adjustedRandIndex(df$subType_new, df$cluster_new)
    
    return(data.frame(setLabel, normalAdjustedRandIndex, basalAdjustedRandIndex))
    
}

# Compute the two adjusted Rand indices for the regulated events
adjustedRandIndicesReg <- computeAdjustedRandIndices(regulatedPSIs, setLabel="reg")

# Randomly draw 100 times a similar number of non-regulated events and
#   compute the two adjusted Rand indices
set.seed(123)
randomNonRegSets <- lapply(1:100, function(i){
    nonRegEventIDs <- nonRegulatedPSIs$eventID %>% unique
    nonRegEventIDs[sample(x=1:length(nonRegEventIDs),
                          size=regulatedPSIs$eventID %>% unique %>% length)]})

adjustedRandIndicesNonReg <- lapply(randomNonRegSets, function(eventIDs){
    subSetNonRegulatedPSIs <- nonRegulatedPSIs %>% dplyr::filter(eventID %in% eventIDs)
    return(computeAdjustedRandIndices(subSetNonRegulatedPSIs, setLabel="nonReg"))
}) %>% bind_rows()

adjustedRandIndices <- rbind(adjustedRandIndicesReg, adjustedRandIndicesNonReg)

# Save
saveRDS(adjustedRandIndices, paste0(tcgaDir,"/rds_files/adjustedRandIndices.rds"))
    