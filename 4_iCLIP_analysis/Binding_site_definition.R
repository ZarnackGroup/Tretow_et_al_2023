###########################################################################
#                                                                         #
# In this R script HNRNPH binding sites are defined based on the PureCLIP #
#   output and crosslink events loaded from BigWig-Files.                 #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

dataDir <- paste0(projectDir, "/data")
iClipDir <- paste0(projectDir, "/4_iCLIP_analysis")

###############################
#          Load data          #
###############################

# Load PureCLIP crosslink sites as GRanges object
PureClipSites <- import(paste0(dataDir,
                                "/iCLIP/PureCLIP_HNRNPH_crosslink_sites.bed"),
                            extraCols=c(scoreInfo = "character"))

# Load the raw crosslink events from BigWig-Files as RleList and combine
#   the + and - strand in a list
crosslinkSites <- list(
    "+" = import(paste0(dataDir,"/iCLIP/HNRNPH_Xlinks_merge_plus.bw"),
                 as = "Rle"),
    "-" = import(paste0(dataDir,"/iCLIP/HNRNPH_Xlinks_merge_minus.bw"),
                 as = "Rle")
)

#################################################
#          Definition of binding sites          #
#################################################

# The following steps were performed:
#   1. Merge near-by PureCLIP sites into regions if there are <= 3 nucleotides
#      between them (min.gapwidth = 4)
#   2. Remove regions with width < 2
#   3. Center detection and extension
#       a. In each left region define a binding site by indentifying the center
#          (position with highest number of crosslink events) and extend it
#          by +/- 2 nts (=> width of 5)
#       b. Remove the binding site and the 4 nts up- and downstream from the
#          region, which generates zero, one or two new regions
#       c. Repeat steps 3a and 3b for the new regions until no region is left
#   4. Keep binding sites where:
#       a. the center is a PureCLIP site
#       b. the center has the highest number of crosslink events
#       c. there are at least 3 positions with a crosslink event


###########
#  Step 1 #
###########

regions <- reduce(PureClipSites, min.gapwidth = 4)

###########
#  Step 2 #
###########

regions <- regions[width(regions) > 1]

###########
#  Step 3 #
###########

# As the + and - strand need to be discriminated, the procedure is implemented
#   in a function that is use once for the + and once for the - strand
detectCenterAndExtend <- function(regions, crosslinkSites, strand) {
    
    # Extract the remaining regions on the current strand
    regionsRemaining <- regions[strand(regions) == strand]
    
    # Extract the crosslink sites for the current strand
    crosslinkSitesCurrentStrand <- crosslinkSites[[strand]]
    
    # Initialize an empty GRanges object for the binding sites
    bindingSites <- GRanges()
    
    # Set the window size 
    windowSize <- 2
    
    # Runs as long as there are regions left
    while(TRUE) {
        
        #No regions left to check for binding sites
        if (length(regionsRemaining) == 0) {
            break
        }
        
        #Create a matrix of crosslink events in the remaining regions
        crosslinkSitesInRemainingRegions <- as.matrix(
            crosslinkSitesCurrentStrand[regionsRemaining])
        
        # Define the ties method for the the + and - strand
        tiesMethod = ifelse(strand == "+", "first", "last")
        
        # Set NAs to -Inf
        crosslinkSitesInRemainingRegions[
            is.na(crosslinkSitesInRemainingRegions)] <- -Inf
        
        #Identify for each region (rows in the matrix) the position with the
        #   highest number of crosslink events
        maxPositionIndice <- max.col(crosslinkSitesInRemainingRegions,
                                     ties.method = tiesMethod)
        
        #Create a 5 nt binding site centered at the max position
        bindingSite <- regionsRemaining
        start(bindingSite) <- start(bindingSite) + maxPositionIndice - 1
        end(bindingSite) <- start(bindingSite)
        bindingSite <- bindingSite + windowSize
        
        #Store the new binding site
        bindingSites <- c(bindingSites, bindingSite)
        
        
        #Remove the binding sites and the 4 nts up- and downstream from the 
        #   remaining regions. This leads to zero, one or two new regions for
        #   each previous region
        regionToRemove <- as(bindingSite+4, "GRangesList")
        regionsRemaining <- unlist(psetdiff(regionsRemaining,
                                                   regionToRemove))
    }
    return(bindingSites)
}

# Apply the function on the + and - strand
bindingSites <- c(
    detectCenterAndExtend(regions, crosslinkSites, "+"),
    detectCenterAndExtend(regions, crosslinkSites, "-")
)

###########
#  Step 4 #
###########

#######
#  a  #
#######

# the center is a PureCLIP site
bindingSites <- bindingSites[(bindingSites - 2) %in% PureClipSites]

#######
#  b  #
#######

# Separate + and - strand binding sites
bindingSitesPlus <- bindingSites[strand(bindingSites) == "+"]
bindingSitesMinus <- bindingSites[strand(bindingSites) == "-"]

# Get the crosslink events in the binding sites as matrix (rows = binding sites)
crosslinkSitesInBindingSitesPlus <- as.matrix(crosslinkSites[["+"]][bindingSitesPlus])
crosslinkSitesInBindingSitesMinus <- as.matrix(crosslinkSites[["-"]][bindingSitesMinus])

# Determine the maximum number of crosslink events in each binding site
posWithMaxCrosslinksPlus <- apply(crosslinkSitesInBindingSitesPlus, 1, max)
posWithMaxCrosslinksMinus <- apply(crosslinkSitesInBindingSitesMinus, 1, max)

# Keep only binding sites where the center (position 3) has the maximum
#   number of crosslink events
bindingSitesPlus <- bindingSitesPlus[posWithMaxCrosslinksPlus == crosslinkSitesInBindingSitesPlus[,3]]
bindingSitesMinus <- bindingSitesMinus[posWithMaxCrosslinksMinus == crosslinkSitesInBindingSitesMinus[,3]]

# Re-combine the strands
bindingSites <- c(bindingSitesPlus, bindingSitesMinus)

# Remove unneccesary objects
rm(bindingSitesPlus, bindingSitesMinus,
   crosslinkSitesInBindingSitesPlus, crosslinkSitesInBindingSitesMinus,
   posWithMaxCrosslinksPlus, posWithMaxCrosslinksMinus)

#######
#  c  #
#######

# Separate + and - strand binding sites
bindingSitesPlus <- bindingSites[strand(bindingSites) == "+"]
bindingSitesMinus <- bindingSites[strand(bindingSites) == "-"]

# Get the crosslink events in the binding sites as matrix (rows = binding sites)
crosslinkSitesInBindingSitesPlus <- as.matrix(crosslinkSites[["+"]][bindingSitesPlus])
crosslinkSitesInBindingSitesMinus <- as.matrix(crosslinkSites[["-"]][bindingSitesMinus])

# Keep only binding sites that have at leat 3 (>2) positions with crosslink 
#   events
bindingSitesPlus <- bindingSitesPlus[apply((crosslinkSitesInBindingSitesPlus > 0),1,sum) > 2]
bindingSitesMinus <- bindingSitesMinus[apply((crosslinkSitesInBindingSitesMinus > 0),1,sum) > 2]

# Re-combine the strands
bindingSites <- c(bindingSitesPlus, bindingSitesMinus)

# Add the PureCLIP score of the binding site center 
bindingSites$score <- PureClipSites$score[match(bindingSites-2, PureClipSites)]

# Add seqlengths
seqlevels(bindingSites) <- paste0("chr", c(1:22,"X","Y"))
seqlengths(bindingSites) <- seqlengths(Hsapiens)[1:24]

seqlevels(crosslinkSites$`+`) <- paste0("chr", c(1:22,"X","Y"))
seqlengths(crosslinkSites$`+`) <- seqlengths(Hsapiens)[1:24]

seqlevels(crosslinkSites$`-`) <- paste0("chr", c(1:22,"X","Y"))
seqlengths(crosslinkSites$`-`) <- seqlengths(Hsapiens)[1:24]


# Store the binding sites in a RDS-File
saveRDS(bindingSites, paste0(iClipDir, "/rds_files/HNRNPH_binding_sites.rds"))

# Also store the crosslinkSites in a RDS-File
saveRDS(crosslinkSites, paste0(iClipDir,"/rds_files/HNRNPH_crosslink_sites.rds"))
