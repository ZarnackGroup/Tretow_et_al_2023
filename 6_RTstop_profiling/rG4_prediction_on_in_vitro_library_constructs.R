###############################################################
#                                                             #
# In this R script RNA G-quadruplexes (rG4s) are predicted in #
#   the constructs of the in vitro library                    #
#                                                             #
###############################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(pqsfinder)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(janitor)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

RTstopDir <- paste0(projectDir, "/6_RTstop_profiling")

###############################
#          Load data          #
###############################

constructInformation <- readRDS(paste0(RTstopDir,
                                       "/rds_files/constructInformation.rds"))
constructInformation <- clean_names(constructInformation, case = "lower_camel")

#################################
#          Adjust data          #
#################################

rss <- DNAStringSet(constructInformation$seq %>%
                        setNames(., constructInformation$constructId)) %>%
    RNAStringSet()
    
##################################
#          Predict rG4s          #
##################################

# For the prediction rG4s with a maximum length of 30, loop length of 1 to 7
#   and no defects (mismatches or bulges) were taken into account. In addition,
#   only rG4s with a minimum score of at least 25 were considered.

pqsfinderResults_noOverlaps <- lapply(rss, function(rs){
  pqsfinder(rs,
            strand="+", 
            max_len = 30L, min_score=25L,
            loop_min_len = 1L, loop_max_len = 7L,
            max_defects=0 # no bulges or mismatches
  )
})

# The predicted rG4s are stored in RDS-Files
saveRDS(pqsfinderResults_noOverlaps, paste0(RTstopDir, "/rds_files/pqsfinderResults_noOverlaps.rds"))

