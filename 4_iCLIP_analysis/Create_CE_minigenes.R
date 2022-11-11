###########################################################################
#                                                                         #
# In this R script minigenes for the CE events with high-quality fittings #
#   and the non-regulated CE events are created. Minigenes are made up by #
#   They are made up the three involved exons.                            #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

majiqDir <- paste0(projectDir,"/2_MAJIQ_AS_analysis")
drcDir <- paste0(projectDir, "/3_Dose_response_curve_fitting")
iClipDir <- paste0(projectDir, "/4_iCLIP_analysis")

###############################
#          Load data          #
###############################

# Load the dose-response curve fittings of the CE events. The drc object
#   provides information like the fitted parameters and other information
fittedCEs <- readRDS(paste0(drcDir,"/rds_files/fittedCEs.rds"))

# Load the data.frame with regulated CE events and the additional information
#   and subset to the events with a high-quality fitting (pseudoR2 > 0.75).
#   This data.frame has one row per CE event and provides the reference junction
#   that was determined for each CE event.
regulatedCEs <- readRDS(paste0(majiqDir, "/rds_files/regulatedCEs.rds")) %>%
    dplyr::filter(event_id %in% names(fittedCEs))

# Load the orignial MAJIQ Modulizer output, extract the data.frame of the
#   CE events and subset to those with a high-qualitity fitting. This data.frame
#   has four rows per event and provides the coordinates of the three exons.
allCEs <- readRDS(paste0(majiqDir, "/rds_files/MAJIQ_binaryEvents.rds"))
allCEs <- allCEs$cassette %>%
    dplyr::filter(event_id %in% names(fittedCEs))

##################################
#          Computations          #
##################################

# The split() creates for each CE event a 4 row data.frame as list entry
regulatedMiniGenes <- lapply(split(allCEs, allCEs$event_id), function(df){
    
    eventID <- df$event_id[1]
    
    # Select columns of interest. The distinct removes one redundant row (C1_A
    #   and C2_A produce identical rows after the select()). In addition, the
    #   "spliced_with_coord" column is split into start and end columns
    df <- df %>%
        dplyr::select(module_id, gene_id, gene_name, event_id, seqid, strand,
                      spliced_with, spliced_with_coord) %>%
        distinct %>%
        dplyr::rename(chr=seqid, exon=spliced_with) %>%
        extract(spliced_with_coord,
                into=c("start", "end"), "([[:digit:]]+)-([[:digit:]]+)")
    
    # Here a lot of additional information is added.
    df$refJx <- regulatedCEs %>%
        dplyr::filter(event_id == eventID) %>%
        pull(junction_name)
    df$nH <- fittedCEs[[eventID]]$nH
    df$min <- fittedCEs[[eventID]]$min
    df$max <- fittedCEs[[eventID]]$max
    df$ec50 <- fittedCEs[[eventID]]$ec50
    df$pseudoR2 <- fittedCEs[[eventID]]$pseudoR2
    df$hillCat <- fittedCEs[[eventID]]$hillCat
    
    # Also the dPSI quantifications are added, which are extracted from
    #   the origData slot of the fitted model
    df <- cbind(df,
                fittedCEs[[eventID]]$origData %>%
                    tibble::rownames_to_column() %>%
                    dplyr::select(c(1,3)) %>%
                    pivot_wider(., names_from = "rowname", values_from="dPSI"))
    
    # In rare cases the start of an exon might be unknown. In this case, MAJIQ
    #   sets it to -1 in the "spliced_with_coord" column, which leads to a 
    #   a value of 1 in the start column. This case is filtered out. 
    df <- df %>% dplyr::filter(start != 1)
    
    #In rare cases the end of an exon might be unknown. This cases is filtered
    #   out.
    df <- df %>% dplyr::filter((!is.na(start)) & (!is.na(end)))
    
    # The events where one of the two filters above apply have less than
    #   3 rows in the data.frame. These are turned into an empty GRanges object.
    #   The data.frames of the other events can be turned into a GRanges object
    #   with 3 ranges (one for each exon).
    if(nrow(df) == 3){
        gr <- df %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    } else {
        gr <- GRanges()
    }
    
    return(gr)
}) %>% as(., "GRangesList")

# Remove the CE events that do not have 3 exons 
regulatedMiniGenes <- regulatedMiniGenes[lengths(regulatedMiniGenes)==3]

#Add seqlengths
regulatedMiniGenes <- lapply(regulatedMiniGenes, function(gr){
    seqlevels(gr) <- paste0("chr", c(1:22,"X","Y"))
    seqlengths(gr) <- seqlengths(Hsapiens)[1:24]
    return(gr)
})

# The regulated minigenes are stored in a RDS-File.
saveRDS(regulatedMiniGenes, paste0(iClipDir,"/rds_files/regulatedMiniGenes.rds"))



###############################
#          Load data          #
###############################

# Load the data.frame with non-regulated CE events and the additional information
#   This data.frame has one row per CE event and provides the reference junction
#   that was determined for each CE event.
nonregulatedCEs <- readRDS(paste0(majiqDir, "/rds_files/nonregulatedCEs.rds"))

# Load the orignial MAJIQ Modulizer output, extract the data.frame of the
#   CE events and subset to the non-regulated CE events. This data.frame
#   has four rows per event and provides the coordinates of the three exons.
allCEs <- readRDS(paste0(majiqDir, "/rds_files/MAJIQ_binaryEvents.rds"))
allCEs <- allCEs$cassette %>%
    dplyr::filter(event_id %in% nonregulatedCEs$event_id)

# The split() creates for each CE event a 4 row data.frame as list entry
nonregulatedMiniGenes <- lapply(split(allCEs, allCEs$event_id),function(df){
    
    eventID <- df$event_id[1]
    
    # Select columns of interest. The distinct removes one redundant row (C1_A
    #   and C2_A produce identical rows after the select()). In addition, the
    #   "spliced_with_coord" column is split into start and end columns
    df <- df %>%
        dplyr::select(module_id, gene_id, gene_name, event_id, seqid, strand,
                      spliced_with, spliced_with_coord) %>%
        distinct %>%
        dplyr::rename(chr=seqid, exon=spliced_with) %>%
        extract(spliced_with_coord,
                into=c("start", "end"), "([[:digit:]]+)-([[:digit:]]+)")
    df$refJx <- nonregulatedCEs %>%
        dplyr::filter(event_id == eventID) %>%
        pull(junction_name)
    
    # In rare cases the start of an exon might be unknown. In this case, MAJIQ
    #   sets it to -1 in the "spliced_with_coord" column, which leads to a 
    #   a value of 1 in the start column. This case is filtered out. 
    df <- df %>% dplyr::filter(start != 1)
    
    #In rare cases the end of an exon might be unknown. This cases is filtered
    #   out.
    df <- df %>% dplyr::filter((!is.na(start)) & (!is.na(end)))
    
    # The events where one of the two filters above apply have less than
    #   3 rows in the data.frame. These are turned into an empty GRanges object.
    #   The data.frames of the other events can be turned into a GRanges object
    #   with 3 ranges (one for each exon).
    if(nrow(df) == 3){
        gr <- df %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    } else {
        gr <- GRanges()
    }
    
    return(gr)
}) %>% as(., "GRangesList")

# Remove the CE events that do not have 3 exons 
nonregulatedMiniGenes <- nonregulatedMiniGenes[lengths(nonregulatedMiniGenes)==3]

#Add seqlengths
nonregulatedMiniGenes <- lapply(nonregulatedMiniGenes, function(gr){
    seqlevels(gr) <- paste0("chr", c(1:22,"X","Y"))
    seqlengths(gr) <- seqlengths(Hsapiens)[1:24]
    return(gr)
})

# The regulated minigenes are stored in a RDS-File.
saveRDS(nonregulatedMiniGenes, paste0(iClipDir,"/rds_files/nonregulatedMiniGenes.rds"))