###########################################################################
#                                                                         #
# In this R script a GRanges object for junctions of interest is          #
#   created. The object will be used to subset the                        #
#   RangedSummarizedExperiment object containing the counts of all        #
#   junctions quantified in BRCA samples. There are always two junctions  #
#   per CE and ALE event.                                                 #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(GenomicRanges)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

majiqDir <- paste0(projectDir, "/2_MAJIQ_AS_analysis")
drcDir <- paste0(projectDir, "/3_Dose_response_curve_fitting")
tcgaDir <- paste0(projectDir, "/8_TCGA_analysis")

###########################################
#          Load and adjust data          #
##########################################

MAJIQ_binaryEvents <- readRDS(paste0(majiqDir,"/rds_files/MAJIQ_binaryEvents.rds"))

regulatedCEs <- readRDS(paste0(majiqDir,"/rds_files/regulatedCEs.rds"))
regulatedALEs <- readRDS(paste0(majiqDir,"/rds_files/regulatedALEs.rds"))
nonRegulatedCEs <- readRDS(paste0(majiqDir,"/rds_files/nonregulatedCEs.rds"))

fittedCEs <- readRDS(paste0(drcDir,"/rds_files/fittedCEs.rds"))
fittedALEs <- readRDS(paste0(drcDir,"/rds_files/fittedALEs.rds"))

# Subset the regulated CE events to those with a good fitting (pseudo R2 >= 0.75)
regulatedCEs <- regulatedCEs %>% dplyr::filter(event_id %in% names(fittedCEs))
regulatedALEs <- regulatedALEs %>% dplyr::filter(event_id %in% names(fittedALEs))

#########################################################
#          Create the junction GRanges objects          #
#########################################################

# Create the GRanges object for the inclusion (IJ) and skipping junctions (SJ)
jxGRangesCEs <- lapply(1:nrow(regulatedCEs), function(i){

    # General information
    chr <- regulatedCEs$seqid[i]
    strand <- regulatedCEs$strand[i]
    #Meta information
    eventID <- regulatedCEs$event_id[i]
    geneName <- regulatedCEs$gene_name[i]
    geneID <- regulatedCEs$gene_id[i]
    hillCat <- fittedCEs[[eventID]]$hillCat
    # IJ - Information is stored in the regulatedCEs data.frame
    start <- regulatedCEs$junction_coord[i] %>%
        strsplit("-", fixed=T) %>%
        sapply(., "[[", 1)
    end <- regulatedCEs$junction_coord[i] %>%
        strsplit("-", fixed=T) %>%
        sapply(., "[[", 2)
    jxType <- "IJ"
    jxName <- regulatedCEs$junction_name[i]
    # SJ - Information needs to be extracted from MAJIQ_binaryEvents
    jxName <- c(jxName,
                ifelse(regulatedCEs$junction_name[i] == "C1_A", "C1_C2", "C2_C1"))
    start <- c(start,
               MAJIQ_binaryEvents$cassette %>%
                   dplyr::filter(event_id == eventID & junction_name == jxName[2]) %>%
                   pull(junction_coord) %>%
                   strsplit("-", fixed=T) %>%
                   sapply(., "[[", 1))
    end <- c(end,
             MAJIQ_binaryEvents$cassette %>%
                 dplyr::filter(event_id == eventID & junction_name == jxName[2]) %>%
                 pull(junction_coord) %>%
                 strsplit("-", fixed=T) %>% 
                 sapply(., "[[", 2))
    jxType <- c(jxType, "SJ")
    
    df <- data.frame(chr, strand, start, end, jxType, jxName, eventID,
                     geneName, geneID, hillCat)
    
    return(df)
}) %>% bind_rows() %>% makeGRangesFromDataFrame(., keep.extra.columns = T)

jxGRangesALEs <- lapply(1:nrow(regulatedALEs), function(i){
    
    # General information
    chr <- regulatedALEs$seqid[i]
    strand <- regulatedALEs$strand[i]
    #Meta information
    eventID <- regulatedALEs$event_id[i]
    geneName <- regulatedALEs$gene_name[i]
    geneID <- regulatedALEs$gene_id[i]
    hillCat <- fittedALEs[[eventID]]$hillCat
    # IJ - Information is stored in the regulatedALEs data.frame
    start <- regulatedALEs$junction_coord[i] %>%
        strsplit("-", fixed=T) %>%
        sapply(., "[[", 1)
    end <- regulatedALEs$junction_coord[i] %>%
        strsplit("-", fixed=T) %>%
        sapply(., "[[", 2)
    jxType <- "IJ"
    jxName <- regulatedALEs$junction_name[i]
    # SJ - Information needs to be extracted from MAJIQ_binaryEvents
    jxName <- c(jxName, "Distal")
    start <- c(start,
               MAJIQ_binaryEvents$alternate_last_exon %>%
                   dplyr::filter(event_id == eventID & junction_name == jxName[2]) %>%
                   pull(junction_coord) %>%
                   strsplit("-", fixed=T) %>%
                   sapply(., "[[", 1))
    end <- c(end,
             MAJIQ_binaryEvents$alternate_last_exon %>%
                 dplyr::filter(event_id == eventID & junction_name == jxName[2]) %>%
                 pull(junction_coord) %>%
                 strsplit("-", fixed=T) %>%
                 sapply(., "[[", 2))
    jxType <- c(jxType, "SJ")
    
    df <- data.frame(chr, strand, start, end, jxType, jxName, eventID,
                     geneName, geneID, hillCat)
    
    return(df)
}) %>% bind_rows() %>% makeGRangesFromDataFrame(., keep.extra.columns = T)

jxGRangesNonRegCEs <- lapply(1:nrow(nonRegulatedCEs), function(i){
    
    # General information
    chr <- nonRegulatedCEs$seqid[i]
    strand <- nonRegulatedCEs$strand[i]
    #Meta information
    eventID <- nonRegulatedCEs$event_id[i]
    geneName <- nonRegulatedCEs$gene_name[i]
    geneID <- nonRegulatedCEs$gene_id[i]
    hillCat <- NA
    # IJ - Information is stored in the nonRegulatedCEs data.frame
    start <- nonRegulatedCEs$junction_coord[i] %>%
        strsplit("-", fixed=T) %>%
        sapply(., "[[", 1)
    end <- nonRegulatedCEs$junction_coord[i] %>%
        strsplit("-", fixed=T) %>%
        sapply(., "[[", 2)
    jxType <- "IJ"
    jxName <- nonRegulatedCEs$junction_name[i]
    # SJ - Information needs to be extracted from MAJIQ_binaryEvents
    jxName <- c(jxName,
                ifelse(nonRegulatedCEs$junction_name[i] == "C1_A", "C1_C2", "C2_C1"))
    start <- c(start,
               MAJIQ_binaryEvents$cassette %>%
                   dplyr::filter(event_id == eventID & junction_name == jxName[2]) %>%
                   pull(junction_coord) %>%
                   strsplit("-", fixed=T) %>%
                   sapply(., "[[", 1))
    end <- c(end, MAJIQ_binaryEvents$cassette %>%
                 dplyr::filter(event_id == eventID & junction_name == jxName[2]) %>%
                 pull(junction_coord) %>%
                 strsplit("-", fixed=T) %>%
                 sapply(., "[[", 2))
    jxType <- c(jxType, "SJ")
    
    df <- data.frame(chr, strand, start, end, jxType, jxName, eventID, geneName, geneID, hillCat)
    
    return(df)
}) %>% bind_rows() %>% makeGRangesFromDataFrame(., keep.extra.columns = T)

# There is a CE event, where the skipping leads to a NMD isoform. We also
#   want to monitor this event in the BRCA samples
jxGrangesHNRNPH1nmdIsoform <-GRanges(
    seqnames = "chr5",
    strand = "-",
    ranges=IRanges(start=c(179619407, 179618323, 179618323),end=c(179620892, 179619269, 179620892)),
    jxType=c("IJ1", "IJ2", "SJ"),
    jxName <- NA,
    eventID = "HNRNPH1nmdIsoform",
    geneName="HNRNPH1",
    geneID="ENSG00000169045",
    hillCat <- NA)

# Add the eventType to each GRanges object
jxGRangesCEs$eventType <- "regulatedCE"
jxGRangesALEs$eventType <- "regulatedALE"
jxGRangesNonRegCEs$eventType <- "nonRegulatedCE"
jxGrangesHNRNPH1nmdIsoform$eventType <- "HNRNPH1nmdIsoform"

# Combine the four objects into a single one
jxGRanges <- c(jxGRangesCEs, jxGRangesALEs,
               jxGRangesNonRegCEs, jxGrangesHNRNPH1nmdIsoform)

# Save the object as RDS-File
saveRDS(jxGRanges, paste0(tcgaDir,"/rds_files/jxGRanges.rds"))
