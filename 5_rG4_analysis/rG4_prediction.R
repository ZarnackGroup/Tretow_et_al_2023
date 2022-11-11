##########################################################################
#                                                                        #
# In this R script RNA G-quadruplexes (rG4s) are predicted in            #
#   regulated and non-regulated CE events and around all HNRNPH          #
#   binding sites. This script runs for multiple hours due to the        #
#   prediction around the > 400.000 HNRNPH binding sites                 #
#                                                                        #
##########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)
library(pqsfinder)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

iClipDir <- paste0(projectDir, "/4_iCLIP_analysis")
rG4Dir <- paste0(projectDir, "/5_rG4_analysis")

###############################
#          Load data          #
###############################

regulatedMiniGenes <- readRDS(paste0(iClipDir,
                                     "/rds_files/regulatedMiniGenes.rds"))
nonRegulatedMiniGenes <- readRDS(paste0(iClipDir,
                                        "/rds_files/nonregulatedMiniGenes.rds"))

bindingSites <- readRDS(paste0(iClipDir,
                               "/rds_files/HNRNPH_binding_sites.rds"))

#################################
#          Adjust data          #
#################################

# Reduce regulated and non-regulated minigenes to a single range 
#   starting at the start of the leftmost exon and ending at the end of the
#   rightmost exon
regulatedMiniGenes <- lapply(regulatedMiniGenes, function(gr){
    event_id <- gr$event_id %>% unique
    hillCat <- gr$hillCat %>% unique
    # range() reduces the three ranges into one range
    gr <- range(gr)
    gr$event_id <- event_id
    gr$hillCat <- hillCat
    return(gr)
}) %>% as(., "GRangesList") %>% unlist

nonRegulatedMiniGenes <- lapply(nonRegulatedMiniGenes, function(gr){
    event_id <- gr$event_id %>% unique
    # range() reduces the three ranges into one range
    gr <- range(gr)
    gr$event_id <- event_id
    gr$hillCat <- "nonReg"
    return(gr)
}) %>% as(., "GRangesList") %>% unlist

# Create a window of +/-50 nt around the binding site center
bindingSites <- bindingSites+48

# Get the sequences
regulatedMiniGenesRSS <- RNAStringSet(getSeq(Hsapiens, regulatedMiniGenes))
nonRegulatedMiniGenesRSS <- RNAStringSet(getSeq(Hsapiens, nonRegulatedMiniGenes))
bindingSitesRSS <- RNAStringSet(getSeq(Hsapiens, bindingSites))

##################################
#          Predict rG4s          #
##################################

# For the prediction rG4s with a maximum length of 30, loop length of 1 to 7
#   and no defects (mismatches or bulges) were taken into account. In addition,
#   only rG4s with a minimum score of at least 25 were considered.

regulatedMiniGenesG4s <- lapply(regulatedMiniGenesRSS, function(rss){
  pqsfinder(rss,
            strand="+", 
            max_len = 30L, min_score=25L,
            loop_min_len = 1L, loop_max_len = 7L,
            max_defects=0 # no bulges or mismatches
  )
})

nonRegulatedMiniGenesG4s <- lapply(nonRegulatedMiniGenesRSS, function(rss){
  pqsfinder(rss,
            strand="+", 
            max_len = 30L, min_score=25L,
            loop_min_len = 1L, loop_max_len = 7L,
            max_defects=0 # no bulges or mismatches
  )
})

bindingSitesG4s <- lapply(bindingSitesRSS, function(rss){
    pqsfinder(rss,
              strand="+", 
              max_len = 30L, min_score=25L,
              loop_min_len = 1L, loop_max_len = 7L,
              max_defects=0 # no bulges or mismatches
    )
})

#####################################################################
#          Turn rG4s coordinates into genomic coordinates           #
#####################################################################

regulatedMiniGenesG4s <- lapply(names(regulatedMiniGenesG4s), function(event_id){

    # Fetch the predicted G4s and turn them into a GRanges object
    G4gr <- regulatedMiniGenesG4s[[event_id]] %>% as(., "GRanges")
      
    # If there was no predicted G4 return an empty GRanges object
    if(length(G4gr)==0){
        return(GRanges())
    }
    
    # Extract the Input RNA sequence (RNAString class)
    RNAseq <- subject(regulatedMiniGenesG4s[[event_id]]) 
     
    # Excise the rG4 sequences from the RNA sequence
    G4seqs <- sapply(1:length(G4gr),function(i){
        RNAseq[start(G4gr)[i]:end(G4gr)[i]] %>% as.character
    })
    
    # Fetch the minigene to extract information of the genomic location
    MGgr <- regulatedMiniGenes[regulatedMiniGenes$event_id == event_id]
    strand(G4gr) <- strand(MGgr)
    seqlevels(G4gr) <- seqlevels(MGgr)
    seqnames(G4gr) <- seqnames(MGgr)
    
    # Place the rG4s to the correct genomic location
    if(strand(G4gr) %>% as.character == "+"){
        G4gr <- shift(G4gr, start(MGgr)-1)
    } else {
        G4gr <- shift(G4gr, - (width(MGgr)+1))
        new_start <- end(G4gr) * -1
        new_end <- start(G4gr) * -1
        end(G4gr) <- new_end
        start(G4gr) <- new_start
        G4gr <- shift(G4gr, start(MGgr)-1)
      }
    
    # Add the sequence of the rG4
    G4gr$sequence <- G4seqs
    return(G4gr)
}) %>% setNames(., regulatedMiniGenes$event_id)


nonRegulatedMiniGenesG4s <- lapply(names(nonRegulatedMiniGenesG4s), function(event_id){
    
    # Fetch the predicted G4s and turn them into a GRanges object
    G4gr <- nonRegulatedMiniGenesG4s[[event_id]] %>% as(., "GRanges")
    
    # If there was no predicted G4 return an empty GRanges object
    if(length(G4gr)==0){
        return(GRanges())
    }
  
    # Extract the Input RNA sequence (RNAString class)
    RNAseq <- subject(nonRegulatedMiniGenesG4s[[event_id]])
    
    # Excise the rG4 sequences from the RNA sequence
    G4seqs <- sapply(1:length(G4gr),function(i){
        RNAseq[start(G4gr)[i]:end(G4gr)[i]] %>% as.character
    })
    
    # Fetch the minigene to extract information of the genomic location
    MGgr <- nonRegulatedMiniGenes[nonRegulatedMiniGenes$event_id == event_id]
    strand(G4gr) <- strand(MGgr)
    seqlevels(G4gr) <- seqlevels(MGgr)
    seqnames(G4gr) <- seqnames(MGgr)
    
    # Place the rG4s to the correct genomic location
    if(strand(G4gr) %>% as.character == "+"){
        G4gr <- shift(G4gr, start(MGgr)-1)
    } else {
        G4gr <- shift(G4gr, -(width(MGgr)+1))
        new_start <- end(G4gr) * -1
        new_end <- start(G4gr) * -1
        end(G4gr) <- new_end
        start(G4gr) <- new_start
        G4gr <- shift(G4gr, start(MGgr)-1)
        }
    
    # Add the sequence of the rG4
    G4gr$sequence <- G4seqs
  return(G4gr)
}) %>% setNames(., nonRegulatedMiniGenes$event_id)


bindingSitesG4s <- lapply(1:length(bindingSites), function(BS){
    
    # Fetch the predicted G4s and turn them into a GRanges object
    G4gr <- bindingSitesG4s[[BS]] %>% as(., "GRanges")
    
    # If there was no predicted G4 return an empty GRanges object
    if(length(G4gr)==0){
        return(GRanges())
    }
    
    # Extract the Input RNA sequence (RNAString class)
    RNAseq <- subject(bindingSitesG4s[[BS]]) #minigene sequence
    
    # Excise the rG4 sequences from the RNA sequence
    G4seqs <- sapply(1:length(G4gr),function(i){
        RNAseq[start(G4gr)[i]:end(G4gr)[i]] %>% as.character
    })
    
    # Fetch the bindg site window to extract information of the genomic location
    BSgr <- bindingSites[BS]
    strand(G4gr) <- strand(BSgr)
    seqlevels(G4gr) <- seqlevels(BSgr)
    seqnames(G4gr) <-seqnames(BSgr)
    
    # Place the rG4s to the correct genomic location
    if(strand(G4gr) %>% as.character == "+"){
        G4gr <- shift(G4gr, start(BSgr)-1)
    } else {
        G4gr <- shift(G4gr, -(width(BSgr)+1))
        new_start <- end(G4gr) * -1
        new_end <- start(G4gr) * -1
        end(G4gr) <- new_end
        start(G4gr) <- new_start
        G4gr <- shift(G4gr, start(BSgr)-1)
        
    }
    
    # Add the sequence of the rG4
    G4gr$sequence <- G4seqs
    return(G4gr)
})

# Unlist
bindingSitesG4s <- bindingSitesG4s %>% as(., "GRangesList") %>% unlist

# Remove duplicated rG4s, which arise from close by binding sites that have
#   overlapping windows in which the same rG4 can be predicted
bindingSitesG4s <- unique(bindingSitesG4s)

#Add seqlengths
regulatedMiniGenesG4s <- lapply(regulatedMiniGenesG4s, function(gr){
    seqlevels(gr) <- paste0("chr", c(1:22,"X","Y"))
    seqlengths(gr) <- seqlengths(Hsapiens)[1:24]
    return(gr)
})

nonRegulatedMiniGenesG4s <- lapply(nonRegulatedMiniGenesG4s, function(gr){
    seqlevels(gr) <- paste0("chr", c(1:22,"X","Y"))
    seqlengths(gr) <- seqlengths(Hsapiens)[1:24]
    return(gr)
})

seqlevels(bindingSitesG4s) <- paste0("chr", c(1:22,"X","Y"))
seqlengths(bindingSitesG4s) <- seqlengths(Hsapiens)[1:24]


# The predicted rG4s are stored in RDS-Files
saveRDS(regulatedMiniGenesG4s, paste0(rG4Dir,
                                      "/rds_files/regulatedMiniGenesG4s.rds"))
saveRDS(nonRegulatedMiniGenesG4s, paste0(rG4Dir,
                                         "/rds_files/nonRegulatedMiniGenesG4s.rds"))
saveRDS(bindingSitesG4s, paste0(rG4Dir,
                                "/rds_files/bindingSitesG4s.rds"))
