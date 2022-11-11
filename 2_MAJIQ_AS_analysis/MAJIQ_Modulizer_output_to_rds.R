###########################################################################
#                                                                         #
# In this R script the 14 MAJIQ Modulizer TSV-Files are transfered into   #
#   a list of data.frames (one for each TSV-File), which is stored as     #
#   RDS-File.                                                             #
#                                                                         #
###########################################################################

###############################
#          Libraries          #
###############################

library(tidyverse)

#################################
#          Directories          #
#################################

projectDir <- "/Users/mariokeller/projects/HNRNPH_project/Tretow_et_al_2023"

datarDir <- paste0(projectDir, "/data")
majiqDir <- paste0(projectDir,"/2_MAJIQ_AS_analysis")


# Create the list by reading TSV-Files as data.frames
MAJIQ_binaryEvents <- list(
    cassette=read.table(paste0(datarDir,"/MAJIQ_Modulizer/cassette.tsv"),
                        header=T, sep="\t"),
    alternative_intron=read.table(paste0(datarDir,"/MAJIQ_Modulizer/alternative_intron.tsv"),
                        header=T, sep="\t"),
    alt5prime=read.table(paste0(datarDir,"/MAJIQ_Modulizer/alt5prime.tsv"),
                        header=T, sep="\t"),
    alt3prime=read.table(paste0(datarDir,"/MAJIQ_Modulizer/alt3prime.tsv"),
                        header=T, sep="\t"),
    alt3and5prime=read.table(paste0(datarDir,"/MAJIQ_Modulizer/alt3and5prime.tsv"),
                        header=T, sep="\t"),
    alternate_first_exon=read.table(paste0(datarDir,"/MAJIQ_Modulizer/alternate_first_exon.tsv"),
                        header=T, sep="\t"),
    alternate_last_exon=read.table(paste0(datarDir,"/MAJIQ_Modulizer/alternate_last_exon.tsv"),
                        header=T, sep="\t"),
    p_alt5prime=read.table(paste0(datarDir,"/MAJIQ_Modulizer/p_alt5prime.tsv"),
                        header=T, sep="\t"),
    p_alt3prime=read.table(paste0(datarDir,"/MAJIQ_Modulizer/p_alt3prime.tsv"),
                        header=T, sep="\t"),
    p_alternate_first_exon=read.table(paste0(datarDir,"/MAJIQ_Modulizer/p_alternate_first_exon.tsv"),
                        header=T, sep="\t"),
    p_alternate_last_exon=read.table(paste0(datarDir,"/MAJIQ_Modulizer/p_alternate_last_exon.tsv"),
                        header=T, sep="\t"),
    mutually_exclusive=read.table(paste0(datarDir,"/MAJIQ_Modulizer/mutually_exclusive.tsv"),
                        header=T, sep="\t"),
    tandem_cassette=read.table(paste0(datarDir,"/MAJIQ_Modulizer/tandem_cassette.tsv"),
                        header=T, sep="\t"),
    multi_exon_spanning=read.table(paste0(datarDir,"/MAJIQ_Modulizer/multi_exon_spanning.tsv"),
                        header=T, sep="\t")
    )

# Store the list as RDS-File
saveRDS(MAJIQ_binaryEvents, paste0(majiqDir, "/rds_files/MAJIQ_binaryEvents.rds"))



