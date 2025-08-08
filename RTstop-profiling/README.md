# RTstop-profiling

NEWS

Update 14-02-2022 (by Mirko)
* Added Construct classification in 03_PeakDifferencesFinal folder
  * Important new files that are located in the /data/ folder are:
    * constructInfo2.rds -> basically the orignal constructInfo table with the classififcation status added (use this to extend shiny app)
    * peaks_expressed.rds -> GRangesObject of all peaks defined on all constructs, with the respective classification of the construct added
    * G4Peaks.rds -> basically a subset on the peaks_expressed GRangesObject, selection for "Up" and "Mixed" group and only significant Peaks. Additionally the KCl and NaCl ratios (read-through/ peak signal) is added
    * g4peaklist.xlsx -> the same as G4peaks object, but as excel
  * Please don't kill me for the stupid naming of files and folder xD
* Added Construct classification based on signal shape in 04_PeakClassification folder
  * we will probably never use this
  * just sits here
* Added 05_ConstructClassification folder
  * This is a placeholder where I will set up the final Region/ Construct/ Peak classification and filtering scheme based on our last meeting once I'm back
  * This is empty atm
