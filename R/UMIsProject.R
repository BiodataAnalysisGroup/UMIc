


########## Loading libraries ########## 
library(optparse)
library(getopt)

rm(list = ls())

# Source files relative to the current script
scriptPath = get_Rscript_filename()
functions = paste(dirname(scriptPath), "functions.R", sep="/")
workflows = paste(dirname(scriptPath), "casesWorkflows.R", sep="/")
source(functions)
source(workflows)

########## Inputs ##########

parser <- OptionParser()
parser <- add_option(parser, "--pairedData", action="store_true", default=FALSE,
                     help="The reads are paired-end")

parser <- add_option(parser, "--UMIlocation", default="R1 & R2",
                     help="UMI location ('R1' or 'R1 & R2')")

parser <- add_option(parser, "--UMIlength", default=12,
                     help="Length of the UMI")

parser <- add_option(parser, "--sequenceLength", default=150,
                     help="Length of the sequnce")

parser <- add_option(parser, "--countsCutoff", default=1,
                     help="Minimum read count per UMI, for initial data cleaning")

parser <- add_option(parser, "--UMIdistance", default=4,
                     help="Maximum UMI distance for UMI merging")

parser <- add_option(parser, "--sequenceDistance", default=100,
                     help="Maximum sequence distance for UMI correction")

parser <- add_option(parser, "--inputsFolder",
                     help="Folder with FastQ files to parse")

parser <- add_option(parser, "--outputsFolder",
                     help="Folder to write the output to")


########## Parse the command line options ##########
args <- parse_args(parser)
pairedData = args$pairedData
UMIlocation = args$UMIlocation
UMIlength = args$UMIlength
sequenceLength = args$sequenceLength
countsCutoff = args$countsCutoff
UMIdistance = args$UMIdistance
sequenceDistance = args$sequenceDistance
inputsFolder = args$inputsFolder
outputsFolder = args$outputsFolder


########## Loading libraries (slow) ##########
library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(stringdist)
library(pryr)

########## Run the appropriate scenario ##########
if (pairedData & UMIlocation == "R1"){   #case 1 -- paired data and UMI only in Read1
  
  inputFiles <- list.files(inputsFolder, pattern = "fastq") 
  
  while(length(inputFiles) > 0){
    
    file1 <- inputFiles[1]
    
    commonPart <- as.character(str_split(file1,"R1", simplify = T))
    commonPart <- commonPart[length(commonPart)]
    
    file2 <- str_replace(file1,paste0("R1",commonPart),paste0("R2",commonPart))
    
    filepath1 <- paste0(inputsFolder,"/",file1)
    filepath2 <- paste0(inputsFolder,"/",file2)
    
    pairedR1(filepath1, 
             filepath2, 
             outputsFolder, 
             UMIlength, 
             UMIdistance, 
             sequenceLength, 
             sequenceDistance, 
             countsCutoff)
    
    inputFiles <- inputFiles[str_detect(inputFiles,paste0(file1,"|",file2), negate = T)]
  }
 
  
} else if (pairedData & UMIlocation == "R1 & R2"){   #case 2 -- paired data and UMI in Read1 and Read2
  
  inputFiles <- list.files(inputsFolder, pattern = "fastq") 
  
  while(length(inputFiles) > 0){
    
    file1 <- inputFiles[1]
    
    commonPart <- as.character(str_split(file1, "R1", simplify = T))
    commonPart <- commonPart[length(commonPart)]
    
    file2 <- str_replace(file1, 
                         paste0("R1", commonPart), 
                         paste0("R2", commonPart))
    
    filepath1 <- paste0(inputsFolder, "/", file1)
    filepath2 <- paste0(inputsFolder, "/", file2)
    
    cat(c("Files:", file1, file2, "\n"))

    pairedR1R2(filepath1, 
               filepath2, 
               outputsFolder, 
               UMIlength, 
               UMIdistance, 
               sequenceLength, 
               sequenceDistance, 
               countsCutoff)
    
    
    inputFiles <- inputFiles[str_detect(inputFiles, paste0(file1, "|", file2), negate = T)]
  }

  
} else if (!pairedData){  #case 3 -- single data
  
  inputFiles <- list.files(inputsFolder, pattern = "fastq") 
  
  while(length(inputFiles) > 0){
    
    file1 <- inputFiles[1]
    
    filepath1 <- paste0(inputsFolder,"/",file1)

    single(
        filepath1, 
        outputsFolder, 
        UMIlength, 
        UMIdistance, 
        sequenceLength, 
        sequenceDistance, 
        countsCutoff
    )
    
    inputFiles <- inputFiles[str_detect(inputFiles,file1, negate = T)]
  }
 
}  

