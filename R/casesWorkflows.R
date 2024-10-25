pairedR1 <- function(filepath1, 
                     filepath2, 
                     outputsFolder, 
                     UMIlength, 
                     UMIdistance, 
                     sequenceLength, 
                     sequenceDistance,
                     countsCutoff){
  
  dir.create(outputsFolder)
  file.create(paste0(outputsFolder, "/extra_info.txt")) 
  memory <- c() 
  memory[1] <-  mem_used()
  #read input files
  reads1 <- readFastq(filepath1)
  reads2 <- readFastq(filepath2)
  
  start.time <- Sys.time()
  # reads1 = reads1[1:300,]
  # reads2 = reads2[1:300,]
  
  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
  
  #separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength+1, sequenceLength)
  
  #File 2
  seq2 <- as.data.table(sread(reads2))
  ids2 <- as.data.table(reads2@id)
  
  full2 <- cbind(seq2,ids2)
  names(full2) <- c("read", "id")
   
  #keep intermediate information 
  partIDs <- as.data.table(cbind(UMI = full$UMI, ID1 = full$id, ID2 = full2$id))
  partIDs <- partIDs[!duplicated(UMI),]
  partIDs <- partIDs[order(partIDs$UMI, decreasing = TRUE), ]
  
  intermediate.table <- full[,.(count = .N),by=UMI,]
  intermediate.table <- intermediate.table[order(intermediate.table$UMI, decreasing = TRUE), ]
  intermediate.table$ID1 <- partIDs$ID1
  intermediate.table$ID2 <- partIDs$ID2
  intermediate.table <- intermediate.table[which(intermediate.table$count >= countsCutoff),]
  
  rm(partIDs)
  memory[2] <- mem_used()
  intermediate.table = intermediate.table[order(intermediate.table$count, decreasing = TRUE), ]
  
  #data preparation
  #file 1
  full <- separate(full, id, c("id1", "id2"), " ", remove = T)
  full <- select(full, read, id1, UMI)
  colnames(full) <- c("read", "id", "UMI")
  
  quality <- as(quality(reads1), "matrix")
  quality = as.data.table(quality)
  quality = quality[,(UMIlength+1):sequenceLength]
  quality$id <- full$id
  
  full2<- separate(full2,id, c("id1", "id2")," ",remove = T)
  full2 <- select(full2,read, "id1")
  colnames(full2) <- c("read", "id")
  
  quality2 <- as(quality(reads2), "matrix")
  quality2 = as.data.table(quality2)
  quality2$id <- full2$id
  
  rm(ids, ids2, seq, seq2, reads1, reads2)
  
  
  #first consensus
  # result_mean = list()
  # 
  # for(i in c(1:nrow(intermediate.table))){
    
    result_mean = groupingPairedR1(intermediate.table, # $UMI[i], 
                                   # intermediate.table$count[i], 
                                   full, 
                                   quality, 
                                   full2, 
                                   quality2, 
                                   UMIlength)
    
  # }
  
  # result_mean <- bind_rows(result_mean)
  memory[3] <- mem_used()
  #UMI correction
  newUMIs <- UMIcorrectionPairedR1(intermediate.table,
                                   result_mean,
                                   sequenceDistance, 
                                   UMIdistance,
                                   outputsFolder)
  
  
  rm(intermediate.table)
  memory[4] <- mem_used()
  #final consensus
  # consensus_mean = list()
  
  # for(i in newUMIs$UMI){
    
    consensus_mean = groupingFinalPairedR1(newUMIs, # i, 
                                           full, 
                                           quality, 
                                           full2, 
                                           quality2,
                                           result_mean, 
                                           UMIlength)
    
    memory[5] <- mem_used()
    
    consensus_mean = consensus_mean[order(consensus_mean$UMI), ]
    result_mean    = result_mean[order(result_mean$UMI), ]
    newUMIs$key = str_sub(newUMIs$UMI, 1, UMIlength)
    newUMIs        = newUMIs[order(newUMIs$key), ]
  # }
  
  # consensus_mean <- bind_rows(consensus_mean)

  #produce Outputs 

  end.time <- Sys.time()
  
  total.time.secs <- difftime(end.time,start.time,units = "secs")
  total.time.mins <- difftime(end.time,start.time,units = "mins")
  
  line = paste0("Total time in secs: ", total.time.secs)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Total time in mins: ", total.time.mins)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used at start: ", memory[1])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after data cleaning: ", memory[2])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after first consensus: ", memory[3])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after UMI merging: ", memory[4])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after second consensus: ", memory[5])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  #File1
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1),
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(newUMIs$ID1," ",consensus_mean$UMI)))


  fileSplit <- as.data.table(str_split(filepath1,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  part <- fileSplit[1]
  file.create(output)
  writeFastq(file, output, mode = "a")

  #File2
  file <- ShortReadQ(DNAStringSet(consensus_mean$read2),
                     FastqQuality(consensus_mean$quality2),
                     BStringSet(paste0(newUMIs$ID2," ",consensus_mean$UMI)))

  fileSplit <- as.data.table(str_split(filepath2,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")


  output.csv <- as.data.table(cbind(UMI = consensus_mean$UMI,
                                    UMIs = newUMIs$UMI,
                                    counts = newUMIs$Counts,
                                    read1 = consensus_mean$read1,
                                    quality1 = consensus_mean$quality1,
                                    read2 = consensus_mean$read2,
                                    quality2 = consensus_mean$quality2,
                                    ID1 = paste0(newUMIs$ID1," ",consensus_mean$UMI),
                                    ID2 = paste0(newUMIs$ID2," ",consensus_mean$UMI)))

  write.table(output.csv,paste0(outputsFolder,"/",part,"_summary_table.csv"), sep ="\t", row.names = F)

  remove(part, file, output, fileSplit, output.csv)
  
}

pairedR1R2 <- function(filepath1, 
                       filepath2, 
                       outputsFolder, 
                       UMIlength, 
                       UMIdistance, 
                       sequenceLength, 
                       sequenceDistance,
                       countsCutoff){
  
  dir.create(outputsFolder, showWarnings = FALSE)
  file.create(paste0(outputsFolder, "/extra_info.txt"))
  memory <- c() 
  memory[1] <-  mem_used()
  
  # read input files
  reads1 <- ShortRead::readFastq(filepath1)
  reads2 <- ShortRead::readFastq(filepath2)
  
  start.time <- Sys.time()
  
   #reads1 <- reads1[1:1000]
   #reads2 <- reads2[1:1000]

  # File 1
  seq <- data.table::as.data.table(ShortRead::sread(reads1))
  ids <- data.table::as.data.table(reads1@id)
  
  full <- base::cbind(seq, ids)
  names(full) <- c("seq", "id")
 
  # separate UMI and read
  full$UMI <- base::substring(full$seq, 1, UMIlength)
  full$read <- base::substring(full$seq, UMIlength + 1, sequenceLength)
  
  # File 2
  seq2 <- data.table::as.data.table(ShortRead::sread(reads2))
  ids2 <- data.table::as.data.table(reads2@id)
  
  full2 <- cbind(seq2, ids2)
  names(full2) <- c("seq", "id")
  
  # separate UMI and read
  full2$UMI <- base::substring(full2$seq, 1, UMIlength)
  full2$read <- base::substring(full2$seq, UMIlength + 1, sequenceLength)
  
  # keep the information of the read2 UMI, in the read1 also
  full$UMI12 <- paste0(full$UMI, full2$UMI)
  
  # keep intermediate information 
  partIDs <- data.table::as.data.table(cbind(UMI12 = full$UMI12, 
                                             ID1 = full$id, 
                                             ID2 = full2$id))
  
  intermediate.table <- partIDs[,.(count = .N), by = UMI12]
  
  partIDs <- partIDs[!duplicated(UMI12),]
  partIDs <- partIDs[order(partIDs$UMI12, decreasing = TRUE), ]
  
  intermediate.table <- intermediate.table[order(intermediate.table$UMI12, decreasing = TRUE), ]
  intermediate.table$ID1 <- partIDs$ID1
  intermediate.table$ID2 <- partIDs$ID2
  
  intermediate.table <- intermediate.table[which(intermediate.table$count >= countsCutoff), ]
  
  rm(partIDs)
  memory[2] <-  mem_used()
  # data preparation
  full <- tidyr::separate(full, id, c("id1", "id2"), " ", remove = T)
  full <- dplyr::select(full, read, id1, UMI, UMI12)
  colnames(full) <- c("read", "id", "UMI", "UMI12")
    
  quality <- as(Biostrings::quality(reads1), "matrix")
  quality <- as.data.table(quality)
  quality <- quality[,(UMIlength + 1):sequenceLength]
  quality$id <- full$id
  
  # File2
  full2 <- tidyr::separate(full2, id, c("id1", "id2"), " ", remove = T)
  full2 <- dplyr::select(full2, read, id1, UMI)
  colnames(full2) <- c("read", "id", "UMI")
  
  quality2 <- as(Biostrings::quality(reads2), "matrix")
  quality2 <- as.data.table(quality2)
  quality2 <- quality2[,(UMIlength+1):sequenceLength]
  quality2$id <- full2$id
  
  rm(ids, ids2, seq, seq2, reads1, reads2)
  
  # first consensus for each unique pair
  
  # result_mean = list()
  
  cat("step 1\n")
  
  intermediate.table = setDT(intermediate.table)
  full               = setDT(full)
  full2              = setDT(full2)
  quality            = setDT(quality)
  quality2           = setDT(quality2)
  
  result_mean = groupingPairedR1R2(intermediate.table, 
                                   full,
                                   quality,
                                   full2,
                                   quality2)
  
  memory[3] <-  mem_used()
  
  # for(i in 1:nrow(intermediate.table)){
  #   
  #   result_mean[[i]] = groupingPairedR1R2(intermediate.table$count[i], 
  #                                         full, 
  #                                         quality, 
  #                                         full2, 
  #                                         quality2, 
  #                                         intermediate.table$UMI12[i])
  #   
  # }
  
  # result_mean <- bind_rows(result_mean)
  
  # UMI correction
  # checks both UMI1 and UMI2
  
  cat("step 2\n")
  
  newUMIs <- UMIcorrectionPairedR1R2(intermediate.table,
                                     result_mean,
                                     sequenceDistance, 
                                     UMIdistance, 
                                     UMIlength,
                                     outputsFolder)
  
  rm(intermediate.table)
  memory[4] <-  mem_used()
  # final consensus
  
  cat("step 3\n")
  
  consensus_mean = groupingFinalPairedR1R2(newUMIs,
                                           full,
                                           quality,
                                           full2, 
                                           quality2,
                                           result_mean,
                                           UMIlength)
  memory[5] <-  mem_used()
  
  consensus_mean = consensus_mean[order(consensus_mean$UMI), ]
  result_mean    = result_mean[order(result_mean$UMI), ]
  
  newUMIs$key    = str_sub(newUMIs$UMI, 1, 2 * UMIlength)
  newUMIs        = newUMIs[order(newUMIs$key), ]
  
  # consensus_mean = list()
  # 
  # for(i in c(1:nrow(newUMIs))){
  #   
  #   consensus_mean[[i]] = groupingFinalPairedR1R2(newUMIs$UMI[i], 
  #                                                 full, 
  #                                                 quality, 
  #                                                 full2, 
  #                                                 quality2,
  #                                                 result_mean, 
  #                                                 UMIlength)
  #   
  # }
  # 
  # consensus_mean <- bind_rows(consensus_mean)
  
  # produce Outputs 
  
  
  # File1
  end.time <- Sys.time()
  
  total.time.secs <- difftime(end.time,start.time,units = "secs")
  total.time.mins <- difftime(end.time,start.time,units = "mins")
  
  line = paste0("Total time in secs: ", total.time.secs)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Total time in mins: ", total.time.mins)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  
  line = paste0("Memory used at start: ", memory[1])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after data cleaning: ", memory[2])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after first consensus: ", memory[3])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after UMI merging: ", memory[4])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after second consensus: ", memory[5])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(newUMIs$ID1, " ", substring(consensus_mean$UMI, 1, UMIlength))))
  
  
  fileSplit <- as.data.table(str_split(filepath1, "\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)], "\\."))
  output <- paste0(outputsFolder, "/", fileSplit[1], "_corrected.fastq.gz")
  part <- fileSplit[1]
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  # File2
  file <- ShortReadQ(DNAStringSet(consensus_mean$read2), 
                     FastqQuality(consensus_mean$quality2),
                     BStringSet(paste0(newUMIs$ID2, " ", substring(consensus_mean$UMI, UMIlength + 1, 2 * UMIlength))))
  
  fileSplit <- as.data.table(str_split(filepath2, "\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder, "/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  output.csv <- as.data.table(cbind(UMI = consensus_mean$UMI,
                                    UMIs = newUMIs$UMI,
                                    counts = newUMIs$Counts,
                                    read1 = consensus_mean$read1,
                                    quality1 = consensus_mean$quality1,
                                    read2 = consensus_mean$read2,
                                    quality2 = consensus_mean$quality2,
                                    ID1 = paste0(newUMIs$ID1, " ", substring(consensus_mean$UMI, 1, UMIlength)),
                                    ID2 = paste0(newUMIs$ID2, " ", substring(consensus_mean$UMI, UMIlength + 1, 2 * UMIlength))))
  
   write.table(output.csv,paste0(outputsFolder, "/", part, "_summary_table.csv"), sep = "\t", row.names = F)
  
  remove(part, file, output, fileSplit, output.csv)
  
}

single <- function(filepath1, 
                   outputsFolder, 
                   UMIlength, 
                   UMIdistance, 
                   sequenceLength, 
                   sequenceDistance,
                   countsCutoff){
  
  dir.create(outputsFolder)
  file.create(paste0(outputsFolder, "/extra_info.txt")) 
  
  memory <- c() 
  memory[1] <-  mem_used()
   
  #read input file
  reads1 <- readFastq(filepath1)
  
  start.time <- Sys.time()
  start.all <- proc.time()
  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(id(reads1))
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
  
  #separate UMI and read
  full$UMI  <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength + 1, sequenceLength)
  
  # full <- select(full, read, id, UMI)
  # colnames(full) <- c("read", "id", "UMI")
  
  full = full[, c("read", "id", "UMI"), with = FALSE]
  
  
  #keep intermediate information 
  partIDs <- data.table(
      "UMI" = full$UMI,
      "ID1" = full$id
  ) # as.data.table(cbind(UMI = full$UMI, ID1 = full$id))
  
  partIDs <- partIDs[which(!duplicated(UMI)),]
  partIDs <- partIDs[order(partIDs$UMI, decreasing = TRUE), ]
  
  intermediate.table     <- full[, by = UMI, .(count = .N)]
  intermediate.table     <- intermediate.table[order(intermediate.table$UMI, decreasing = TRUE), ]
  intermediate.table$ID1 <- partIDs$ID1
  intermediate.table     <- intermediate.table[which(intermediate.table$count >= countsCutoff),]
  
  rm(partIDs)
  memory[2] <-  mem_used()
  intermediate.table = intermediate.table[order(intermediate.table$count, decreasing = TRUE), ]
  
  #data preparation
  quality <- as(quality(reads1), "matrix")
  quality <- setDT(as.data.frame(quality))
  quality <- quality[,(UMIlength + 1):sequenceLength]
  quality$id <- full$id

  rm(ids, seq, reads1)
  
  #first consensus
  result_mean = groupingSingle(
      intermediate.table, # $UMI[i], 
      # intermediate.table$count[i], 
      full, 
      quality,  
      UMIlength
  )
  
  memory[3] <-  mem_used()
  #UMI correction
  newUMIs <- UMIcorrectionSingle(
      intermediate.table,
      result_mean,
      sequenceDistance, 
      UMIdistance,
      outputsFolder
  )
  
  newUMIs <- newUMIs[which(!is.na(UMI))]
  
  rm(intermediate.table)
  memory[4] <-  mem_used()
  #final consensus
  consensus_mean = groupingFinalSingle(
      newUMIs, # i, 
      full, 
      quality, 
      result_mean, 
      UMIlength
  )
  
  memory[5] <-  mem_used()
  
  consensus_mean = consensus_mean[order(consensus_mean$UMI), ]
  result_mean    = result_mean[order(result_mean$UMI), ]
  
  newUMIs$key = str_sub(newUMIs$UMI, 1, UMIlength)
  newUMIs        = newUMIs[order(newUMIs$key), ]
  #first consensus
  #result_mean = list()
  
  #for(i in c(1:nrow(intermediate.table))){
    
   # result_mean[[i]] = groupingSingle(intermediate.table$UMI[i], intermediate.table$count[i], full, quality, UMIlength)
    
  #}
  
  #result_mean <- bind_rows(result_mean)
  
  #UMI correction
  #newUMIs <- UMIcorrectionSingle(intermediate.table,result_mean,UMIdistance, sequenceDistance)
  #rm(intermediate.table)
  
  #final consensus
  #consensus_mean = list()
  
  #for(i in newUMIs$UMI){
    
   # consensus_mean[[i]] = groupingFinalSingle(i, full, quality,result_mean, UMIlength)
    
  #}
  
  #consensus_mean <- bind_rows(consensus_mean)
  
  #produce Outputs 
  
  end.time <- Sys.time()
  end.all  <- proc.time()
  
  total.time.secs <- difftime(end.time,start.time,units = "secs")
  total.time.mins <- difftime(end.time,start.time,units = "mins")
  time.fifference <- end.all - start.all
  
  line = paste0("Total time in secs: ", total.time.secs)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Total time in mins: ", total.time.mins)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("User time: ", time.fifference[1])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("System time: ", time.fifference[2])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Elapsed time: ", time.fifference[3])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used at start: ", memory[1])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
 
  line = paste0("Memory used after data cleaning: ", memory[2])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after first consensus: ", memory[3])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after UMI merging: ", memory[4])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  line = paste0("Memory used after second consensus: ", memory[5])
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  #File1
  file <- ShortReadQ(
      DNAStringSet(consensus_mean$read1), 
      FastqQuality(consensus_mean$quality1),
      BStringSet(paste0(newUMIs$ID1, " ", consensus_mean$UMI))
  )
  
  
  fileSplit <- setDT(str_split(filepath1, "\\/"))
  fileSplit <- setDT(str_split(fileSplit[nrow(fileSplit)], "\\."))
  
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  
  part <- fileSplit[1]
  file.create(output)
  writeFastq(file, output, mode = "a")
 
  output.csv <- data.table(
      "UMI"      = consensus_mean$UMI,
      "UMIs"     = newUMIs$UMI,
      "counts"   = newUMIs$Counts,
      "read1"    = consensus_mean$read1,
      "quality1" = consensus_mean$quality1,
      "ID1"      = paste0(newUMIs$ID1," ",consensus_mean$UMI)
  )
  
  fwrite(
      output.csv,
      paste0(outputsFolder, "/", part,"_summary_table.csv"), 
      sep = ",", row.names = F, quote = TRUE
  )
 
  remove(part, file,output, fileSplit, output.csv)
  
}