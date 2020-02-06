###
#   File name : MakeRawCounts.R
#   Author    : Hyunjin Kim
#   Date      : Feb 5, 2020
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Do alignment, make Bam files, and make raw count matricies originally from fastq.gz files
#
#   * THIS CODE SHOULD BE RUN ON LINUX
#
#   * THIS IS HUMAN DATA, HENCE HG38 WILL BE USED AS REFERENCE
#
#   Instruction
#               1. Source("MakeRawCounts.R")
#               2. Run the function "makeRCnt" - specify the input directory (fastq.gz) and output directory
#               3. The Bam files and raw counts will be generated under the output directory
#
#   Example
#               > source("The_directory_of_MakeRawCounts.R/MakeRawCounts.R")
#               > makeRCnt(fastqgzPath="/mnt/c/Research/CUMC/Giovanni/Alexander/For_GEO/firadazer/raw_data_files/",
#                          referencePath="/mnt/e/Reference/hg38.fa",
#                          referenceIdxPath="/mnt/e/Reference/hg38.index",
#                          readType=c("single-end", "paired-end"),
#                          outputDir="/mnt/c/Research/CUMC/Giovanni/Alexander/For_GEO/")
###

makeRCnt <- function(fastqgzPath="//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/af_lab/hk2990/Ferrari_data/raw_data_files/",
                     referencePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Hyunjin/Reference/hg38.fa",
                     referenceIdxPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Hyunjin/Reference/hg38.index",
                     readType=c("single-end", "paired-end"),
                     outputDir="//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/af_lab/hk2990/Ferrari_data/") {
  
  ### load library
  if(!require(Rsubread, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Rsubread", version = "3.8")
    require(Rsubread, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db", version = "3.8")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  ### build the index of the reference genome
  ### already built, therefore, no need to run this line
  # buildindex(basename=referenceIdxPath, reference=referencePath)
  
  ### get a list of fastq files
  fastqFiles <- list.files(fastqgzPath, recursive = TRUE)
  fastqFiles <- fastqFiles[which(endsWith(fastqFiles, ".fq.gz"))]
  
  ### create new directory for Bam files
  dir.create(paste0(fastqgzPath, "../bam_files/"), showWarnings = FALSE)
  
  ### iteratively perform alignment
  if(readType[1] == "paired-end") {
    sample_names <- unique(c(sapply(fastqFiles[grep("R1", fastqFiles)],
                                    function(x) {
                                      strsplit(x, "R1", TRUE)[[1]][1]                                      
                                    }),
                             sapply(fastqFiles[grep("R2", fastqFiles)],
                                    function(x) {
                                      strsplit(x, "R2", TRUE)[[1]][1]
                                    })))
    sample_names <- sapply(sample_names, function(x) substr(x, 1, nchar(x)-1), USE.NAMES = FALSE)
    
    for(samp in sample_names) {
      align(index=referenceIdxPath,
            readfile1=paste0(fastqgzPath, samp, "_R1.fastq.gz"),
            readfile2=paste0(fastqgzPath, samp, "_R2.fastq.gz"),
            output_file=paste0(fastqgzPath, "../bam_files/", samp, ".bam"),
            nthreads = 4)
    }
  } else if(readType[1] == "single-end") {
    for(i in 1:length(fastqFiles)) {
      align(index=referenceIdxPath,
            readfile1 = paste0(fastqgzPath, fastqFiles[i]),
            readfile2 = NULL,
            output_file = paste0(fastqgzPath, "../bam_files/", strsplit(fastqFiles[i], ".fastq", fixed = TRUE)[[1]][1], ".bam"),
            nthreads = 4)
    }
  } else {
    stop("The \"readType\" parameter should be \"paired-end\" or \"single-end\".")
  }
  
  ### get a list of created bam files
  bamFiles <- list.files(paste0(fastqgzPath, "../bam_files/"))
  bamFiles <- bamFiles[which(endsWith(bamFiles, ".bam"))]
  
  ### iteratively get counts from the bam files
  rawCnt <- featureCounts(files = paste0(fastqgzPath, "../bam_files/", bamFiles),
                          annot.inbuilt = "hg38",
                          isPairedEnd = TRUE)
  rawCnt <- rawCnt$counts
  
  ### numerize the raw counts (now they are characters)
  rawCnt <- as.data.frame(rawCnt)
  rawCnt[1:ncol(rawCnt)] <- lapply(rawCnt[1:ncol(rawCnt)], function(x) as.numeric(as.character(x)))
  
  ### set column names for rawCnt and annotate gene symbols
  colnames(rawCnt) <- c("1106_AV", "1204_AV", "1354_AV", "1722_AV", "1770_AV", "1815_AV")
  map_eg_symbol <- mappedkeys(org.Hs.egSYMBOL)
  list_eg2symbol <- as.list(org.Hs.egSYMBOL[map_eg_symbol])
  rawCnt <- cbind(Entrez_ID=rownames(rawCnt),
                  Gene_Symbol=as.character(list_eg2symbol[rownames(rawCnt)]),
                  rawCnt)
  
  ### sample information
  sample_info <- c(rep("Stenosis", 3), rep("Sclerosis", 3))
  names(sample_info) <- colnames(rawCnt)[-c(1,2)]
  
  ### write out the result
  write.table(rawCnt, file = paste0(outputDir, "raw_counts.txt"),
              sep = "\t", row.names = FALSE)
  save(list = c("rawCnt", "sample_info"), file = paste0(outputDir, "raw_counts.rda"))
  
}
