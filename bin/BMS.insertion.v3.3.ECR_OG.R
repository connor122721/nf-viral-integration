
# List all installed packages
installed_packages <- installed.packages()

# Extract the names of the installed packages
package_names <- installed_packages[, "Package"]

# Load all installed packages
lapply(package_names, library, character.only = TRUE)


  ###############################
  ## BMS.Insertion.v2.2.6.ECR  ##
  ## MODIFIED: 01/03/2024      ##
  ## MODIFIED BY: Eric Rouchka ##
  ## ORIGINAL CODE BY:         ##
  ##    Steven Nadakal         ##
  ##    Justin Kos             ##
  ###############################
  
  #################################################################
  ## Global Variables -- May need to change based on the SMRTSeq ##
  ##                     Run type                                ##
  #################################################################
  #REF_SIZE_DIR <- "/home/kyinbre/data/Genomes/chrSizes/"      ## Location of reference chromosome sizes
  options(java.parameters = "-Xmx8000m")
  REF_SIZE_DIR <-  "/home/stnada01/scratch/forConnor/sizeReferences/" #E:\\SMRTCap\\"
  coordinateDir <- "/home/stnada01/scratch/forConnor/hg38/"
  insertType <- "HIV"       ## Others: BMS, SIV, ...
  refGenome  <- "Hs"    ## Others: RheMac, ...
  excludeGeneList <- c("NGFR", "CD24P2", "IRES3") ## Genes to exclude from analysis #Edited 14 April 2025 for V. Simon. Remove list as needed.
  tolerance <- 10           ## Tolerance to make identical calls in the context of indels
  exciseFlag = FALSE       ## Flag for excising same reads -- set to TRUE for BMS
  hostTolerance <- 100000    ## Max number of bases difference in host genome to filter 
  minReadLength <- 3000      ## Min number of bases for total read length to assume true (instrument dependent)
  maxPCRSlippage <- 10       ## Max number of bases for PCR slippage
  templateToRemove <- tools::file_path_sans_ext(list.files(pattern = "\\.csv$"))
  ADD_CHR_FLAG = TRUE      ## True to add chr to chromosome names; False to keep as normal
  USE_ENSEMBL <- FALSE       ## Identifies whether to download genomes from Ensembl or locally
  #RDS_DIR <- "/home/kyinbre/data/"
  #RDS_DIR <- "E:\\SMRTCap\\Mmul_10_Annotations\\"
  #biocPackages <- c("GenomicRanges", "rtracklayer", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
   #                 "ensembldb", "EnsDb.Hsapiens.v86", "annotatr", "AnnotationHub",
   #                 "S4Vectors", "org.Hs.eg.db")
  #biocPackageVersions <- c("1.52.0", "1.60.1", "3.17.0", "2.24.1", "2.99.0", "1.26.0", "3.8.0", "0.38.2", "3.17.0")
  
  #RPackages <- c("ggplot2", "dplyr", "stringr", "readxl", "ggforce", "purrr",
  #               "stringi", "chromoMap", "randomcoloR", "data.table", 
  #               "patchwork", "GenomicRanges", "ensembldb", "plyranges",
  #               "S4Vectors", "AnnotationHub", "annotatr", "tidyr", "xlsx", "ids", "tictoc", "circlize")
  #RPackageVersions <- c("3.4.3", "1.1.3", "1.5.0", "1.4.3", "0.4.1", "1.0.2",
  #                      "1.7.12", "4.1.1", "1.1.0.1", "1.14.8",
  #                      "1.1.3", "1.52.0", "2.24.1", "1.20.0",
  #                      "0.38.2", "3.8.0", "1.26.0", "1.3.0", "0.6.5", "1.0.1", "1.2", "0.4.15")
  
  finalNames <- c("seqnames", "start", "end", "width", "strand", 
                  "chromosome", "fivePrimeInsertion","threePrimeInsertion", "fivePrimeShear", "threePrimeShear", 
                  "Flanking", "slippage", "passQualityCheck", "first.failedStep", "first.READ", 
                  "Clones", "Clone.ID", "PCR.Dups", "PCR.ID", "INSERT_LEN",
                  "Human.Gene.Directionality", "subject", "is.5PI.na", "is.3PI.na", "annotation", "second.seqnames",
                  "second.start", "second.end", "second.width", "second.strand", "second.id",
                  "second.type", "second.region", "second.symbol", "second.counts")
  
  cat("Genes being excluded for integration analysis are:\n", excludeGeneList)
  
  #-----------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------
  getReferenceSizeFile <- function(REF_SIZE_DIR, refGenome) {
  
     #################################################
     ## Get the name of the file for the chromosome ##
     ## sizes based on the reference genome         ##
     #################################################
      tic("getReferenceSizeFile")
  
     refSizeFile = "hg38.size.references.txt"  ## Hs by default
     if(refGenome == "RheMac") {
        refSizeFile = "Mmul_10.size.references.txt"
     }
     refSizeFile = paste(REF_SIZE_DIR, refSizeFile, sep="/")
     elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
     elapsed_time
     return(refSizeFile)
  }
  #-----------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------
  getChromosomeOrder <- function(ref.Genome) {
    #############################################
    ## gets the chromosomal order based on the ##
    ## reference genome                        ##
    #############################################
    
    tic("getChromosomeOrder")
    chrOrder <- c((1:22), "X", "Y", "M")   ## Human by default
    if(ref.Genome == "RheMac") {
      chrOrder <- c((1:20), "X", "Y", "M") ## Rhesus macaque
    }
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(chrOrder)
  }
  #-----------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------
  getReferenceInsertSize <- function(iType) { 
    
    ###########################################
    ## Returns the size of the reference for ##
    ## the inserted viral sequence           ##
    ## Values                                ##
    ##    HIV: 9710                          ##
    ##    BMS: 9382                          ##
    ##    SIV: 10535                         ##
    ###########################################
    
    tic("getReferenceInsertSize") ######PUT IN NAMES
    refSize <- 9710                     ## HIV by default
    if(insertType == "BMS") {
      refSize <- 9382                  ## Bristol Myers Squibb insert
    }
    if(insertType == "SIV") {
      refSize <- 10535                 ## SIV virus isolate 239 ASM319262v1
    }
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(refSize)
  }
  #-----------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------
  #installAndLoadPackage <- function(pkg, version) {
    # Check if package is installed
   # if (!requireNamespace(pkg, quietly = TRUE)) {
      # If not installed, install the package
    #  install.packages(pkg, dependencies = TRUE)
    #} else {
      # If installed, check if the version matches
     # installed_version <- packageVersion(pkg)
      #if (installed_version != version) {
        # If version doesn't match, reinstall the package
       # install.packages(pkg, dependencies = TRUE)
      #}
    #}
    # Load the package
    #library(pkg, character.only = TRUE)
  #}
  
  
  #-----------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------
  #loadLibraries <- function(pkgList, pkgVersions) { 
    ####################################################
    ## Loads the required packages ##
    ####################################################
   # cat("Entering loadLibraries...\n")
    
    #for (i in 1:length(pkgList)) { 
     # currPkg <- pkgList[i]
      #version <- pkgVersions[i]
      
      #installAndLoadPackage(currPkg, version)
      #cat("Loaded package", currPkg, "version", version, "\n")
    #}
    
    #old version!
    #cat("Entering loadLibraries...\n")
    #for(i in 1:length(pkgList1)) { 
    #  currPkg <- pkgList1[i]
    #  library(currPkg, character.only = TRUE)
    #}
    #for(i in 1:length(pkgList2)) { 
    #  currPkg <- pkgList2[i]
    #  library(currPkg, character.only = TRUE)
    #}
 # }
  #-----------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------
 # InstallBioconductorPackages <- function (pkgList) {
      
      #########################################
      ## Check to see if bioconductor exists ##
      ## and install if it does not          ##
      #########################################
  #    tic("InstallBioconductorPackages")
   #   cat("Entering InstallBioconductorPackages...\n")
    #  if(!require("BiocManager", quietly = TRUE)) { 
     #    install.packages("BiocManager")
     # }
     # library("BiocManager")
     # for(i in 1:length(pkgList)) { 
      #   currPkg <- pkgList[i]
  #       if(!require(currPkg, character.only=TRUE, quietly = TRUE)) { 
  #          BiocManager::install(currPkg, update = FALSE, ask = TRUE)
  #       }
  #    }
  #    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
  #    elapsed_time
  #}
  #-----------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------
 # InstallRPackages <- function(pkgList) { 
      
     ###########################################################
     ## checks to see if each CRAN R package is installed and ##
     ## installs it if it is not found                        ##
     ###########################################################
    #tic("InstallRPackages")
  #   cat("Entering InstallRPackages...\n")
   #  for(i in 1:length(pkgList)) { 
    #    currPkg <- pkgList[i]
     #   if(!require(currPkg, character.only=TRUE, quietly = TRUE)) { 
      #     install.packages(currPkg, character.only=TRUE, update = FALSE, ask = TRUE)
      #  }
     #}
     #elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
     #elapsed_time
  #}
  #-----------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
 getReadsFromCSVFiles <- function() {
    ##########################################################################################
    ## Import excel files (converted to csv) into a list of dataframes for speed and lapply ##
    ## purposes -- excel files cannot be read if they are currently opened                  ##
    ##########################################################################################
    tic("getReadsFromCSVFiles")
    cat("importing files getReadsFromCSVFiles...\n")
    file.list            <- list.files(pattern='*.csv')
    names                <- list.files(pattern='*.csv')
    raw.read.list        <- lapply(file.list, fread)
    names(raw.read.list) <- names
    ##################################################################
    ## Add a column to each dataframe that contains the name of the ##
    ## source excel sheet/data file -- "Run ID"                     ##
    ##################################################################
    cat("creating dataframe\n")
    reads <- Map(cbind, raw.read.list, name=names(raw.read.list))
    reads <- bind_rows(raw.read.list, .id = "Run_ID")
    setDF(reads)
    setattr(reads, "class", c("tbl_df", "tbl", "data.frame"))
    reads[reads==""]<-NA
    input_stringL <- reads$LEFT_FLANK
    input_stringR <- reads$RIGHT_FLANK
    # Split the input string by ':' and '-'
    #stn edits
    if (all(is.na(input_stringL)) && all(is.na(input_stringR))) {
      stop("There are no flank coordinates. \n")
    }
    if (!all(is.na(input_stringL)) && !all(is.na(input_stringR))) {
      split_parts <- unlist(strsplit(input_stringR, "[:-]"))
    }
    if (all(is.na(input_stringL)) | all(is.na(input_stringR))) {
      if (all(is.na(input_stringL))) {
        cat("There are no left flank coordinates.\n")
        split_parts <- unlist(strsplit(input_stringR, "[:-]"))
      } else if (all(is.na(input_stringR))) {
        cat("There are no right flank coordinates.\n")
      split_parts <- unlist(strsplit(input_stringL, "[:-]"))}}
    #else (split_parts <- unlist(strsplit(input_stringL, "[:-]")))
    # Join the split parts with a tab separator
    output_string <- paste(split_parts, collapse = "\t")
    # Print the formatted output
    cat(output_string, "\n")
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(reads)
  } 


  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  Partition <- function(dataframe) {
  
    #MAPPING QUESTION -- WHAT TO DO WHEN FLANKS OVERLAP? ADD IN FLAG
    
    #######################################################
    ## Filters into Left Flanked, Right Flanked and Both ##
    ## Flanked subsets and adds chromosomes              ##
    #######################################################
    cat("*************************************************\n")
    cat("____________________________Entering Partition...\n")
    cat("*************************************************\n")
    cat("Starting read count:", nrow(dataframe), "\n")
    L <- dplyr::filter(dataframe,!is.na(LEFT_FLANK) & is.na(RIGHT_FLANK))  ## Only has left flank
    R <- dplyr::filter(dataframe,is.na(LEFT_FLANK) & !is.na(RIGHT_FLANK))  ## Only has right flank
    B <- dplyr::filter(dataframe,!is.na(LEFT_FLANK) & !is.na(RIGHT_FLANK)) ## Has both flanks
    N <- dplyr::filter(dataframe,is.na(LEFT_FLANK) & is.na(RIGHT_FLANK)) #has no flanks?
    cat(nrow(L), "5' Flanked Reads.\n") 
    cat(nrow(R), "3' Flanked Reads.\n") 
    cat(nrow(B), "Dually Flanked Reads.\n") 
    cat(nrow(N), "No Flanked Reads.\n") 
    drops <- c("Left Prior", 'Left Omit', 'Right Prior', 'Right Omit', 'lchromosome', 'rchromosome')
   
    ################################################################ 
    ## Adding in the chromosomes for each subset --               ##
    ## the Both subset incorporates a sanity check to make sure   ##
    ## the chromosome assigned to each flank agrees.              ## 
    ## str_split_fixed splits coordinates by text features, with  ##
    ## the omit column splitting off unnecessary text             ##
    ################################################################ 
    
    #######################################
    ## ADDRESSING THE READS THAT FAIL QC ##
    #######################################
    
    if(nrow(N) > 0) {  ## Check to make sure there are some Left reads ##
      N$passQualityCheck <- FALSE
      N$failedStep <- "cleanReads"
      N$Human.Gene.Directionality <-  "Nonsense"
    } else {
      ## No left reads -- create blank columns ##
      cat("There are no non-flanked reads worth addressing.")
    }
  
    #####################################
    ## ADD LEFT CHROMOSOME INFORMATION ##
    #####################################
    if(nrow(L) > 0) {  ## Check to make sure there are some Left reads ##
       L$passQualityCheck <- TRUE
       L$chromosome <- word(L$LEFT_FLANK,1,sep = ":")
       L$chromosome <- sapply(L$chromosome , function(x) gsub("chr", "",  x))
       L[c('Left Prior', 'fivePrimeInsertion')] <- str_split_fixed(L$LEFT_FLANK, pattern = '-', n = 2)
       L[c('Left Omit', 'fivePrimeShear')] <- str_split_fixed(L$`Left Prior`, pattern = ':', n = 2)
    } else {
       ## No left reads -- create blank columns ##
       L$chromosome       <- character(0)
       L$`Left Prior`     <- character(0)
       L$fivePrimeInsertion       <- character(0)
       L$`Left Omit`      <- character(0)
       L$fivePrimeShear <- character(0)
    }
    L <- L[ , !(names(L) %in% drops)] 
    L$fivePrimeShear <- as.numeric(L$fivePrimeShear)
    L$fivePrimeInsertion <- as.numeric(L$fivePrimeInsertion)
    if(nrow(L) > 0) {
       L$Human.Gene.Directionality <-  ifelse(L$fivePrimeShear <= L$fivePrimeInsertion, yes = "Sense", no = "Antisense")
    } else {
       L$Human.Gene.Directionality <- character(0)
    }
    cat("Done with Left\n")
    ###################################### 
    ## ADD RIGHT CHROMOSOME INFORMATION ## 
    ###################################### 
    if(nrow(R) > 0) {
       R$passQualityCheck <- TRUE
       R$chromosome <- word(R$RIGHT_FLANK,1,sep = ":")
       R$chromosome  <- sapply(R$chromosome , function(x) gsub("chr", "",  x))
       #R[c('Beginning', 'End')] <- str_split_fixed(R$RIGHT_FLANK, pattern = '-', n = 2)
       R[c('Right Prior', 'threePrimeShear')] <- str_split_fixed(R$RIGHT_FLANK, pattern = '-', n = 2)
       R[c('Right Omit', 'threePrimeInsertion')] <- str_split_fixed(R$`Right Prior`, pattern = ':', n = 2)
    } else {
      ## No right reads -- create blank columns ##
       R$chromosome       <- character(0)
       R$`Right Prior`     <- character(0)
       R$threePrimeShear       <- character(0)
       R$`Right Omit`      <- character(0)
       R$threePrimeInsertion <- character(0)   
    }
    R <- R[ , !(names(R) %in% drops)] 
    R$threePrimeInsertion <- as.numeric(R$threePrimeInsertion)
    R$threePrimeShear <- as.numeric(R$threePrimeShear)
    if(nrow(R) > 0) {
    R$Human.Gene.Directionality <-  ifelse(R$threePrimeInsertion <= R$threePrimeShear, yes = "Sense", no = "Antisense")
    } else {
       R$Human.Gene.Directionality <- character(0)
    }
    cat ("Done with Right\n")
    
    #################################### 
    ## ADD BOTH CHROMOSME INFORMATION ##
    #################################### 
    if(nrow(B) > 0) { 
       B$passQualityCheck <- TRUE
       B$lchromosome <- word(B$LEFT_FLANK,1,sep = ":")
       B$rchromosome <- word(B$RIGHT_FLANK,1,sep = ":")
       B$chromosome <- ifelse(B$lchromosome == B$rchromosome, B$lchromosome,"ERROR")
       B$chromosome  <- sapply(B$chromosome , function(x) gsub("chr", "",  x))
       B[c('Left Prior', 'fivePrimeInsertion')] <- str_split_fixed(B$LEFT_FLANK, pattern = '-', n = 2)
       B[c('Left Omit', 'fivePrimeShear')] <- str_split_fixed(B$`Left Prior`, pattern = ':', n = 2)
       B[c('Right Prior', 'threePrimeShear')] <- str_split_fixed(B$RIGHT_FLANK, pattern = '-', n = 2)
       B[c('Right Omit', 'threePrimeInsertion')] <- str_split_fixed(B$`Right Prior`, pattern = ':', n = 2)
    }
    else {
      ## No both flanking reads -- create blank columns ##
      B$chromosome        <- character(0)
      B$lchromosome       <- character(0)
      B$`Left Prior`      <- character(0)
      B$fivePrimeInsertion        <- character(0)
      B$`Left Omit`       <- character(0)
      B$fivePrimeShear  <- character(0)  
      B$rchromosome       <- character(0)
      B$`Right Prior`     <- character(0)
      B$threePrimeShear       <- character(0)
      B$`Right Omit`      <- character(0)
      B$threePrimeInsertion <- character(0)   
    }
    B <- B[ , !(names(B) %in% drops)]
    B$fivePrimeShear <- as.numeric(B$fivePrimeShear)
    B$fivePrimeInsertion <- as.numeric(B$fivePrimeInsertion)
    B$threePrimeInsertion <- as.numeric(B$threePrimeInsertion)
    B$threePrimeShear <- as.numeric(B$threePrimeShear)
    
    
    cat("Identifying orientation...\t")
    if (nrow(B) > 0) { 
    row_avg_fivePrime <- rowMeans(B[, c("fivePrimeShear", "fivePrimeInsertion")], na.rm = TRUE)
    row_avg_threePrime <- rowMeans(B[, c("threePrimeShear", "threePrimeInsertion")], na.rm = TRUE)
    
    B$Human.Gene.Directionality <- ifelse(
      row_avg_fivePrime < row_avg_threePrime, 
      yes = "Sense", no = "Antisense")
    } else {
       cat("was this the error?")
      #B$Human.Gene.Directionality <- character(0)
    }
    cat("Organizing flank coordinates according to orientation...\n")
    
    if (nrow(B) > 0) {
      row_avg_fivePrime <- rowMeans(B[, c("fivePrimeShear", "fivePrimeInsertion")], na.rm = TRUE)
      row_avg_threePrime <- rowMeans(B[, c("threePrimeShear", "threePrimeInsertion")], na.rm = TRUE)
      
      for (i in seq_along(B$Human.Gene.Directionality)) {
        #to say this part of the code didn't break my brain and cost me years of my life expectancy would be a bold lie. 
        if (row_avg_fivePrime[i] <= row_avg_threePrime[i]) {
          # Sense orientation
          B$fivePrimeShear[i] <- pmin(B$fivePrimeShear[i], B$fivePrimeInsertion[i])
          B$fivePrimeInsertion[i] <- pmax(B$fivePrimeShear[i], B$fivePrimeInsertion[i])
          B$threePrimeShear[i] <- pmax(B$threePrimeShear[i], B$threePrimeInsertion[i])
          B$threePrimeInsertion[i] <- pmin(B$threePrimeShear[i], B$threePrimeInsertion[i])
        } else {
          # Antisense orientation
          B$fivePrimeShear[i] <- pmax(B$fivePrimeShear[i], B$fivePrimeInsertion[i])
          B$fivePrimeInsertion[i] <- pmin(B$fivePrimeShear[i], B$fivePrimeInsertion[i])
          B$threePrimeShear[i] <- pmin(B$threePrimeShear[i], B$threePrimeInsertion[i])
          B$threePrimeInsertion[i] <- pmax(B$threePrimeShear[i], B$threePrimeInsertion[i])
        }
      }
    }
    
    cat("Done with both\n")
    sumAll <- nrow(L) + nrow(R) + nrow(B) + nrow(N)
    cat("The sum total read count while partitioning is:", sumAll, "\n")
    if(sumAll == 0) { 
       ## There are no insert reads with flanking genomic sequences ##
       ## Therefore, there is nothing else to do with annotations   ##
       cat("There are no reads with flanking sequences.  Processing complete.\n")
       quit()
    }
    list("5'-Flanked" = L, "3'-Flanked" = R, "Dually-Flanked" = B, "Non-Flanked" = N)
  }  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  Concat <- function(list.of.dataframes) {
    #############################################################################
    ## Undoes the list after adding columns to the respective rows such that   ##
    ## it efficiently labels the type of shifting. Also adds an index variable ##
    ## for easy proof reading.                                                 ##
    #############################################################################
    cat("*************************************************\n")
    cat("_______________________________Entering Concat...\n")
    cat("*************************************************\n")
    if(nrow(list.of.dataframes$`5'-Flanked`) == 0) { ## ADD Blank row when no 5'-Flanked reads present
       list.of.dataframes$`5'-Flanked` <- list.of.dataframes$`5'-Flanked` %>% add_row(Run_ID="BLANK_ROW")
    }
    if(nrow(list.of.dataframes$`3'-Flanked`) == 0) {  ## Add Blank row when no 3'-Flanked reads present
      list.of.dataframes$`3'-Flanked` <- list.of.dataframes$`3'-Flanked` %>% add_row(Run_ID="BLANK_ROW")
    }
    if(nrow(list.of.dataframes$`Dually-Flanked`) == 0) {  ## Add Blank row when no Dually-Flanked reads present
      list.of.dataframes$`Dually-Flanked` <- list.of.dataframes$`Dually-Flanked` %>% add_row(Run_ID="BLANK_ROW")
    }
    if(nrow(list.of.dataframes$`Non-Flanked`) == 0) {  ## Add Blank row when no Dually-Flanked reads present
      list.of.dataframes$`Non-Flanked` <- list.of.dataframes$`Non-Flanked` %>% add_row(Run_ID="BLANK_ROW")
    }
    Map(cbind, list.of.dataframes, name=names(list.of.dataframes))
    cat("   reads mapped\n")
    reads <- bind_rows(list.of.dataframes, .id = "Flanking")
    cat("   reads bound\n")
    readcount <- nrow(reads)
    cat("Total number of reads:", readcount, "\n")
    
    tibble::rowid_to_column(reads, var = "INDEX")
    
  } 
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  Excise <- function(dataframe, t) {
    #######################################################################
    ## IF THE EXPERIMENT IS ONE THAT HAS INTRINSIC NOISE --              ##
    ## SUCH AS THE REMOVAL OF VECTORS -- UNCOMMENT THE EXCISE FUNCITON   ##
    ## These are the conditions to omit!                                 ##
    ## EXAMPLE:                                                          ##
    ## CHANGE CONDITIONS AS NEEDED.                                      ##
    ## AS OF NOW, THESE ARE THE BMS NOISE REMOVAL                        ##
    #######################################################################
  
    dataframe <- dataframe[!(dataframe$chromosome == 7 & ((dataframe$fivePrimeShear >= (55156531-t) & dataframe$fivePrimeInsertion <= (55173067+t))| 
                (dataframe$threePrimeInsertion >= (55156531-t) & dataframe$threePrimeShear <= (55173067+t)))) |
                  (dataframe$chromosome == 6 & ((dataframe$fivePrimeShear >= (73521000-t) & dataframe$fivePrimeInsertion <= (73521229+t))|
                (dataframe$threePrimeInsertion >= (73521000-t) & dataframe$threePrimeShear <= (73521229+t)))) |
                  (dataframe$chromosome == 1 & ((dataframe$fivePrimeShear >= (7920838-t) & dataframe$fivePrimeInsertion <= (7933201+t))|
                (dataframe$threePrimeInsertion >= (7920838-t) & dataframe$threePrimeShear <= (7933201+t))))   |
                  (dataframe$chromosome == 24 & ((dataframe$fivePrimeShear >= (1282704-t) & dataframe$fivePrimeInsertion <= (1282769+t))|
                (dataframe$threePrimeInsertion >= (1282704-t) & dataframe$threePrimeShear <= (1282769+t))))   |
                  (dataframe$chromosome == 23 & ((dataframe$fivePrimeShear >= (1282704-t) & dataframe$fivePrimeInsertion <= (1282769+t))|
                (dataframe$threePrimeInsertion >= (1282704-t) & dataframe$threePrimeShear <= (1282769+t))))   |
                  (dataframe$chromosome == 7 & ((dataframe$fivePrimeShear >= (22512030-t) & dataframe$fivePrimeInsertion <= (225122062+t))|
                (dataframe$threePrimeInsertion >= (22512030-t) & dataframe$threePrimeShear <= (225122062+t))))|
                  (dataframe$chromosome == 5 & ((dataframe$fivePrimeShear >= (14653357-t) & dataframe$fivePrimeInsertion <= (14653389+t))|
                (dataframe$threePrimeInsertion >= (14653357-t) & dataframe$threePrimeShear <= (14653389+t)))) |
                  (dataframe$chromosome == 9 & ((dataframe$fivePrimeShear >= (133019423-t) & dataframe$fivePrimeInsertion <= (133019455+t))|
                (dataframe$threePrimeInsertion >= (133019423-t) & dataframe$threePrimeShear <= (133019455+t)))), ]
    dataframe <- dataframe[!is.na(dataframe$Flanking), ]
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  
  statusCheck <- function(dataframe) {
    ###########################################################
    ##Quick function to calculate losses to cleaning process.##
    ###########################################################
    
    cat("\n\n\nAssessing losses...\n")
    
    goodCount <- length(dataframe$READ[dataframe$passQualityCheck == TRUE])
    badCount <- length(dataframe$READ[dataframe$passQualityCheck == FALSE])
    totalCount <- goodCount + badCount
    percentFailed <- round( x = (badCount/totalCount), digits = 2)*100
    
    
    cat("At this point in time", percentFailed, "percent, or", badCount, "reads,", "of the", totalCount, "reads are not usable for analysis. \n\n\n")
    
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  qualCreate <- function(dataframe) {
    tic("qualCreate")
    cat("*************************************************\n")
    cat("__________________________Entering qualCreate...\n")
    cat("*************************************************\n")
    
    dataframe$fivePrimeGap <- NA 
    dataframe$threePrimeGap <- NA
    sense     <- dataframe %>% dplyr::filter(Human.Gene.Directionality == "Sense")
    antisense <- dataframe %>% dplyr::filter(Human.Gene.Directionality == "Antisense")
    nonsense <- dataframe %>% dplyr::filter(Human.Gene.Directionality == "Nonsense")
    
    cat(nrow(sense), "viral integrations in the sense orientation. \n") 
    cat(nrow(antisense), "viral integrations in the antisense orientation. \n") 
    cat(nrow(nonsense), "viral integrations with an indeterminate constitution at this point in the pipeline.\n") 
    
    #Gap Coordinates are the coordinates corresponding to the shear sites -- a new columm was created solely for clarity and to minimize having to ifelse statements given orientation and flanking.
    sense <- sense %>%
      mutate(fivePrimeGap = ifelse(Flanking == "Dually-Flanked", min(fivePrimeInsertion, fivePrimeShear), 
                               ifelse(Flanking == "3'-Flanked", NA, 
                                      ifelse(Flanking == "5'-Flanked", min(fivePrimeInsertion, fivePrimeShear), NA))))
    sense <- sense %>% 
      mutate(threePrimeGap = ifelse(Flanking == "Dually-Flanked", max(threePrimeShear, threePrimeInsertion), 
                                ifelse(Flanking == "3'-Flanked", max(threePrimeShear, threePrimeInsertion), 
                                       ifelse(Flanking == "5'-Flanked", NA, "error"))))
    
    sensefivegapcount <- length(sense$fivePrimeGap[!is.na(sense$fivePrimeGap)])
    sensethreegapcount <- length(sense$threePrimeGap[!is.na(sense$threePrimeGap)])
    senseLflankcount <- sum(sense$Flanking == "5'-Flanked")
    senseRflankcount <- sum(sense$Flanking == "3'-Flanked")
    senseBflankcount <- sum(sense$Flanking == "Dually-Flanked")
    
    cat("There are", sensefivegapcount, "and", sensethreegapcount, "non-NAs for the 5' and 3' insertion coordinates in the sense orientation, respectively.\n")
    cat("Of these,", senseLflankcount, "are flanked only at the 5' end,", senseRflankcount, "only at the 3' end,", "and", senseBflankcount, "are flanked on both ends.\n\n")
    
    #Gap Coordinates are the coordinates corresponding to the shear sites -- a new columm was created solely for clarity and to minimize having to ifelse statements given orientation and flanking.
    antisense <- antisense %>%
      mutate(fivePrimeGap = ifelse(Flanking == "Dually-Flanked", max(fivePrimeInsertion, fivePrimeShear), 
                               ifelse(Flanking == "3'-Flanked", NA, 
                                      ifelse(Flanking == "5'-Flanked", max(fivePrimeInsertion, fivePrimeShear), "error"))))
    antisense <- antisense %>% 
      mutate(threePrimeGap = ifelse(Flanking == "Dually-Flanked", min(threePrimeShear, threePrimeInsertion), 
                                ifelse(Flanking == "3'-Flanked", min(threePrimeShear, threePrimeInsertion), 
                                       ifelse(Flanking == "5'-Flanked", NA, "error"))))
    antisensefivegapcount <- length(antisense$fivePrimeGap[!is.na(antisense$fivePrimeGap)])
    antisensethreegapcount <- length(antisense$threePrimeGap[!is.na(antisense$threePrimeGap)])
    antisenseLflankcount <- sum(antisense$Flanking == "5'-Flanked")
    antisenseRflankcount <- sum(antisense$Flanking == "3'-Flanked")
    antisenseBflankcount <- sum(antisense$Flanking == "Dually-Flanked")
    cat("There are", antisensefivegapcount, "and", antisensethreegapcount, "non-NAs for the 5' and 3' insertion coordinates in the antisense orientation, respectively.\n")
    cat("Of these,", antisenseLflankcount, "are flanked only at the 5' end,", antisenseRflankcount, "only at the 3' end,", "and", antisenseBflankcount, "are flanked on both ends.\n\n")
    
    
    
    dataframe <- rbind(sense, antisense, nonsense) 
    ##View(dataframe) here for ease of double-checking
    dataframefivegapcount <- length(dataframe$fivePrimeGap[!is.na(dataframe$fivePrimeGap)])
    dataframethreegapcount <- length(dataframe$threePrimeGap[!is.na(dataframe$threePrimeGap)])
    cat("There are", dataframefivegapcount, "and", dataframethreegapcount, "non-NAs for the 5' and 3' insertion coordinates in the dataframe orientation, respectively.\n")
    cat("Exiting qualCreate, returning to cleanReads...\n")
    
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(dataframe) 
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  cleanReads <- function(reads, tolerance, eFlag) {
    tic("cleanReads")
     ################################################################
     ## Cleans up the reads and prepares them for further analysis ##
     ##     STEP 1: PARTITION INTO LEFT, RIGHT, BOTH FLANKING      ##
     ##     STEP 2: CONVERT TO SENSE ORIENTATION                   ## 
     ################################################################
     cat("*************************************************\n")
     cat("____________________________Entering cleanReads..\n")
     cat("*************************************************\n")
     ##########################
     ## Step 1: Partitioning ##
     ##########################
     cat("CLEANING STEP 1: PARTITION INTO LEFT, RIGHT, BOTH FLANKING\n")
     if (is.null(reads) == TRUE) {
       cat("No data. Aborting.\n")
       quit()
     }
     Cleaned.reads <- reads %>% Partition %>% Concat()
     Discarded.reads <- rbind(Cleaned.reads$READ[Cleaned.reads$Run_ID == "BLANK_ROW"], Cleaned.reads$READ[Cleaned.reads$Flanking == "Non-Flanked"] )
     
     cat("Discarded Reads are: \n" )
     for (row in Discarded.reads) {
       cat(row, "\n")
     }
     cat("These reads are suspect because they lack flanks.\n While included for completeness of data, they are not recommended for use in downstream analysis. \n")
     
     
     cat("Done with Partitioning and concatenation\n")
     readlength <- nrow(reads)
     cleanlength <- nrow(Cleaned.reads)
     cat("Original number of reads:", readlength, "\n")
     cat("Number of rows after first cleaning step:", cleanlength, "\n")
     cat("\n")
     cat("Number of nonsensical and/or empty rows:", length(Discarded.reads))
     cat("\n")
      
     Cleaned.reads$chromosome<-factor(Cleaned.reads$chromosome, levels=chr0.order)
     if(eFlag) { 
        Cleaned.reads <- Excise(Cleaned.reads, tolerance)  ## Remove Noise
     }
    
     ####################################################
     ## Creates Gaps, preliminary orientation of gaps. ##
     ####################################################
     Cleaned.reads <- Cleaned.reads %>% rowwise() 
     Cleaned.reads <- qualCreate(Cleaned.reads)
    
     ######################################################
     ## STEP 2:Converts everything to sense orientation  ##
     ######################################################
     cat("Converting to sense orientation. Original designation is maintained, but conversion occurs for ease of analysis.\n")
     Cleaned.reads <- transform(Cleaned.reads,
                                threePrimeGap = ifelse(Flanking == "Dually-Flanked", pmax(threePrimeGap, fivePrimeGap), threePrimeGap),
                                fivePrimeGap = ifelse(Flanking == "Dually-Flanked", pmin(threePrimeGap, fivePrimeGap), fivePrimeGap)
     )
     cat("Is everything in the sense orientation?\n")
     if (all(Cleaned.reads$fivePrimeGap[Cleaned.reads$Flanking == "Dually-Flanked"] <= Cleaned.reads$threePrimeGap[Cleaned.reads$Flanking == "Dually-Flanked"])) {
       cat("Conditions met -- directionality is standardized.\n")
     } else {
       cat("Directionality not entirely standardized -- something has gone awry.\n")
       cat("Please check input.\n")
       cat("Aborting.\n")
       quit()
     } 
     
     statusCheck(Cleaned.reads)
     
     #goodCount <- length(Cleaned.Reads$READ[Cleaned.reads$passQualityCheck == TRUE])
     #badCount <- length(Cleaned.Reads$READ[Cleaned.reads$passQualityCheck == FALSE])
     #percentFailed <- round( x = (badCount/goodCount), digits = 2)
     #totalCount <- goodCount + badCount
     
     #cat("At this point in time", percentFailed, "reads of", totalCount, "have failed. \n")
     
     elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
     elapsed_time
     return(Cleaned.reads)  
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  findMappingErrors <- function(Cleaned.reads) { 
    tic("findMappingErrors")
    ####################################################################################################################################
    ## Assesses whether dually-flanked reads have reads that overlap to any degree, which would be a nonsensical result.              ##
    ## Overlapping insertion sites would leave no room for integration to occur and contrast with proposed mechanisms of integration. ##
    ####################################################################################################################################
    cat("*************************************************\n")
    cat("____________________Entering findMappingErrors...\n")
    cat("*************************************************\n\n")
    
    cat("This identifies reads that have seemingly mistaken mapping. This is typified with insertion sites exceed the range that they should\n (entering into the other flanks) and flank coordinates that are identical between shear and insert locations.\n\n\n")
    nonBothReads <- subset(Cleaned.reads, !Flanking == "Dually-Flanked")
    bothReads <- subset(Cleaned.reads, Flanking == "Dually-Flanked")
    cat("checking bothReads.\n")
    if (nrow(bothReads) > 0) {
      bothReads$failedStep <- NA}
    cat("checking nonBothReads.\n")
    if (nrow(nonBothReads) > 0) {
      nonBothReads$failedStep <- NA}    
    
    nrow(Cleaned.reads) == (nrow(nonBothReads) + nrow(bothReads))
   
    cat("Checking dually-flanked.\n")
    for (i in seq_along(bothReads)) {
    if (!is.na(bothReads$Human.Gene.Directionality[i])) {
      if (bothReads$Human.Gene.Directionality[i] == "Sense") {
        if (((bothReads$fivePrimeShear[i] <= bothReads$fivePrimeInsertion[i]) & (bothReads$fivePrimeInsertion[i] > bothReads$threePrimeInsertion[i])) |
            (bothReads$fivePrimeShear[i] == bothReads$fivePrimeInsertion[i]) |
            (bothReads$threePrimeShear[i] == bothReads$threePrimeInsertion[i])) {
          bothReads$passQualityCheck[i] <- FALSE
          bothReads$failedStep[i] <- "findMappingErrors"
        }
      } else if (bothReads$Human.Gene.Directionality[i] == "Antisense") {
        if (((bothReads$threePrimeShear[i] <= bothReads$threePrimeInsertion[i]) & (bothReads$threePrimeInsertion[i] > bothReads$fivePrimeInsertion[i])) |
            (bothReads$fivePrimeShear[i] == bothReads$fivePrimeInsertion[i]) |
            (bothReads$threePrimeShear[i] == bothReads$threePrimeInsertion[i])) {
          bothReads$passQualityCheck[i] <- FALSE
          bothReads$failedStep[i] <- "findMappingErrors"
        }
      }
    }
  }
    
    cat("Double-checking. \n")
    for (i in seq_along(nonBothReads)) {
      if (i <= nrow(nonBothReads)) {  # Ensure i doesn't exceed the number of rows
        if ((is.na(nonBothReads$fivePrimeShear[i]) & !is.na(nonBothReads$fivePrimeInsertion[i])) ||
            (is.na(nonBothReads$threePrimeShear[i]) & !is.na(nonBothReads$threePrimeInsertion[i]))) {
          nonBothReads$passQualityCheck[i] <- FALSE
          nonBothReads$failedStep[i] <- "findMappingErrors"
        } else if ((!is.na(nonBothReads$fivePrimeShear[i]) && !is.na(nonBothReads$fivePrimeInsertion[i]) &&
                nonBothReads$fivePrimeShear[i] == nonBothReads$fivePrimeInsertion[i]) ||
               (!is.na(nonBothReads$threePrimeShear[i]) && !is.na(nonBothReads$threePrimeInsertion[i]) &&
                nonBothReads$threePrimeShear[i] == nonBothReads$threePrimeInsertion[i])) {
	nonBothReads$passQualityCheck[i] <- FALSE
     	nonBothReads$failedStep[i] <- "findMappingErrors"
	  }    
	} else {
        cat("Error: row count exceeds the number of rows.\n")
        break  # Exit the loop if i exceeds the number of rows
      	  }
    }
    
    
    failedbothreads <- subset(bothReads, bothReads$failedStep == "findMappingErrors")
    failednonbothreads <- subset(nonBothReads, nonBothReads$failedStep == "findMappingErrors")
    failedmapping <- rbind(failednonbothreads, failedbothreads)
    failedmapping <- failedmapping$READ
    failedmappingcount <- length(failedmapping)
    
    
    cat("Discarded Reads are: \n" )
    for (row in failedmapping) {
      cat(row, "\n")
    }
    cat("These are the reads that have nonsensical coordinate mapping. These are not included in downstream analysis, but are maintained in the dataset for completeness' sake.\n")
    
  
    Cleaned.reads <- rbind(bothReads, nonBothReads)
    cleancount <- nrow(Cleaned.reads)
    cat("In context, this means that", failedmappingcount, "reads of", cleancount, "are mapped in questionable way and have been excluded." )
    
    statusCheck(Cleaned.reads)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(Cleaned.reads)
  }
    
    
  #--------------------------------------------------------------------------------------------------------
  findGapDifferences <- function(Cleaned.reads, insert.refSize) { 
    tic("findGapDifferences")
    #########################################################################################
    ## Defines difference as the difference between the gaps' extremes (if Dually-Flanked),  ##
    ## difference between the gap positions + insert reference (if uni-flanked)            ##
    #########################################################################################
    cat("*************************************************\n")
    cat("___________________Entering findGapDifferences...\n")
    cat("*************************************************\n")
    #totalReadLength should correspond to the total read length of the sequenced DNA -- flank(s) + insert.
    
    validReads <- subset(Cleaned.reads, passQualityCheck == TRUE)
    invalidReads <- subset(Cleaned.reads, passQualityCheck == FALSE)
    invalidReads <- invalidReads %>% mutate(slippage = NA, totalReadLength = NA)
    
    Cleaned.reads <- validReads %>% 
      mutate(totalReadLength = ifelse(Flanking == "Dually-Flanked", (abs(validReads$threePrimeGap - validReads$fivePrimeGap) + validReads$INSERT_LEN), 
                                      ifelse(Flanking == "3'-Flanked" & Human.Gene.Directionality == "Sense", (abs(validReads$threePrimeShear - validReads$threePrimeInsertion) + validReads$INSERT_LEN), 
                                             ifelse(Flanking == "5'-Flanked" & Human.Gene.Directionality == "Sense", abs(validReads$fivePrimeShear - validReads$fivePrimeInsertion) + validReads$INSERT_LEN, 
                                                    ifelse(Flanking == "3'-Flanked" & Human.Gene.Directionality == "Antisense", (abs(validReads$threePrimeInsertion - validReads$threePrimeShear) + validReads$INSERT_LEN),
                                                           ifelse(Flanking == "5'-Flanked" & Human.Gene.Directionality == "Antisense", (abs(validReads$fivePrimeShear - validReads$fivePrimeInsertion) + validReads$INSERT_LEN), "error"))))))
                  
    Cleaned.reads$totalReadLength <- as.numeric(Cleaned.reads$totalReadLength)
    Cleaned.reads <- Cleaned.reads %>% 
      mutate(slippage = ifelse(Flanking == "Dually-Flanked", abs(Cleaned.reads$threePrimeInsertion - Cleaned.reads$fivePrimeInsertion), NA))
    
    tooslipped <- nrow(Cleaned.reads[Cleaned.reads$slippage > maxPCRSlippage & Cleaned.reads$Flanking == "Dually-Flanked", ])
    cat("There are", tooslipped, "dually-flanked reads that have insertion points beyond the designated 10 base pair threshold... \n")
    cat("This may be because of mismapping.\n")
    concerningreads <- Cleaned.reads$READ[Cleaned.reads$slippage > maxPCRSlippage & Cleaned.reads$Flanking == "Dually-Flanked"]
    
    cat("\nRead IDs of concern are: \n \n" )
    for (row in concerningreads) {
      cat(row, "\n")
    }
    cat("These reads are suspect because integration coordinates are too far away. While included for completeness of data, they are not recommended for use in downstream analysis. \n")
    
    cat("Updating quality check status...\t")
    Cleaned.reads$passQualityCheck[(Cleaned.reads$READ %in% concerningreads)] <- FALSE
    Cleaned.reads$failedStep[(Cleaned.reads$READ %in% concerningreads)] <- "findGapDifferences"
    
    cat(colnames(validReads))
    cat(colnames(invalidReads))
    
    cat("Adding flank spans...\n")
    Cleaned.reads <- rbind(Cleaned.reads, invalidReads)
    Cleaned.reads <- Cleaned.reads %>% 
      mutate(fivePrime_span = abs(Cleaned.reads$fivePrimeInsertion - Cleaned.reads$fivePrimeShear)) %>% 
      mutate(threePrime_span = abs(Cleaned.reads$threePrimeShear - Cleaned.reads$threePrimeInsertion))
    
    statusCheck(Cleaned.reads)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(Cleaned.reads)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  getPCRDuplicates <- function(dataframe, flanking) {
  
    
    #####################################################################################################################
    #####################################################################################################################
    ##The function that actually quantifies how many PCR Duplicates there are in a subsetted dataframe.                ##
    ##After flanking is pre-set, it is used to characterize which grouping parameters are used within the function.    ##
    ##Grouping then occurs such that a sum and unique ID can be mutated into the appropriate column.                   ##
    ##Called by markPCRDuplicates.                                                                                     ##
    #####################################################################################################################
    #####################################################################################################################
    
    cat("*************************************************\n")
    cat("_______________________Entering getPCRDuplicates.\n")
    cat("*************************************************\n\n")
    
    name_df <- deparse(substitute(dataframe))
    name_flanking <- deparse(substitute(flanking))
    
    cat("Counting PCR duplicates for the reads of", name_df, "to quantify noise generated in sequencing preparation....\n")
    cat("ReadLength is not being uses as a metric to determine PCR duplicates right now.\n")
    
    config_both <- list(
      groupingColumns = c("fivePrimeShear", "fivePrimeInsertion", "threePrimeInsertion", "threePrimeShear", "totalReadLength")
    )
    config_left <- list(
      groupingColumns = c("fivePrimeShear", "fivePrimeInsertion", "totalReadLength")
    )
    config_right <- list(
      groupingColumns = c("threePrimeInsertion", "threePrimeShear", "totalReadLength")
    )
    config_non <- list(
      groupingColumns = c("INSERT_LEN")
    )  
    
  
    
    if (flanking == "Dually-Flanked") {
      config <- config_both
      cat("Working on dually-flanked reads.\n")
    } else if (flanking == "5'-Flanked") {
      config <- config_left
      cat("Working on 5'-flanked reads.\n")
    } else if (flanking == "3'-Flanked"){
      config <- config_right
      cat("Working on 3'-flanked reads.\n")
    }
      else if (flanking == "Non-Flanked") {
      cat("These are non-flanked reads.\n")
      cat("Because of the general dearth of data for this subset, assumptions about matching values in INSERT_LEN\n are being used to characterize reads as PCR duplicates. \n")
      config <- config_non
    } else {
      stop("Invalid 'flanking' value.")
    }
    
    groupedReads <- dataframe %>% 
      group_by(across(all_of(config$groupingColumns))) %>% 
      mutate(PCR.Dups = dplyr::n(), PCR.ID = random_id(bytes = 5))
  
    
    cat("Exiting getPCRDuplicates for the", name_df, "subset of the broader Cleaned.reads dataframe.\n\n")  
    return(groupedReads)
    
    
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
   markPCRDuplicates <- function(Cleaned.reads) {
    
     ########################################
     ## Find reads that are PCR duplicates ##
     ########################################
     tic("markPCRDuplicates")
     cat("*************************************************\n")
     cat("____________________Entering markPCRDuplicates...\n")
     cat("*************************************************\n")
     Cleaned.reads$PCR.Dups <- 1
     
     cat("Splitting Cleaned.reads into two dataframes based on whether they have passed the quality check or not.\n")
     cat("Respectively, validReads or invalidReads.\n")
     
     startCount <- nrow(Cleaned.reads)
     
     cat("Converting antisense reads to sense directionality to ensure maximum PCR Duplicate capture.\n")
     Cleaned.reads[which(Cleaned.reads$Human.Gene.Directionality == "Antisense"), 
          c("fivePrimeShear", "fivePrimeInsertion", "threePrimeInsertion", "threePrimeShear")] <- 
       rev(Cleaned.reads[which(Cleaned.reads$Human.Gene.Directionality == "Antisense"), 
                 c("fivePrimeShear", "fivePrimeInsertion", "threePrimeInsertion", "threePrimeShear")])
     
     validReads <- subset(Cleaned.reads, passQualityCheck == TRUE)
     invalidReads <- subset(Cleaned.reads, passQualityCheck == FALSE)
  
     ## Merge the reads if their chromosome, fivePrimeShear/insertion, and threePrimeInsertion/shear are the same ##
     ## These reads are likely PCR duplicates                                                         ##
     #Cleaned.reads <- Cleaned.reads %>% group_by(fivePrimeShear, fivePrimeInsertion, threePrimeInsertion, threePrimeShear) %>% 
                                        #mutate(PCR.Dups = sum(PCR.Dups)) #%>% 
                                        #distinct(chromosome, fivePrimeShear, fivePrimeInsertion, 
                                                 #threePrimeInsertion, threePrimeShear, .keep_all = TRUE) #might be destroying data
     
     
     cat("This is for ease of quality control. The following numbers should be used as reference,\nand they should correspond to the numbers that follow. \n\n")
     bothvalnumber <- nrow(validReads[validReads$Flanking == "Dually-Flanked", ])
     rightvalnumber <- nrow(validReads[validReads$Flanking == "3'-Flanked", ])
     leftvalnumber <- nrow(validReads[validReads$Flanking == "5'-Flanked", ])
     bothinvalnumber <- nrow(invalidReads[invalidReads$Flanking == "Dually-Flanked", ])
     rightinvalnumber <- nrow(invalidReads[invalidReads$Flanking == "3'-Flanked", ])
     leftinvalnumber <- nrow(invalidReads[invalidReads$Flanking == "5'-Flanked", ])
     noninvalnumber <- nrow(invalidReads[invalidReads$Flanking == "Non-Flanked", ])
     usefulnumber <- nrow(validReads)
     uselessnumber <- nrow(invalidReads)
     
     cat("There were", startCount, "reads at the start of the analysis. Of those,\n\t", usefulnumber, "remain viable for downstream analysis, while", uselessnumber, "are already being excluded.\n")
     
     cat("Of those that are still viable", bothvalnumber, "are reads flanked on both sides while", bothinvalnumber, "are flanked on both sides on the inviable pool.\n")
     cat("Of those", leftvalnumber, "are reads flanked only on the 5' end while", leftinvalnumber,"are 5' flanked in the inviable pool.\n")
     cat("Of those", rightvalnumber, "are reads flanked only on the 3' end while", rightinvalnumber, "are 3' flanked in the inviable pool.\n\n")
     
     
     validReadsDualFlanked <- validReads %>%
       dplyr::filter(Flanking == "Dually-Flanked") 
     validReads5Flanked <- validReads %>%
       dplyr::filter((Flanking == "5'-Flanked" & Human.Gene.Directionality == "Sense") | 
                       (Flanking == "3'-Flanked" & Human.Gene.Directionality == "Antisense")) #(sense and 5') or antisense and 3'
     validReads3Flanked <- validReads %>%
       dplyr::filter((Flanking == "5'-Flanked" & Human.Gene.Directionality == "Antisense") | 
                       (Flanking == "3'-Flanked" & Human.Gene.Directionality == "Sense")) 
     
     vBF <- getPCRDuplicates(validReadsDualFlanked, "Dually-Flanked")
     validIntermediateCount <- vBF 
     v5F <- getPCRDuplicates(validReads5Flanked, "5'-Flanked")
     validIntermediateCount <- rbind(validIntermediateCount, v5F)
     v3F <- getPCRDuplicates(validReads3Flanked, "3'-Flanked")
     validIntermediateCount <- rbind(validIntermediateCount, v3F)
     
     cat("Containing only the valid reads, PCR duplicates have been quantified for", nrow(validIntermediateCount), "reads.\n")
     cat("This compares to the", nrow(validReads), "reads of all validReads.\n\n")
     
     invalidReadsDualFlanked <- invalidReads %>%
       dplyr::filter(Flanking == "Dually-Flanked") 
     invalidReads5Flanked <- invalidReads %>%
       dplyr::filter((Flanking == "5'-Flanked" & Human.Gene.Directionality == "Sense") | 
                       (Flanking == "3'-Flanked" & Human.Gene.Directionality == "Antisense"))
     invalidReads3Flanked <- invalidReads %>%
       dplyr::filter((Flanking == "5'-Flanked" & Human.Gene.Directionality == "Antisense") | 
                       (Flanking == "3'-Flanked" & Human.Gene.Directionality == "Sense")) 
     invalidReadsNonFlanked <- invalidReads %>%
       dplyr::filter(Flanking == "Non-Flanked") 
     
     iBF <- getPCRDuplicates(invalidReadsDualFlanked, "Dually-Flanked")
     invalidIntermediateCount <- iBF
     i5F <- getPCRDuplicates(invalidReads5Flanked, "5'-Flanked")
     invalidIntermediateCount <- rbind(invalidIntermediateCount, i5F)
     i3F <- getPCRDuplicates(invalidReads3Flanked, "3'-Flanked")
     invalidIntermediateCount <- rbind(invalidIntermediateCount, i3F)
     iNF <- getPCRDuplicates(invalidReadsNonFlanked, "Non-Flanked")
     invalidIntermediateCount <- rbind(invalidIntermediateCount, iNF)
     
     cat("Containing only the invalid reads, clones have been quantified for", nrow(invalidIntermediateCount), "reads.\n")
     cat("This compares to the", nrow(invalidReads), "reads of all validReads.\n\n")
     
     
     
     #both_reads <- validReads %>%
      # filter(Flanking == "Dually-Flanked") %>% 
      # group_by(fivePrimeShear, fivePrimeInsertion, threePrimeInsertion, threePrimeShear, totalReadLength) %>% 
      # mutate(PCR.Dups = sum(PCR.Dups))
     #bothpcrdupcount <- nrow(both_reads[both_reads$PCR.Dups >1, ])
     #cat("Of", nrow(both_reads), "dually-flanked reads, there are", bothpcrdupcount, "reads associated with greater than 1 PCR duplicates.\n")
     #both_reads <- both_reads %>% ungroup()
     
     #left_reads <- validReads %>%
      # filter(Flanking == "5'-Flanked") %>% 
      # group_by(fivePrimeShear, fivePrimeInsertion, totalReadLength) %>% 
      # mutate(PCR.Dups = sum(PCR.Dups))
     #leftpcrdupcount <- nrow(left_reads[left_reads$PCR.Dups >1, ])
     #cat("Of", nrow(left_reads), "5'-flanked reads, there are", leftpcrdupcount, "reads associated with greater than 1 PCR duplicates.\n")
     #left_reads <- left_reads %>% ungroup()
     
     #right_reads <- validReads %>%
      # filter(Flanking == "3'-Flanked") %>% 
      # group_by(threePrimeInsertion, threePrimeShear, totalReadLength) %>% 
      # mutate(PCR.Dups = sum(PCR.Dups))
     #rightpcrdupcount <- nrow(right_reads[right_reads$PCR.Dups >1, ])
     #cat("Of", nrow(right_reads), "3'-flanked reads, there are", rightpcrdupcount, "reads associated with greater than 1 PCR duplicates.\n")
     #right_reads <- right_reads %>% ungroup()
     
     #non_reads <- validReads %>%
      # filter(Flanking == "Non-Flanked") 
     
     totalreadcount <- nrow(invalidIntermediateCount) + nrow(validIntermediateCount) 
     cat("There are", totalreadcount, "reads accounted for.\n")
     
     Cleaned.reads <- rbind(validIntermediateCount, invalidIntermediateCount)
     cat("Reverting antisense reads to antisense directionality.\n")
     Cleaned.reads[which(Cleaned.reads$Human.Gene.Directionality == "Antisense"), 
                   c("fivePrimeShear", "fivePrimeInsertion", "threePrimeInsertion", "threePrimeShear")] <- 
       rev(Cleaned.reads[which(Cleaned.reads$Human.Gene.Directionality == "Antisense"), 
                         c("fivePrimeShear", "fivePrimeInsertion", "threePrimeInsertion", "threePrimeShear")])
     statusCheck(Cleaned.reads)
     elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
     elapsed_time
     return(Cleaned.reads)
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  writePhyloMasterFrame <- function(Cleaned.reads) {
    tic("writePhyloMasterFrame")
    ###########################################
    ## Creates the PhyloMasterFrame.csv file ##
    ###########################################
    cat("*************************************************\n")
    cat("________________Entering writePhyloMasterFrame...\n")
    cat("*************************************************\n")
    ## The following line is not used, but kept in ##
    PhyloWClones <- Cleaned.reads %>% group_by(fivePrimeInsertion, threePrimeInsertion) %>% 
      distinct(fivePrimeInsertion, threePrimeInsertion, .keep_all = TRUE) 
    
    #write.xlsx(x = as.data.frame(Cleaned.reads), file = "PhyloMasterFrame.xlsx", col.names = TRUE, row.names = T, showNA =  FALSE)
    write.table(PhyloWClones, file = "PhyloMasterFrame.txt", sep = "\t", row.names = FALSE)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  removeBadReads <- function(Cleaned.reads, maxDiff, minDiff) { 
    tic("removeBadReads")
    cat("*************************************************\n")
    cat("_________________________Entering removeBadReads.\n")
    cat("*************************************************\n\n")
    
    cat("The current maximal limit, set in the code as hostTolerance, is", hostTolerance, "bases.\nAnything above this will be removed from downstream analysis.\n")
    cat("Likewise, the current minimal limit, set in the code as minReadLength, is", minReadLength, "bases.\nAnything below this is unlikely to have mapped correctly; PacBio instrumentation does not yet support readlengths this small.\n")
    cleancount <- nrow(Cleaned.reads)
    cat("removeBadReads is starting with an input of", cleancount, "reads.\n")
  
    
    validReads <- subset(Cleaned.reads, passQualityCheck == TRUE)
    invalidReads <- subset(Cleaned.reads, passQualityCheck == FALSE)
    ##View(failedReads) here for ease of double-checking
    ##View(validReads)
    validReadsCorrected <- validReads[(validReads$totalReadLength <= maxDiff & validReads$totalReadLength >= minDiff), ]
    ##View(validReadsCorrected)
    
    difference <- nrow(validReads) - nrow(validReadsCorrected)
    
    validReads <- subset(Cleaned.reads, passQualityCheck == TRUE)
    ##View(validReadsCorrected) here for ease of double-checking
    badReads <- validReads[!(validReads$totalReadLength <= maxDiff & validReads$totalReadLength >= minDiff), ]#appending quality check info
    badReads$passQualityCheck <- FALSE ####################################### CREATING WEIRD NA ROWS (due to antisense not getting a totalreadlength?)
    ##View(badReads)
    badReads$failedStep <- "removeBadReads"
    
    #getting counts for report.
    highCount <- nrow(badReads$READ[badReads$totalReadLength > maxDiff])
    lowCount <- nrow(badReads$READ[badReads$totalReadLength < minDiff])
    
    cat("\nRead IDs of concern are: \n \n" )
    for (row in badReads$READ) {
      cat(row, "\n")
    }
    cat("\tAgain, these are of dubious quality because they are too long or short to be valid. These are still included in the output, but not used for downstream analysis.\n\n")
    ##View(badReads) here for ease of double-checking
    badcount <- nrow(badReads)
    
    
    ########################################################################
    ## REMOVE THE FOLLOWING READS -- THE AREA IS GREATER THAN 1MB IN SIZE ##
    ## SHOULD WRITE A FUNCTION TO CALCULATE THESE AT SOME POINT IN TIME   ##
    ########################################################################
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/133937/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/2557018/ccs")   
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/2688601/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/7472827/ccs")   
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/8127194/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/17434289/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/17827328/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/19005922/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/19269654/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/21234022/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/30083085/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/36571170/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/41092042/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/47385334/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/68225388/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/120259137/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/132186983/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/136383206/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/147065410/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/165216958/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/31720860/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/42469112/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/42732288/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/55837869/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/89393647/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/127468094/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/138545946/ccs")
    Cleaned.reads <- subset(Cleaned.reads, READ != "m64152e_220512_221341/150800061/ccs")
    
    Cleaned.reads <- rbind(validReadsCorrected, invalidReads, badReads)
    cleancount <- nrow(Cleaned.reads)
    failcount <- length(Cleaned.reads$READ[Cleaned.reads$passQualityCheck == FALSE]) #Inelegant, but the number is the same. 
    usefulnumber <- cleancount - failcount
    
    cat("There were", cleancount, "reads at the start of the analysis. Of those,", usefulnumber, "remain viable for downstream analysis.\n")
    cat(failcount, "are the number of reads that have dubious quality. These are maintained in the dataframe for completeness, but excluded for further analysis.\n")
    
    
    cat("There are", highCount, "reads that had lengths greater than", hostTolerance, "bases.\n Likewise, there were", lowCount, "reads with lenths less than", minReadLength, "bases.\n", "This means that there are", badcount, "reads above this threshold or below", minDiff, "bases. \n")
    statusCheck(Cleaned.reads)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    
    return(Cleaned.reads)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  getClones <- function(dataframe, flanking) {
    
    
    ##############################################################################################################################################
    ##############################################################################################################################################
    ##The function that actually quantifies how many clones there are in a Cleaned.reads subsetted dataframe.                                   ##
    ##Essentially, the dataframe is copied, and the copy is de-duplicated. This is because PCR duplicates always added on to clone counts.      ##
    ##De-duplicating could not be done by hierarchical group_by use because you cannot subgorup by the same variables in a sensical way.        ##
    ##De-duplicated data is then grouped, counted, and ungrouped. These data are then left-joined with the original dataframe (not de-dupped).  ##
    ##based on insertion coordinates, with the same insertion coordinates being given the same clone counts.                                    ##
    ##Called by markClonalExpansion.                                                                                                            ##
    ##############################################################################################################################################
    ##############################################################################################################################################
    
    cat("*************************************************\n")
    cat("______________________________Entering getClones.\n")
    cat("*************************************************\n\n")
    
    name_df <- deparse(substitute(dataframe))
    name_flanking <- deparse(substitute(flanking))
    
    cat("Counting Clones for the reads of", name_df, "after identifying unique reads....\n")
    cat("totalReadLength not currently being used to identify and separate clones.\n")
    
      
    config_both <- list(
      distinctColumns = c("fivePrimeShear", "fivePrimeInsertion", "threePrimeInsertion", "threePrimeShear", "totalReadLength"),
      groupByColumns = c("fivePrimeInsertion", "threePrimeInsertion"),
      directionReference = c("threePrimeInsertion")
    )
    
    config_left <- list(
      distinctColumns = c("fivePrimeShear", "fivePrimeInsertion", "totalReadLength"),
      groupByColumns = c("fivePrimeInsertion"),
      directionReference = c("fivePrimeInsertion")
    )
    
    config_right <- list(
      distinctColumns = c("threePrimeInsertion", "threePrimeShear", "totalReadLength"),
      groupByColumns = c("threePrimeInsertion"),
      directionReference = c("threePrimeInsertion")
    )
    
    #Conditionally sets what is done based on input.
    if (flanking == "Dually-Flanked") {
      config <- config_both
      cat("Working on dually-flanked reads.\n")
    } else if (flanking == "5'-Flanked") {
      config <- config_left
      cat("Working on 5'-flanked reads.\n")
    } else if (flanking == "3'-Flanked"){
      config <- config_right
      cat("Working on 3'-flanked reads.\n")
    }
      else if (flanking == "Non-Flanked") {
        cat("These are non-flanked reads.\n")
        cat("It is impossible to determine clonality from these data. Assuming broadly and assigning each read as a unique clone.\n")
        dataframe$Clones <- 1
        groupedReads <- dataframe %>% 
          group_by(INSERT_LEN) %>% 
          mutate(Clone.ID = random_id(bytes = 6))
        return(groupedReads)
        quit()
    } else {
    stop("Invalid 'flanking' value.")
    }
    
    readsDedup <- dataframe %>%
      distinct(across(all_of(config$distinctColumns))) %>% #, ~ (.keep_all = TRUE))) %>%  ##~ (.keep_all = TRUE)))  %>% #created problems after update
      mutate(Clones = 1) %>%  # site of errors if clones don't get reported accurately
      group_by(across(all_of(config$groupByColumns))) %>%
      mutate(Clones = sum(Clones), Clone.ID = random_id(bytes = 7))
    
    
    #used for reporting the number  
    CloneSum <- readsDedup %>%  
      group_by(across(all_of(config$groupByColumns))) %>% 
      summarize(Clones = sum(Clones)/dplyr::n()) 
      ##View(bothCloneSum) ###here for ease of double checking
      CloneSum <- sum(CloneSum$Clones) 
      cat("There are", CloneSum, "counted clones for the", name_df, "subset of the broader Cleaned.reads dataframe.\n")
      
    sharedColumns <- intersect(names(dataframe), names(readsDedup))
    cat(sharedColumns, "\n\n")
    sharedColumns <- setdiff(sharedColumns, c("INDEX", "READ", "HUMAN_ALTS", "UNMAPPED", "PCR.Dups", "PCR_ID")) 
    cat("subtracted:", "\n", sharedColumns, "\n")
    readsIntermediate <- dataframe %>% 
      left_join(y = readsDedup, by = sharedColumns)
    ##View(readsIntermediate)   for quality control
    readsIntermediate <- readsIntermediate %>%
      dplyr::select(-ends_with(".y")) %>% #made select dplyr::select
      rename_at(vars(ends_with(".x")), ~gsub("\\.x$", "", .))
    ##View(readsIntermediate)  
    
    output <- readsIntermediate %>%
      arrange(desc(config$directionReference)) %>% 
      group_by(across(all_of(c(sharedColumns, "Clones"))))  %>%
      fill(Clones, Clone.ID, .direction = "downup") %>%
      ungroup()
    ##View(output) # for quality control
    
    cat("Exiting getClones for the", name_df, "subset of the broader Cleaned.reads dataframe.\n\n")  
    return(output)
  
  
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  
  findOddFails <- function(Cleaned.reads) {
    tic("findOddFails")
    ###################################################################################
    ## Finds odd scenario of high quality reads with no human flanks mapped.         ##
    ###################################################################################
    
    cat("******************************************\n")
    cat("__________________Entering findOddFails...\n")
    cat("****************************************\n\n")
    
    
    cat("Splitting Cleaned.reads into two dataframes based on whether they have passed the quality check or not.\n")
    cat("Respectively, validReads or invalidReads.\n")
    validReads <- subset(Cleaned.reads, passQualityCheck == TRUE)
    invalidReads <- subset(Cleaned.reads, passQualityCheck == FALSE)
    
    Pre.C.E.Count <- nrow(Cleaned.reads)
    cat("Prior to looking for flankless but high quality integrations, there were", Pre.C.E.Count, "reads.\n")
    
    stillValids <- validReads %>% 
      dplyr::filter(!is.na(chromosome))
    
    badReads <- validReads %>% 
      dplyr::filter(is.na(chromosome)) %>%
      mutate(passQualityCheck = FALSE) %>%
      mutate(failedStep = "findOddFails")
    
    Cleaned.reads <- rbind(stillValids, invalidReads, badReads)
    cleancount <- nrow(Cleaned.reads)
    failcount <- length(Cleaned.reads$READ[Cleaned.reads$passQualityCheck == FALSE]) #Inelegant, but the number is the same. 
    usefulnumber <- cleancount - failcount
    
    cat("There were", cleancount, "reads at the start of the analysis. Of those,", usefulnumber, "remain viable for downstream analysis.\n")
    cat(failcount, "are the number of reads that have dubious quality. These are maintained in the dataframe for completeness, but excluded for further analysis.\n")
    
    statusCheck(Cleaned.reads)
    
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(Cleaned.reads)
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  standardizeCloneIDs <- function(validIntermediateCount) {
    ####################################################################################
    ## Integrations that are the same but on different cells may shear differently.   ##
    ## This function identifies reads with different flanking patterns and calculates ##
    ## the appropriate clone count and assigns a new clone ID to these disparate reads##
    ## Called by markClonalExpansion.                                                 ##
    ####################################################################################
    
    if (nrow(validIntermediateCount) == 0) {
      return(validIntermediateCount)
    }
    
    cat("Entering standardizeCloneIDs...\n")
    cat("Identifying clonal expansion events that sheared and sequenced differently.\n")
    cat("Calculating updated clone counts. \n")
    identifier <- validIntermediateCount[c("Flanking", "chromosome", "PCR.ID", "Clone.ID", "Clones", "fivePrimeInsertion", "threePrimeInsertion")]
    
    identifier <- identifier %>%
      # this standardizes the clone number but not ID -- I need it to not do that yet -- if clones don't have amended clone counts, I can standardize directionality, then use reduce_ranges?
      distinct(Flanking, Clone.ID, chromosome, Clones, fivePrimeInsertion, threePrimeInsertion) %>%
      group_by(chromosome, fivePrimeInsertion) %>%
      mutate(
        IDsPer5 = ifelse(is.na(fivePrimeInsertion), 0, n_distinct(Clone.ID)),
        CloneIDsList = ifelse(is.na(fivePrimeInsertion), list(), list(unique(Clone.ID[!duplicated(Clone.ID)])))
      ) %>%
      group_by(chromosome, threePrimeInsertion) %>%
      mutate(
        IDsPer3 = ifelse(is.na(threePrimeInsertion), 0, n_distinct(Clone.ID)),
        CloneIDsList = ifelse(is.na(threePrimeInsertion), CloneIDsList, list(union(unique(CloneIDsList[[1]]), unique(Clone.ID))))
      ) %>% 
      group_by(chromosome, fivePrimeInsertion) %>%
      mutate(
        CloneIDsList = ifelse(is.na(fivePrimeInsertion), CloneIDsList, list(unique(unique(CloneIDsList[[1]]), unique(Clone.ID))))
      )
    
    
    ##Modified from v2.2.4 -- this accounts for the mapping artifacts.
    identifier$missingValue <- NA
    identifier <- identifier %>%
      rowwise() %>%
      mutate(
        missingValue = if (Flanking != "Dually-Flanked") {
          if (is.na(fivePrimeInsertion)) {"fivePrimeInsertion"} 
          else if (is.na(threePrimeInsertion)) {"threePrimeInsertion"} 
          else {NA}
        } else {NA}
      )  
    cat("Separating Dual-Flanked from Uni-Flanked.\n")
    onlyDual <- subset(identifier, Flanking == "Dually-Flanked")
    cat("Identifying and flipping any nonsensical reads in onlyDual.\n")
    nonflipped <- subset(onlyDual, onlyDual$fivePrimeInsertion <= onlyDual$threePrimeInsertion)
    nonflipped$flippedOrientation <- FALSE
    flipped <- subset(onlyDual, onlyDual$fivePrimeInsertion > onlyDual$threePrimeInsertion)
    flipped <- flipped %>% 
      mutate(flippedOrientation = TRUE) %>%
      dplyr::rename(fivePrimeInsertion = threePrimeInsertion, threePrimeInsertion = fivePrimeInsertion)
    cat("Standardizing direction in onlyDual.\n")
    onlyDual <- rbind(flipped, nonflipped)
    cat("Filling in missing values of Uni-Flanked -- needed for conversion into GRanges datatype.\n")
    withoutDual <- subset(identifier, Flanking != "Dually-Flanked")
    withoutDual <- withoutDual %>%
      mutate(
        fivePrimeInsertion = if_else(missingValue == "fivePrimeInsertion", (threePrimeInsertion-1), fivePrimeInsertion),
        threePrimeInsertion = if_else(missingValue == "threePrimeInsertion", (fivePrimeInsertion+1), threePrimeInsertion)
      ) 
    withoutDual$flippedOrientation <- FALSE
    #reClonedReducedRanges <- rbindlist(list(onlyDual, withoutDual), fill = TRUE) #this stopped working -- all columns were included
    identifier <- rbind(onlyDual, withoutDual) 
    identifier <- identifier %>%  
      dplyr::relocate(fivePrimeInsertion, .before = threePrimeInsertion) %>%
      dplyr::relocate(chromosome, fivePrimeInsertion, threePrimeInsertion, .before = Flanking)
    
    cat("Reformating done.\n")
    cat("Modifying list structure to rid of duplicates.\n")
    cat("Addressing overlapping but non-identical ranges.\n")
    cat("For example, if read 1 exists between X and Y, while read 2 exists between Y and Z  -- these are clearly the same clone but with mapping issues. \n")
    
    identifier <- makeGRangesFromDataFrame(identifier, seqnames.field =  "chromosome", start.field = "fivePrimeInsertion", end.field =  "threePrimeInsertion", keep.extra.columns = TRUE)
    #test <- group_by_overlaps(identifier, identifier) #ALTERNATIVE APPROACH
    #test <- as.data.frame(test) #ALTERNATIVE APPROACH
    #print(test$CloneIDsList.subject) other option! #ALTERNATIVE APPROACH
    
    if (length(reduce_ranges(identifier)) < length(identifier)) {
      numberOfClones <- length(reduce_ranges(identifier))
      numberOfUniqueCloneIDs <- length(unique(identifier$Clone.ID))
      if (numberOfClones < numberOfUniqueCloneIDs) {
        cat("There are clone IDs to concatenate. \n"
        )
        identifier <- reduce_ranges(identifier, CloneIDsList = c(Clone.ID)) #Works when data isn't simple
      }
      if (numberOfClones == numberOfUniqueCloneIDs) {
        cat("There are not clone IDs to concatenate.")
        random_ids <- sapply(1:numberOfClones, function(x) random_id(bytes = 7))
        identifier <- reduce_ranges(identifier, CloneIDsList = random_ids)
      }
      
    }
    
    removeDuplicates <- function(x) {
      if (is.vector(x)) {
        return(unique(x))
      } else {
        return(x)
      }
    }
    
    # Apply the function to each element of the column
    identifier$CloneIDsList <- lapply(identifier$CloneIDsList, removeDuplicates)
  
    cat("test 1")
    identifier <- as.data.frame(identifier)
    identifier <- identifier %>%
      unnest(CloneIDsList)
    cat("Grouping to unlist and categorize Clone.IDs that correspond to the same clone.\n")
    identifier <- identifier %>%
      group_by(seqnames, start, end, width, strand) %>%
      summarize(CloneIDsList = toString(unique(unlist(CloneIDsList))))
    cat("Creating identifier dataframe that allows for correspondance of temporary Clone.IDs to correct, permanent Clone.IDs.\n\n")
    identifier <- identifier %>%
      rowwise() %>%
      mutate(newCloneID = random_id(bytes = 7))
    
    identifier <- identifier[c("CloneIDsList", "newCloneID")]
    
    identifier <- identifier %>%
      separate_rows(CloneIDsList, sep = ", ") %>%
      mutate_all(~gsub(" ", "", .))
    
    cat("Joining by temporary Clone.ID to assign new Clone.ID.\n")
    #still need to do this, clean up below. Validate that data wasn't lost.
    
    cat("Creating the Clone Count reference dataframe that accounts for mapping artifacts. Then, re-keying it with the new Clone.IDs. \n")
    cloneCountDF <- validIntermediateCount %>%
      group_by(Clones, Clone.ID, PCR.ID) %>% #PCR.ID was added
      summarise(Clones = mean(Clones)) %>% #this was the thing that was changed 
      as.data.frame()
    
    cloneCountDF <- cloneCountDF %>%
      left_join(identifier, by = c("Clone.ID" = "CloneIDsList")) %>%
      mutate(Clone.ID = ifelse(is.na(newCloneID), Clone.ID, newCloneID)) %>%
      select(-newCloneID)
    
    cloneCountDF <- cloneCountDF %>%
      group_by(Clone.ID) %>% 
      summarise(updatedClones = mean(Clones)) %>%
      as.data.frame()
    
    cat("Updating intermediate with the correct Clone.IDs; joining with updated counts.\n")
    validIntermediate <- validIntermediateCount %>%
      left_join(identifier, by = c("Clone.ID" = "CloneIDsList")) %>%
      mutate(Clone.ID = ifelse(is.na(newCloneID), Clone.ID, newCloneID)) %>%
      select(-newCloneID)
    
    validIntermediate <- validIntermediate %>%
      left_join(cloneCountDF, by = c("Clone.ID")) %>%
      mutate(Clones = updatedClones) %>%
      select(-updatedClones)
    
    return(validIntermediate) #error earlier -- validintermediatecount was returned
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  markClonalExpansion <- function(Cleaned.reads) {
    tic("markClonalExpansion")
    ###################################################################################
    ## Collapse clonal expansion clones where the insert hits the exact same         ##
    ## begin and end locations in the host genome -- PCR duplicates have already     ##
    ## been filtered, defined by sequences that are the exact same from begin to end ##
    ###################################################################################
    
    cat("*************************************************\n")
    cat("__________________Entering markClonalExpansion...\n")
    cat("*************************************************\n\n")
    #Pre.C.E.Count <- nrow(Cleaned.reads)
    #Cleaned.reads <- Cleaned.reads %>% group_by(fivePrimeInsertion, threePrimeInsertion) %>% mutate(Clones = sum(Clones)) 
    
    cat("Splitting Cleaned.reads into two dataframes based on whether they have passed the quality check or not.\n")
    cat("Respectively, validReads or invalidReads.\n")
    
    #all of this is based on a swapping of antisense data to collapse clones at the same site but on reverse complement sequenced templates. 
    #because a join operation is done late betwer readsDedup and the input df, I need to separately swap the switched columns back for each?
    
    cat("By nature of PacBio Sequencing, either strand of DNA can be amplified. As a result, integrations at the same site might be reported from either DNA.\n")
    cat("Temporary conversion of of antisense reads to sense directionality allows for identification and quantification of reads that integrate at the same site \n")
    cat("but are detected by opposite strands.\n")
    #Cleaned.reads$flippedOrientation <- FALSE
    #Cleaned.reads[Cleaned.reads$Human.Gene.Directionality == "Antisense", c("fivePrimeShear", "fivePrimeInsertion", "threePrimeInsertion", "threePrimeShear")] <- 
    #  Cleaned.reads[Cleaned.reads$Human.Gene.Directionality == "Antisense", c("threePrimeShear", "threePrimeInsertion", "fivePrimeInsertion", "fivePrimeShear")]
    #Cleaned.reads$flippedOrientation <- ifelse(Cleaned.reads$Human.Gene.Directionality == "Antisense", TRUE, FALSE)
   
    
    swappedReads <- Cleaned.reads %>%
      mutate(
        flippedOrientation = FALSE,
        tempFivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeShear, NA),
        tempFivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeInsertion, NA),
        fivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", threePrimeShear, fivePrimeShear),
        fivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", threePrimeInsertion, fivePrimeInsertion),
        threePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeShear, threePrimeShear),
        threePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeInsertion, threePrimeInsertion),
        flippedOrientation = ifelse(Human.Gene.Directionality == "Antisense", TRUE, FALSE)
      ) %>%
      dplyr::select(-starts_with("temp")) %>%
      dplyr::ungroup()
    rownames(swappedReads) <- NULL
    
    validReads <- subset(swappedReads, passQualityCheck == TRUE)
    invalidReads <- subset(swappedReads, passQualityCheck == FALSE)
    #invalidReads <- invalidReads %>% mutate(Clones = 1)
  
      
    Pre.C.E.Count <- nrow(swappedReads)
    
    cat("Prior to counting clones, there were", Pre.C.E.Count, "reads.\n")
    
    validReadsDualFlanked <- validReads %>%
      dplyr::filter(Flanking == "Dually-Flanked") 
    validReads5Flanked <- validReads %>%
      dplyr::filter((Flanking == "5'-Flanked" & flippedOrientation == FALSE) |
                      (Flanking == "3'-Flanked" & flippedOrientation == TRUE))
    validReads3Flanked <- validReads %>%
      dplyr::filter((Flanking == "3'-Flanked" & flippedOrientation == FALSE) |
                      (Flanking == "5'-Flanked" & flippedOrientation == TRUE))
    
    #This is the section that refers to the getClones function and actually counts clones in each category. 
  
    vBF <- getClones(validReadsDualFlanked, "Dually-Flanked") #THESE ARE BEING PROBLEMATIC
    validIntermediateCount <- vBF 
    v5F <- getClones(validReads5Flanked, "5'-Flanked")
    validIntermediateCount <- rbind(validIntermediateCount, v5F, make.row.names = TRUE)
    v3F <- getClones(validReads3Flanked, "3'-Flanked")
    validIntermediateCount <- rbind(validIntermediateCount, v3F, make.row.names = TRUE)
    validIntermediateCount <- standardizeCloneIDs(validIntermediateCount)
    
    
    cat("Containing only the valid reads, clones have been quantified for", nrow(validIntermediateCount), "reads.\n")
    cat("This compares to the", nrow(validReads), "reads of all validReads.\n\n")
    
    invalidReadsDualFlanked <- invalidReads %>%
      dplyr::filter(Flanking == "Dually-Flanked") 
    invalidReads5Flanked <- invalidReads %>%
      dplyr::filter((Flanking == "5'-Flanked" & flippedOrientation == FALSE) |
                      (Flanking == "3'-Flanked" & flippedOrientation == TRUE))
    invalidReads3Flanked <- invalidReads %>%
      dplyr::filter((Flanking == "3'-Flanked" & flippedOrientation == FALSE) |
                      (Flanking == "5'-Flanked" & flippedOrientation == TRUE))
    invalidReadsNonFlanked <- invalidReads %>%
      dplyr::filter(Flanking == "Non-Flanked") 
    
    iBF <- getClones(invalidReadsDualFlanked, "Dually-Flanked")
    invalidIntermediateCount <- iBF
    i5F <- getClones(invalidReads5Flanked, "5'-Flanked")
    invalidIntermediateCount <- rbind(invalidIntermediateCount, i5F, make.row.names = TRUE)
    i3F <- getClones(invalidReads3Flanked, "3'-Flanked")
    invalidIntermediateCount <- rbind(invalidIntermediateCount, i3F, make.row.names = TRUE)
    iNF <- getClones(invalidReadsNonFlanked, "Non-Flanked")
    invalidIntermediateCount <- rbind(invalidIntermediateCount, iNF, make.row.names = TRUE)
    
    cat("Containing only the invalid reads, clones have been quantified for", nrow(invalidIntermediateCount), "reads.\n")
    cat("This compares to the", nrow(invalidReads), "reads of all validReads.\n\n")
  
    #Switches antisense coordinates back to original designation
    cat("row names somehow get added at this step where swapping is undone. \n")
    #validIntermediateCount <- validIntermediateCount %>%
    #  mutate(
    #    tempFivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeShear, NA),
    #    tempFivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeInsertion, NA),
    #    fivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", threePrimeShear, fivePrimeShear),
    #    fivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", threePrimeInsertion, fivePrimeInsertion),
    #    threePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeShear, threePrimeShear),
    #    threePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeInsertion, threePrimeInsertion)
    #  ) %>%
    #  select(-starts_with("temp"), -flippedOrientation) 
    
    #invalidIntermediateCount <- invalidIntermediateCount %>%
    #  mutate(
    #    tempFivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeShear, NA),
    #    tempFivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeInsertion, NA),
    #    fivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", threePrimeShear, fivePrimeShear),
    #    fivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", threePrimeInsertion, fivePrimeInsertion),
    #    threePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeShear, threePrimeShear),
    #    threePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeInsertion, threePrimeInsertion)
    #  ) %>%
    #  select(-starts_with("temp"), -flippedOrientation) 
    
    groupedreads <- rbind(validIntermediateCount, invalidIntermediateCount, make.row.names = TRUE)
    #groupedreads <- as.data.frame(groupedreads)
    #rownames(groupedreads) <- NULL
    
    
    
    groupedreads <- groupedreads %>%
      mutate(
        tempFivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeShear, NA),
        tempFivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeInsertion, NA),
        fivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", threePrimeShear, fivePrimeShear),
        fivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", threePrimeInsertion, fivePrimeInsertion),
        threePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeShear, threePrimeShear),
        threePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeInsertion, threePrimeInsertion)
      ) %>%
      select(-starts_with("temp"), -flippedOrientation) %>%
      as.data.frame()
    intermediate <- DataFrame(groupedreads)
    groupedreads <- as.data.frame(intermediate)
    
    # Reset the row names to ensure uniqueness
    #rownames(groupedreads) <- groupedreads$INDEX
    ##View(groupedreads)
    
    cat("In total, there are", nrow(groupedreads), "valids reads that correspond to the", Pre.C.E.Count, "reads before counting clones. \n")
    #groupedreads <- groupedreads %>% ungroup()  
    #cat("There are", nrow(groupedreads), "reads after calculating clones.\n")
    #cat("There are", totalclonecount, "total clones accounted for.\n")
    #cat("The sum total of all clones should be equal to reads. This suggests that these results are", validity, "\n")
    Cleaned.reads <- groupedreads
    Post.C.E.Count <- nrow(Cleaned.reads)
    cat("After counting, grouping, and ungrouping this dataset, there were", Post.C.E.Count, "reads. \n")
    
    #ceDiff <- Post.C.E.Count - Pre.C.E.Count
    #cat(ceDiff, "reads were lost. \n")
    
    statusCheck(Cleaned.reads)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(Cleaned.reads)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  removePCRSlips <- function(Cleaned.reads, maxPCRSlippage) {
    tic("removePCRSlips")
    ############################################################
    ## Remove those reads that are within a certain tolerance ##
    ## as PCR slippage events                                 ##
    ## NEED TO WORK ON -- DOES NOT SEEM RIGHT                 ##
    ############################################################
    cat("*************************************************\n")
    cat("_______________________Entering removePCRSlips...\n")
    cat("*************************************************\n")
    
    cat("The mechanism of integration is thought to displace some nucleotides, but generally insert at a fixed point. \nThe tolerance is", maxPCRSlippage, "and represents how far the points of integration can be.\n" )
    validReads <- subset(Cleaned.reads, passQualityCheck == TRUE) #losing NAs for slippage, non-both flanks
    invalidReads <- subset(Cleaned.reads, passQualityCheck == FALSE)
    
    cat(deparse(substitute(Cleaned.reads)), "is first being subsetted into reads based on quality check status. \nReads that have passed the quality metrics thus far are then split based on relationship to the slippage tolerance and flankedness.\n\n")
    
    #finding adequate reads
    nonslippedReads <- dplyr::filter(validReads,  slippage <= maxPCRSlippage) ### LOSING 3000 reads
    nonslippedCount <- nrow(nonslippedReads)
  
    #finding reads of intolerable slippages.
    slippedReads <- dplyr::filter(validReads,  slippage > maxPCRSlippage)
    slippedReads$passQualityCheck <- FALSE
    slippedReads$failedStep <- "removePCRSlips"
    slippedCount <- nrow(slippedReads)
    concerningreads <- slippedReads$READ
    
    #identifying the reads that would have otherwise been lost -- the single flanked reads.
    singleFlankedReads <- validReads[is.na(validReads$slippage), ]
    singleFCount <- nrow(singleFlankedReads)
    cat(nonslippedCount, "reads had slippage values within the tolerable limit.\n" )
    cat(slippedCount, "reads had intolerable slippages. These reads are included in the output, but not used for downstream analysis.\n")
    cat("\nRead IDs of concern are: \n \n" )
    for (row in concerningreads) {
      cat(row, "\n")
    }
    
    cat(singleFCount, "reads are flanked only on one side, and do not have an associated slippage value. These will be appended back into the data.\n")
    
    totalReadCount <- slippedCount + nonslippedCount + nrow(invalidReads) + singleFCount
    
    cat("There were originally", nrow(Cleaned.reads), "reads.\n")
    cat("After finding, isolating, marking, and restitching reads with intolerable slippages, there are", totalReadCount, "reads.\n")
    
    ## is.na(Cleaned.reads$slippage) this is where it duplicates -- need to extract NAs out of cleaned.reads too 
    ## so they don't get added twice -- keep the remove duplicates row below, as a safety.
    Cleaned.reads <- rbind(slippedReads, nonslippedReads, invalidReads, singleFlankedReads)
    #Cleaned.reads <- Cleaned.reads %>% distinct(fivePrimeInsertion, threePrimeInsertion, .keep_all = TRUE)
    ################ idk why, but this duplicated two rows? 
    statusCheck(Cleaned.reads)
    return(Cleaned.reads)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  createBedObject <- function(Cleaned.reads, add_chr) {
    tic("createBedObject")
    #####################################################
    ## Creates a bed file for the reads that are found ##
    ## that can then be used for adding in genome      ##
    ## annotations                                     ##
    #####################################################
    cat("*************************************************\n")
    cat("______________________Entering createBedObject...\n")
    cat("*************************************************\n")
    Clean.bed <- subset(Cleaned.reads, select = c("Flanking", "chromosome", "fivePrimeInsertion", "threePrimeInsertion", "Run_ID", "READ", "PCR.Dups", 
                                                  "Clones", "INSERT_LEN", "Clone.ID", "PCR.ID", "fivePrimeShear", "threePrimeShear", "slippage", "passQualityCheck", "failedStep", "Human.Gene.Directionality"))
    
    write.table(Clean.bed, file="Clean.bed", sep="\t")
    
    #, directionality = "Directionality" -- add in if you want directionality in the bed file
    bed <- subset(Clean.bed, select = c(chromosome = "chromosome", start_position = "fivePrimeInsertion", end_position = "threePrimeInsertion", fivePrimeShear = "fivePrimeShear", threePrimeShear = "threePrimeShear", flanking = "Flanking", slippage = "slippage", passQualityCheck = "passQualityCheck", failedStep = "failedStep", Read = "READ", 
                                        Clones = "Clones", Clone.ID = "Clone.ID", PCR.Dups = "PCR.Dups", PCR.ID = "PCR.ID", Insert.Len = "INSERT_LEN", read = "READ", directionality = "Human.Gene.Directionality")) 
    if(add_chr) {
      ## uncommon -- adds chr to the chromosome name. For the rest of the code to work, add back in the code
      ## that allows for the re-removal of chr
      bed$chromosome <- paste("chr", bed$chromosome, sep = "")  #Uncomment this code and the code below if you want a 
      write.table(x = bed, file = "bed_input.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t") 
    }
    
    bed$subject <- Clean.bed$Run_ID
    #bed$index <- row.names(bed) 
    bed <- bed <- transform(bed,
                            threePrimeInsertion = ifelse(is.na(threePrimeInsertion) | is.na(fivePrimeInsertion),
                                                         threePrimeInsertion,
                                                         pmax(threePrimeInsertion, fivePrimeInsertion)),
                            fivePrimeInsertion = ifelse(is.na(threePrimeInsertion) | is.na(fivePrimeInsertion),
                                                        fivePrimeInsertion,
                                                        pmin(threePrimeInsertion, fivePrimeInsertion)))
    bed$is.5PI.na <- is.na(bed$fivePrimeInsertion)
    bed$is.3PI.na <- is.na(bed$threePrimeInsertion)
    write.table(bed, file="bed.bed", sep="\t", row.names=FALSE, quote=FALSE)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(bed)
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  storeInvalids <- function(bed) {
    tic("storeInvalids")
    cat("*************************************************\n")
    cat("________________________Entering storeInvalids...\n")
    cat("*************************************************\n")
    
    cat("Bioconductor processing does not work on the reads that fail quality control. These reads are not being thrown out,\n but they are being removed temporarily and stored until later.\n")
    
    invalidReads <- bed %>% subset(passQualityCheck == FALSE)
    invalidCount <- nrow(invalidReads)
    
    cat("There are", invalidCount, "reads being stored as a dataframe called invalidReads.\n These will be returned at the end.\n")
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(invalidReads)
    }
  
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  createBedGaps <- function(bed) {
    tic("createBedGaps")
    cat("*************************************************\n")
    cat("______________________Entering createBedGaps...\n")
    cat("*************************************************\n")
    
    cat("Because bioconductor tools are reliant on two coordinate systems when making comparisons, the passing uni-flanked reads are temporarily being given a second coordiate one base pair away.\n")
    validCount <- nrow(subset(bed, passQualityCheck == TRUE))
    cat("Because reads that have failed quality control steps have been set aside for now, there are", validCount, "reads going through the Bioconductor portion of the pipeline.\n")
    
    cat("Subsetting duals...\n")
    dualFlanked <- bed %>% subset(passQualityCheck == TRUE & is.5PI.na == FALSE & is.3PI.na == FALSE)
    cat("Subsetting 5's...\n")
    leftFlanked <- bed %>% 
      subset(passQualityCheck == TRUE & is.5PI.na == FALSE & is.3PI.na == TRUE) 
    leftFlanked$threePrimeInsertion <- (leftFlanked$fivePrimeInsertion) #+1
    
    cat("Subsetting 3's...\n")
    rightFlanked <- bed %>% subset(passQualityCheck == TRUE & is.5PI.na == TRUE & is.3PI.na == FALSE)
    rightFlanked$fivePrimeInsertion <- (rightFlanked$threePrimeInsertion) #-1
    
    validReads <- rbind(dualFlanked, leftFlanked, rightFlanked)
    ##View(validReads) #Here for quality control purposes
    validCount <- nrow(validReads)
    cat("After processing,", validCount, "reads remain.\n")
    
    if (validCount != 0) {
      bedgaps <- makeGRangesFromDataFrame(validReads, seqnames.field = "chromosome", start.field="fivePrimeInsertion", end.field="threePrimeInsertion", keep.extra.columns = TRUE)
    } else if (validCount == 0) {
      cat("There are no reads that pass quality control at this point in time.\n")
      bedgaps <- GRanges()
      mcols(bedgaps) <- validReads
      print(bedgaps)
    }
    
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(bedgaps)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  reduceRedundantRanges <- function(bed) {
    tic("reduceRedundantRanges")
    ####################################################################
    ## Reduces redundant ranges likely associated with strand effects ##
    ## Duplicate insertions 1bp apart are likely the same -- account  ##
    ## for slippage with +-5 bp difference                            ##
    ####################################################################
    cat("Entering reduceRedundantRanges...\n")
    overlaps <- makeGRangesFromDataFrame(bed, seqnames.field = "chromosome", start.field="fivePrimeInsertion", end.field="threePrimeInsertion", keep.extra.columns = TRUE)
    overlaps <-reduce(overlaps, min.gapwidth = 5L, with.revmap = TRUE)
    overlaps <- as_tibble(overlaps)
    at <- overlaps$revmap
    at <- rapply(at, length, how = "list")
    at <- lapply(at, function(i) as.data.frame(unlist(i)))
    at <- bind_rows(at)
    overlaps$wiggle.count <- at$`unlist(i)`
    overlaps <- makeGRangesFromDataFrame(overlaps, keep.extra.columns = TRUE)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(overlaps)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  addMetadata <- function(bedgaps, overlaps) { 
    tic("addMetadata")
    #####################################################################################
    ## The very circular process of getting relevant metadata back into the new ranges ##
    ## created by reduce, also tallies clone count and pcr duplicates                  ##
    #####################################################################################
    cat("*************************************************\n")
    cat("__________________________Entering addMetadata...\n")
    cat("*************************************************\n")
    overlapPairs <- findOverlapPairs(subject = overlaps, query = bedgaps)
    overlapPairs <- as.data.frame(overlapPairs)
    duplicated_cols <- duplicated(t(overlapPairs))
    #overlapPairs <- ovelapPairs[!duplicated_cols]
    names(overlapPairs) <- gsub("first.X.", "", names(overlapPairs), fixed = TRUE)
    overlapPairs <- overlapPairs %>% group_by(seqnames, second.X.start, second.X.end, second.X.width, strand, first.INSERT_LEN, subject) %>%
      summarise(Clones = sum(Clones), PCR.Dups = sum(PCR.Dups), read = READ) 
    bedgaps <- split(overlapPairs, f = overlapPairs$subject)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(bedgaps)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  addGenicAnnotations <- function(bedgaps, annots, chr0.order) {
    tic("addGenicAnnotations")
    ## Adds annotation metadata and performs some ##
    ## quality checks                             ##
    cat("*************************************************\n")
    cat("__________________Entering addGenicAnnotations...\n")
    cat("*************************************************\n")
    
    
    cat("Adding metadata row to genic reads for downstream analysis.\n")
    qual_check <-  as.data.frame(annots, row.names=NULL)
    qual_check$region <- "genic"
    qual_check$feature <- "gene"
    qual_check$id <- row_number(qual_check)
    row.names(qual_check) <- qual_check$id
    qual_check$seqnames <- sapply(qual_check$seqnames, function(x) gsub("chr", "", x)) #SLOWish
    qual_check$seqnames <- factor(qual_check$seqnames, levels=chr0.order)
    #qual_check$id <- row.names(qual_check)
    qual_check <- qual_check[!is.na(qual_check$seqnames), ]
    g <- makeGRangesFromDataFrame(qual_check, keep.extra.columns = TRUE)
    g
    
    
    #extracting names.
    metaCols <- colnames(mcols(bedgaps))
    standardGRCols <- c("seqnames", "ranges", "strand")
    allCols <- c(standardGRCols, metaCols,"i", "id")
    
    
    cat("Mapping and counting genes.\n")
    gene.counts <- countOverlaps(query = g, subject = bedgaps, type = "any", minoverlap = 1)
    gene.counts.df <- lapply(gene.counts, function(i) as.data.frame(i)) 
    gene.counts.df <- gene.counts.df[sapply(gene.counts.df, function(x) any(x$i != 0))]
    gene.counts.df <- Map(cbind, gene.counts.df, id = names(gene.counts.df)) 
    #gene.counts.df <- lapply(gene.counts.df, `rownames<-`, NULL)
    if (length(gene.counts.df) != 0) {
      gene.counts.df <- gene.counts.df %>% bind_rows
    } else if (length(gene.counts.df) == 0) {
      gene.counts.df <- data.frame(matrix(NA, nrow = 0, ncol = length(allCols)))
      colnames(gene.counts.df) <- allCols
    }
    
    
    ##View(gene.counts.df) #here for quality control purposes
    qual_check <- merge(qual_check, gene.counts.df, by = "id") 
    qual_check <- dplyr::rename(qual_check, counts = i)
    qual_check <- subset(qual_check, counts > 0)
    
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(qual_check)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  modifyBedGaps <- function(bedgaps) {
    tic("modifyBedGaps")
    cat("*************************************************\n")
    cat("________________________Entering modifyBedGaps...\n")
    cat("*************************************************\n")
    bedgaps <- as.data.frame(bedgaps) #old version: lapply(bedgaps, function(i) as.data.frame(i))
    bedgaps <- makeGRangesFromDataFrame(bedgaps, keep.extra.columns = TRUE) #lapply(bedgaps, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(bedgaps)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  createTemp <- function(bedgaps) {
    tic("createTemp")
    cat("*************************************************\n")
    cat("___________________________Entering createTemp...\n")
    cat("*************************************************\n")
    bedgaps <- as.data.frame(bedgaps) #lapply(bedgaps, function(i) as.data.frame(i)) this is what it was, and it used to work. Now it doesn't, even though the input is basically the same (but longer?)
    #bedgaps <- lapply(bedgaps, `rownames<-`, NULL)
    #temp <- bedgaps
    #temp <- bind_rows(temp, .id = NULL)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(bedgaps)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  makeGRange <- function(dataframe) {
    ##################################################################################
    ## Function that makes Grange objects with a tolerance for when data is absent. ##
    ##################################################################################
    
    # Filter out rows with all NA values
    df_filtered <- dataframe %>%
      filter(rowSums(is.na(.)) < ncol(.))
    
    if (nrow(df_filtered) != 0) {
      grange <- makeGRangesFromDataFrame(df_filtered, keep.extra.columns = TRUE)
    } else {
      grange <- GRanges()
      mcols(grange) <- dataframe
    } 
    return(grange)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  removeDuplicatedColumns <- function(dataframe) {
    ###############################################################
    ##Removes duplicate columns created in getFullGeneInformation##
    ###############################################################
    
    cat("if you get duplication errors or missing columns, check this. \n")
    #"first.X.is.5PI.na", "first.X.is.3PI.na", add these or anything to the exclusions if stuff is missing. 
    # Columns to exclude from removal
    #added 13 Nov 2024 EDIT
    exclusions <- c("first.X.end", "first.X.Clones", "first.X.Clone.ID", "first.X.PCR.Dups" ,"second.X.strand", "second.X.counts", "first.X.failedStep")
    
    
    cat("checking problematic columns.\n")
    if (length(unique(dataframe$first.X.slippage)) == 1 | all(is.na(dataframe$first.X.slippage))) {
      slippageMissing <- TRUE
      exclusions <- append("first.X.slippage", exclusions)
    }
    if (length(unique(dataframe$first.X.is.3PI.na)) == 1 | all(is.null(dataframe$first.X.is.3PI.na))) {
      slippageMissing <- TRUE
      exclusions <- append("first.X.is.3PI.na", exclusions)
    }
    if (length(unique(dataframe$first.X.is.5PI.na)) == 1 | all(is.null(dataframe$first.X.is.5PI.na))) {
      slippageMissing <- TRUE
      exclusions <- append("first.X.is.5PI.na", exclusions)
    }
    
    cat("Exclusions:", exclusions, "\n")
   
    # Identify duplicated columns excluding the exclusions
    duplicated_cols <- dataframe[, !duplicated(t(dataframe)) | names(dataframe) %in% exclusions]
    
    # Select only the non-duplicated columns and the exclusions
    dataframe <- dataframe %>% dplyr::select(all_of(names(duplicated_cols)))
    return(dataframe)
  }  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  seqTracker <- function(timepoint) {
    ###############################################################
    ##Tracks how many valid integrations have been accounted for.##
    ###############################################################
    cat("Entering seqTracker... \n")
    validCount <- nrow(subset(Cleaned.reads, passQualityCheck == TRUE))
    cat("There are ", validCount, "integrations that have passed the quality check.\n")
  
    if (timepoint == "genic") {
      tracker <- length(unique(full[[1]]$first.READ))
      difference <- validCount - tracker
      cat("After completing integration analysis on the genic coordinates, ", tracker, "integrations have been accounted for. ", difference, " remain.\n")
      return(tracker)
    }
    
    if (timepoint == "intergenic") {
      tracker <- length(unique(full$first.READ)) 
      difference <- validCount - tracker
      cat("After completing integration analysis on the genic and intergeic coordinates, ", tracker, "integrations have been accounted for. ", difference, " remain.\n")
      return(tracker)
    }
  
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  concatNonMatched <- function(colVals) {
    
    ######################################################################################
    ## Collapses reads that map to overlapping genes and concatenates appropriate rows. ##
    ## Called by getFullGeneInformation.                                                ##
    ######################################################################################
    
    uniqueValues <- unique(colVals)
    if (length(uniqueValues) > 1) {
      return(paste(uniqueValues, collapse = "///"))
    } else {
      return(as.character(uniqueValues))
    }
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  getFullGeneInformation <- function(bedgaps, s, runID) { #remember, for testing reads$Run_ID[[1]] was used instead of runID
    tic("getFullGeneInformation")
    #######################################
    ##Gets both ranges onto the dataframe##
    #######################################
    cat("*************************************************\n")
    cat("_______________Entering getFullGeneInformation...\n")
    cat("*************************************************\n\n")
    
    names(bedgaps) <- NULL
    names(s) <- NULL
    pairs <- findOverlapPairs(query = bedgaps, subject = s, type = "any", maxgap = -1, ignore.strand = TRUE)
    genicCount <- length(unique(as.data.frame(pairs)$first.READ)) #this convoluted approach accounts for integrations in overlapped genes that are otherwise counted twice.
    intergenicCount <- integrationCount - genicCount
  
    # total <- lapply(pairs, function(i) as.data.frame(i)) #why does this only work sometimes?
    
    cat("Because overlaps have been calculated, only reads that overlap with the coordinates of a gene correspond are being tested against.\n")
    cat("Integration counts into genic features are being mapped, and the dataframe containing this being generated.\n")
  
    ###############
    ## ADDED ECR ##
    ###############
    
    total <- list(as.data.frame(pairs))
    colnames(total[[1]]) <- sub(paste(runID, ".", sep=""), "", colnames(total[[1]]))
    names(total) <- c(runID) #runID
    cat("Total:\n")
    print(total)
    
    ##View(total[[1]]) 
    #groupedreads <- groupedreads %>%
    #  mutate(
    #    tempFivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeShear, NA),
    #    tempFivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", fivePrimeInsertion, NA),
    #    fivePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", threePrimeShear, fivePrimeShear),
    #    fivePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", threePrimeInsertion, fivePrimeInsertion),
    #    threePrimeShear = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeShear, threePrimeShear),
    #    threePrimeInsertion = ifelse(Human.Gene.Directionality == "Antisense", tempFivePrimeInsertion, threePrimeInsertion)
    #  ) %>%
    #  select(-starts_with("temp"), -flippedOrientation) 
    
    full <- total
    names(full) <- gsub(paste(runID, ".", sep=""), "", names(full), fixed=TRUE)
    ##View(full[[1]]) here for quality control purposes
    cat("Removing Duplicated Columns.\n")
    ##View(full[[1]]) #included for quality control
    
    if (is.null(nrow(full[[1]])) == FALSE | (nrow(full[[1]]) == 0 && intergenicCount != 0)) {
      full <- lapply(full, removeDuplicatedColumns)
      cat("this is a test.\n")
      print(full)
      #full <- full$runID %>% select(all_of(names(duplicated_cols)))
      ## REMOVE Columns explicitly that are now duplicates
      full[[1]]$first.X.id      <- NULL
      full[[1]]$first.X.tx_id   <- NULL
      full[[1]]$first.X.gene_id <- NULL
      full[[1]]$first.X.symbol  <- NULL
      full[[1]]$first.X.type    <- NULL
      full[[1]]$first.X.region  <- NULL
      full[[1]]$first.X.feature <- NULL
      full[[1]]$first.X.counts  <- NULL
      #full[[1]]$first.X.subject <- NULL #used to ensure subject column is in the masterframe.
      full[[1]]$second.X.first.INSERT_LEN <- NULL
      full[[1]]$second.X.subject          <- NULL
      full[[1]]$second.X.Clones           <- NULL
      full[[1]]$second.X.PCR.Dups         <- NULL
      full[[1]]$second.X.read             <- NULL
      full[[1]]$second.X.feature          <- NULL #maybe?
    
      names(full[[1]]) <- gsub(".X", "", names(full[[1]]), fixed=TRUE)
      fullname <- paste(names(full[1]), ".", sep = "")
      fullname <- gsub(" ", ".", fullname)
      ##View(full) #here for debugging
      
      fullcount <- length(unique(full[[1]]$first.READ))
    } else if (is.null(nrow(full[[1]]))) {
      full <- data.frame((matrix(data = NA, ncol = 34)))
      colnames(full) <- finalNames[c(1:24,26:35)]
      full <- list(runID = full)
    }

    if ((nrow(full[[1]])) == 0 && integrationCount == 0) {
      full <- data.frame((matrix(data = NA, ncol = 34)))
      colnames(full) <- finalNames[c(1:24,26:35)]
      full <- list(runID = full)
      cat("No reads pass quality control.\n")
      cat("Printing masterframe of invalid reads.\n")
      write.xlsx(as.data.frame(invalidReads), file="masterframe.xlsx", sheetName="Low Confidence Reads", append=TRUE, row.names=FALSE, showNA = TRUE)
      cat("No reads pass quality control.\nPrinting masterframe of invalid reads.\nTerminating code.\n") 
      return(full) ####
      }

    

    cat("Debug Information:\n")
    cat("Fullname:", fullname, "\n\n")
    cat("Original Column Names:\n") 
    #cat("Modified Column Names:\n")
    for (col_name in names(full[[1]])) {
      cat(col_name, "\n")
    }
    cat("\n")
    
    names(full[[1]]) <- gsub(pattern = fullname, replacement = "", names(full[[1]]), fixed=F)
    cat("Upon entering the analysis, there were", integrationCount, "reads corresponding to the total number of integrations.\n")
    cat("Upon mapping and counting,", genicCount, "these reads correspond to integrations into genes. \nThis is the genic count.\n")
    
    cat("Modified Column Names:\n")
    for (col_name in names(full[[1]])) {
      cat(col_name, "\n")
    }
    
    full$full.symbol[is.na(full$full.symbol)] <- full$full.tx_id[is.na(full$full.symbol)]
    
    
    uncollapsedCount <- nrow(full[[1]])
    uniqueCount <- length(unique(na.omit(full[[1]]$first.READ)))
    difference <- uncollapsedCount - uniqueCount
    cat("There are", uniqueCount, "correspond to unique reads.\n", difference, "are reads integrated into overlapping features.\n")
    cat("If there are integrations where there are overlapping genes, they will be separated by '///' in the relevant columns.\n")
    
    if (nrow(full[[1]]) != 0) {
      full[[1]] <- full[[1]] %>%
        group_by(first.READ) %>%
        summarise(across(everything(), ~concatNonMatched(.))) 
      
      #flip <- list()
      #columns_to_include <- colnames(full[[1]]) #c("second.start", "second.end", "second.width", "second.id", "second.symbol", "second.type")
      #flip[[1]] <- full[[1]] %>%
      #  group_by(first.READ) %>%
      #  summarise(across(all_of(columns_to_include), ~concatNonMatched(.)))
      
      collapsedCount <- nrow(full[[1]])
      uniqueCount <- length(unique(full[[1]]$first.READ))
      difference <- uncollapsedCount - uniqueCount
      
      cat("After collapsing, there are ", collapsedCount, "rows in the dataframe.\n")
      
      full[[1]] <- as.data.frame(full[[1]]) %>% relocate(first.READ, .after = first.failedStep)  #####
      
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(full)
    } else if (nrow(full[[1]]) == 0) {
      
      # so if no reads proceed, this prevents failure of the pipeline
      cat("No reads proceeded. Moving to the next step.\n")
      full[[1]] <- as.data.frame(full[[1]])   #####
      
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(full)
    }
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  getIntergenicInformation <- function(bedgaps, s, temp) {
    tic("getIntergenicInformation")
    
    ##################################################################
    ## returns the reads that didn't map, assumed to be intergenic. ##
    ##finding all that remains by finding the difference.           ##
    ##################################################################
    cat("*************************************************\n")
    cat("_____________Entering getIntergenicInformation...\n")
    cat("*************************************************\n")
    
    tracker <- length(unique(full[[1]]$first.READ)) 
    cat("There are", tracker, "integrations that have passed the quality check.\n")
    
    
    if (integrationCount > tracker) {
      diff <- integrationCount - tracker
      cat(diff, "read(s), of the", integrationCount, " that this pipeline found originally, remain that exist in non-genic contexts.\n")
      cat("Finding all reads that did not map to genic annotations and assessing wheher they map to intergenic annotations.\n")
      cat("Extracting non-mapped reads now.\n")
      #intergenic <- setdiff.Vector(bedgaps, s) #finds all the vectors in bedgaps that do not occur in s. 
      #intergenic <- as.data.frame(intergenic)
      
      highQual <- temp %>% subset(passQualityCheck == TRUE) #WAS USING BED?
      vectorIntergenic <- setdiff(highQual$READ, full[[1]]$first.READ)
      intergenic <- temp[temp$READ %in% vectorIntergenic,]
      #intergenic$strand <- "*"
      #intergenic <- intergenic %>% 
      #  relocate(strand, .after = threePrimeInsertion) %>%
      #  dplyr::rename(
      #  "seqnames" = "chromosome",
      #  "start" = "fivePrimeInsertion",
      #  "end" = "threePrimeInsertion"
      #)
      
      #names(intergenic[1:5]) <- c("seqnames", "start", "end", "width", "strand")
      #intergenic <- left_join(intergenic, temp)
      #intergenic <- split(intergenic, f = intergenic$subject)
      intergenic <- makeGRangesFromDataFrame(intergenic, keep.extra.columns = TRUE, na.rm = FALSE) ###### I THINK I HAVE TO MIRROR over the missing coordinates
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      print(intergenic)
      return(intergenic)
    } else {
      cat("All reads are accounted for, no remaining reads.\n")
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(NULL)
    }
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  addAdditionalAnnotations <- function(refGenome, RDS_DIR) {
    tic("addAdditionalAnnotations")
    
    #importing other, non-gene annotations and counting them
    ## Will need to update for Rhesus genome
    cat("*************************************************\n")
    cat("______________Entering addAdditionalAnnotations...\n")
    cat("*************************************************\n")
  
    #cat("ensure intergenic read column is named appropriately.\n")
    tracker <- length(unique(full[[1]]$first.READ)) + length(unique(as.data.frame(intergenic)$READ.1)) ##STN
    cat(tracker, "integrations of the total are accounted for. In total, there are", integrationCount,"integrations.\n")
    
    if(refGenome == "Hs") {
    	if(USE_ENSEMBL == TRUE) {
         cat("Building other annotations now.\n")
         cat("Using Ensembl.\n")
         annots <- as.list(builtin_annotations())
         annot.list <- c("hg38_genes_intergenic", "hg38_cpg_islands") #, 
                         #  "hg38_cpg_shores", "hg38_cpg_shelves", "hg38_cpg_inter", "hg38_enhancers_fantom", 
                         # "hg38_cpgs")
         complete.list <- c("hg38_genes_intergenic", "hg38_cpg_islands",  "hg38_cpg_shores", "hg38_cpg_shelves", "hg38_cpg_inter", "hg38_enhancers_fantom", "hg38_cpgs")
         cat("Currently only pulling out\n", annot.list, "\nbut", complete.list, "\nis possible.\n")
      
         other.annots <- build_annotations(genome = "hg38", annotations = annot.list)
        }
        if(USE_ENSEMBL == FALSE) {
          cat("Building other annotations now.\n")
          cat("Using table browser.\n")
          other.annots <- read.delim(paste(coordinateDir, "intergenic.bed", sep=""), header = FALSE)
          other.annots$region <- "intergenic"
          other.annots$feature <- NA
          
          other.annots <- other.annots %>% dplyr::rename(
            "seqnames" = "V1",
            "start" = "V2",
            "end" = "V3"
          )
        }
      }
      if(refGenome == "RheMac") {
        other.annots <- readRDS(paste(RDS_DIR, "RheMacOtherRanges.rds", sep="\\"))     
      } 
    
	elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(other.annots)
    
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  updateOtherAnnotations <- function(other.annots, chr0.order) {
    tic("updateOtherAnnotations")
    
    cat("*************************************************\n")
    cat("_______________Entering updateOtherAnnotations...\n")
    cat("*************************************************\n")
    
    tracker <- length(unique(full[[1]]$first.READ)) + length(unique(as.data.frame(intergenic)$READ.1)) ####
    #if (integrationCount > tracker) {
      cat("Formating these non-genic annotations.\n")
      o.annots.check <- as.data.frame(other.annots)
      o.annots.check$seqnames <- sapply(o.annots.check$seqnames, function(x) gsub("chr", "", x))
      o.annots.check$seqnames <- factor(o.annots.check$seqnames, levels=chr0.order)
      o.annots.check <- o.annots.check[!is.na(o.annots.check$seqnames), ]
      row.names(o.annots.check) <- o.annots.check$id
    #}
    #if (integrationCount == tracker) {
    #  cat("Nothing needs to be done.\n")
    #  o.annots.check <- NULL
    #}
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    #print(o.annots.check)
    return(o.annots.check)
  }
  
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  makeOACaGR <- function(o.annots.check) {
    ##########################################
    ##converting o.annots.check toa GRANGE. ##
    ##########################################
    tic("makeOACaGR")
    cat("*************************************************\n")
    cat("___________________________Entering makeOACaGR...\n")
    cat("*************************************************\n\n")
  
    tracker <- length(unique(full[[1]]$first.READ)) + 
      length(unique(as.data.frame(intergenic)$READ.1)) ####
    
    
    
    #if (integrationCount > tracker) {
      cat("Reformating o.annots.check into a GenomicRanges object.\n")
      #print(o.annots.check)
      converted <- makeGRangesFromDataFrame(o.annots.check, keep.extra.columns = TRUE)
      print(converted)
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(converted)
    #}
    #if (integrationCount == tracker) {
    #  cat("Nothing needs to be done.\n")
    #  converted <- NULL
    #  elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    #  elapsed_time
    #  return(converted)
    #}
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  clean_column_names <- function(dataframe) {
    ###################################################################
    ##Cleaning column name artifacts generated by repeated joinings. ##
    ###################################################################
    
    modifiedTemplate <- gsub("-", ".", templateToRemove)
    
    new_names <- gsub(paste0(".X", modifiedTemplate, ".csv."), "", names(dataframe), fixed = TRUE)
    new_names <- gsub(paste0(modifiedTemplate, ".csv."), "", names(dataframe), fixed = TRUE)
    new_names <- gsub(paste0("X.", templateToRemove, ".csv."), "", names(dataframe), fixed = TRUE)
    new_names <- gsub(paste0("X", templateToRemove, ".csv."), "", names(dataframe), fixed = TRUE)
    new_names <- gsub(paste0("X", templateToRemove, ".csv."), "", names(dataframe), fixed = TRUE)
    new_names <- gsub(paste0(templateToRemove, ".csv."), "", new_names, fixed = TRUE)
    new_names <- gsub("X.", "", new_names, fixed = TRUE)
    new_names <- gsub("X", "", new_names, fixed = TRUE)
    
    
    
    names(dataframe) <- new_names
    return(dataframe)
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
   updateIntergenicCounts <- function(other.annots, intergenic) {
    ###########################
    ##Counting and appending!##
    ###########################
    
    tic("updateIntergenicCounts")
    cat("*************************************************\n")
    cat("_______________Entering updateIntergenicCounts...\n")
    cat("*************************************************\n\n")
    
    cat("Counting intergenic annotations.\n")
    cat("Counting, Appending, Cleaning Intergenic annotations.\n")
    tracker <- length(unique(full[[1]]$first.READ)) #+ length(unique(intergenic)) #+ length(unique(intergenic$P2.csv.READ))
    cat("Genic reads accounted for:", tracker, "of", integrationCount, "\n")
    
    #cat("Determining what full and full[[1]] are. \n")
    #print(full)
    #print(full[[1]])
    
    correctNames <- c("first.seqnames", "first.start", "first.end", "first.width", "first.strand", 
                      "first.fivePrimeShear", "first.threePrimeShear", "first.Flanking", "first.slippage",
                      "first.passQualityCheck", "first.failedStep", "first.READ" , "first.Clone.ID",
                      "first.PCR.ID", "first.INSERT_LEN",  "first.is.5PI.na", "first.subject",
                      "first.is.3PI.na", "first.Human.Gene.Directionality", "second.start", "second.end", "second.width", "second.strand", "second.counts", "second.id",
                      "second.symbol", "second.region", "second.type")
    
    if (integrationCount > tracker) {
      cat("Integration Count is greater than the total number of mapped integrations.\n")
      intergenic.counts <- countOverlaps(query = other.annots, subject = intergenic, type = "any")
      intergenic.counts.df <- as.data.frame(intergenic.counts)
      intergenic.counts.df$id <- row.names(intergenic.counts.df)
      intergenic.counts.df <- intergenic.counts.df[intergenic.counts.df$intergenic.counts >= 1,]
      intergenic.counts.df <- intergenic.counts.df %>% dplyr::rename(i = intergenic.counts)
      #print(intergenic.counts.df)
      test <- sum(intergenic.counts.df$i)
      cat(test, "total intergenic integrations.\n")
      ##View(intergenic.counts.df) #for quality control purposes.
      if (test == 0) {
        elapsed_time <- toc(log = TRUE)
        elapsed_time
        print(full) #full[[1]])
        return(full) #full[[1]])
      } else {
        cat("Clearly, intergenic integrations remain. \n")
        intergenicRanges <- as.data.frame(other.annots)
        #print(intergenicRanges)
        intergenicRanges$id <- seq_len(nrow(intergenicRanges))
        intergenic.counts.df <- intergenic.counts.df[intergenic.counts.df$i >= 1,] ######
        #intergenic.counts.df <-  Map(cbind, intergenic.counts.df, id = names(intergenic.counts.df)) 
        #intergenic.counts.df <- mapply(cbind, intergenic.counts.df, "subject"=names(intergenic.counts.df), SIMPLIFY=F)
        intergenic.counts.df <- `rownames<-`(intergenic.counts.df, NULL)
        
        print(head(intergenicRanges))
        intergenic.counts.df <- as.data.frame(merge(intergenicRanges, intergenic.counts.df, by = "id")) 
        cat("Post-merge status check.\n")
        print(intergenic.counts.df)
        intergenic.counts.df <- dplyr::rename(intergenic.counts.df, counts = i)
        intergenic.counts.df <- intergenic.counts.df %>% dplyr::relocate(id, .after = counts)
        cat("Total number of intergenic integration sites:", nrow(intergenic.counts.df),"\n")
        cat("Total number of integrations:", sum(intergenic.counts.df$counts), "\n")
        intergenic.counts <- makeGRangesFromDataFrame(intergenic.counts.df, keep.extra.columns = TRUE)
        
        intergenic.counts.df <- as.data.frame(findOverlapPairs(query = intergenic, subject = intergenic.counts, type = "within", ignore.strand = TRUE)) 
        cat("Printing intergenic.counts.df. \n")
        ##View(intergenic.counts.df) #for quality control purposes.
        columnsToHave <- c("first.X.seqnames", "first.X.start", "first.X.end", "first.X.width", "first.X.strand",
                           "first.X.fivePrimeShear", "first.X.threePrimeShear", "first.X.Flanking", "first.X.slippage",
                           "first.X.passQualityCheck", "first.X.failedStep", "first.X.READ", "first.X.Clone.ID",
                           "first.X.PCR.ID", "first.X.INSERT_LEN", "first.X.is.5PI.na", "first.subject",
                           "first.X.is.3PI.na", "first.X.Human.Gene.Directionality", "second.X.start", "second.X.end", "second.X.width",
                           "second.X.strand", "second.X.counts", "second.X.id",  "second.X.region", "second.X.symbol", "second.X.type")
        
        intergenic.counts.df <- intergenic.counts.df[columnsToHave[1:26]] #[correctNames]
        cat("Printing intergenic.counts.df. \n")
        cat("Number of columns in intergenic.counts.df:", ncol(intergenic.counts.df), "\n")
        print(colnames(intergenic.counts.df))
        
        
        
        #intergenic.counts.df <- intergenic.counts.df %>% dplyr::select(all_of(columnsToHave))
        
        #need to subset based on columnsToHave
        
        
        if (nrow(intergenic.counts.df) > 0) {
          intergenic.counts.df$second.X.symbol <- as.character(NA)
          intergenic.counts.df$second.X.region <- as.character(intergenic.counts.df$second.X.region)
          intergenic.counts.df$second.X.type <- as.character(NA)
          intergenic.counts.df$first.X.Clones <- as.numeric(NA)
          intergenic.counts.df$first.X.PCR.Dups <- as.numeric(NA)
          #colnames(intergenic.counts.df) <- gsub("X\\.", "", colnames(intergenic.counts.df))
          
          
          cat(nrow(full[[1]]), "genic integrations.\n")
          cat(nrow(intergenic.counts.df), "intergenic integrations.\n")
          cat("number of columns in intergenic.counts.df:", ncol(intergenic.counts.df), "\n")
          cat("number of columns in full[[1]]:", ncol(full[[1]]), "\n")
          cat("Renaming.\n")
          colnames(full[[1]]) <- sub("\\.", ".X.", colnames(full[[1]]))
          colnames(intergenic.counts.df) <- sub("first.subject", "first.X.subject", colnames(intergenic.counts.df))
          cat("full[[1]]:\n",colnames(full[[1]]), "\n")
          cat("intergenic.counts.df:\n",colnames(intergenic.counts.df), "\n")
          cat("Binding.\n") ####losing column class correspondance
          
          sharedCols <- intersect(colnames(full[[1]]), colnames(intergenic.counts.df))
          print(sharedCols)
          cat(setdiff(colnames(full[[1]]), colnames(intergenic.counts.df)))
          
          intergenic.counts.df <- intergenic.counts.df[, colnames(full[[1]])]
          #intergenic.counts.df <- full[[1]][, colnames(intergenic.counts.df)]
          
          #columnClasses <- sapply(intergenic.counts.df, class)
          #class(intergenic.counts.df)
          #class(full[[1]])
          
          #test <- full[[1]]
          
          
          for (col in 1:30) {
            if (class(full[[1]][[col]]) != class(intergenic.counts.df[[col]])) {
              # Convert intergenic.counts.df column to match the class of full[[1]]
              target_class <- class(intergenic.counts.df[[col]])
              print(target_class)
              if (target_class == "character") {
                full[[1]][[col]] <- as.character(full[[1]][[col]])
              } else if (target_class == "numeric") {
                full[[1]][[col]] <- as.numeric(full[[1]][[col]])
              } else if (target_class == "integer") {
                full[[1]][[col]] <- as.integer(full[[1]][[col]])
              } else if (target_class == "factor") {
                full[[1]][[col]] <- as.factor(full[[1]][[col]])
              } else if (target_class == "logical") {
                full[[1]][[col]] <- as.logical(full[[1]][[col]])
              } else {
                stop("Unsupported data class:", target_class)
              }
            }
          }
          
          #intergenic.counts.df <- full[[1]][, colnames(intergenic.counts.df)]
          
          full <- dplyr::bind_rows(full[[1]], intergenic.counts.df) ###STNEDIT
          
          print(full)
        }
          
          
          #for (col in columnsToHave) {
      		#    if (!col %in% colnames(intergenic.counts.df)) {
          #		    cat("adding columns.\n")
          #		    intergenic.counts.df[[col]] <- if (nrow(intergenic.counts.df) == 0) character(0) else NA
      		#    }
  	    	#}
  
          
          #colnames(intergenic.counts.df) <- correctNames                                               
          print(intergenic.counts.df)
          
          if (nrow(intergenic.counts.df) == 0) {
      	  	ongoing <- nrow(full)
            difference <- integrationCount - ongoing
            cat(difference, "reads are not yet accounted for.\n")
            ##View(full) #for quality control purposes
          
            elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
            elapsed_time
            return(full)
          }
  
          #View(intergenic.counts.df) #for quality control purposes.
          
          #names(intergenic.counts.df) <- gsub("X.P2.csv.", "", names(intergenic.counts.df), fixed = TRUE)
          #names(intergenic.counts.df) <- gsub("P2.csv.", "", names(intergenic.counts.df), fixed = TRUE)
          #names(intergenic.counts.df) <- gsub("X.", "", names(intergenic.counts.df), fixed = TRUE)
          
          #intergenic.counts.df <- clean_column_names(intergenic.counts.df) #figure out why this doesn't always work?
          
          cat("Correct column names:\n")
          print(correctNames)
          cat("Number of correct column names:", length(correctNames), "\n")
    
            
            #colnames(intergenic.counts.df) <- correctNames
            #if (!is.data.frame(full)) {
    		  	#cat("Expected a data frame in full[[1]], but got something else.")
    		  	#return(full) #this was the edit
		    }
		full <- clean_column_names(full)

    if ((ncol(intergenic.counts.df) - 2) != length(correctNames)) { #3 = clones, pcr.dups, are not represented
  		cat("Mismatch in number of columns: colnames cannot be assigned.\n")
      cat("Length of intergenic.counts.df:", ncol(intergenic.counts.df), "\n")
      cat("Length of correctNames:", ncol(correctNames), "\n")
  		return()
		}
        ##View(intergenic.counts.df) #for quality control purposes
        ##View(full) #for quality control purposes
        #full <- full[[1]] adjusted in 8e5 run -- may uncomment if clean_comment_names on full breaks things down the line.
        ##View(full)
        ##View(intergenic.counts.df)
    cat("full column names\n")
    str(full)     
    cat("intergenic.counts.df column names\n")
		str(intergenic.counts.df)
		
		#for (col_name in colnames(full)) {
  	#		if (col_name %in% colnames(intergenic.counts.df)) {
    #			# Get data type of the column in `full`
    #			target_type <- class(full[[col_name]])
    #
    #			# Coerce `intergenic.counts.df` column to the same type
    #			intergenic.counts.df[[col_name]] <- as(intergenic.counts.df[[col_name]], target_type)
  	#		}
		#}
		
		#intergenic.counts.df <- intergenic.counts.df[, colnames(full)]
		
		#if (!all(colnames(full) == colnames(intergenic.counts.df))) {
  	#		cat("Column names between full and intergenic.counts.df do not match.")
  	#		return()
		#}

    #    full <- plyr::rbind.fill(full, intergenic.counts.df) ####
    #    if (is.null(full) || nrow(full) == 0) {
  	#		cat("Resulting full data frame is NULL or empty after rbind.fill.")
  	#		return()
		#}
        
      full <- full %>% relocate(first.READ, .after = first.strand)
      ongoing <- nrow(full)
      difference <- integrationCount - ongoing
      cat(difference, "reads are not yet accounted for.\n")
      ##View(full) #for quality control purposes
        
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(full)
    
    }
      #return(full)
     
    if (integrationCount == tracker) {
      cat("Nothing needs to be done.\n")
      full <- full[[1]]
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(full)
    }
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  updateOtherAnnotationsWithIntergenicCounts <- function(o.annots.check, full) {
    tic("updateOtherAnnotationsWithIntergenicCounts")
    
    ## Updates the other annotations with the intergenic data
    cat("******************************************************\n")
    cat("Entering updateOtherAnnotationsWithIntergenicCounts...\n")
    cat("******************************************************\n")
    
    tracker <- length(unique(full[[1]]$first.READ)) + length(unique(as.data.frame(intergenic)$READ.1)) ####
    
    if (integrationCount == tracker){
      cat("Nothing to add.\n")
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(o.annots.check)
    } else if (nrow(full) > nrow(tracker)) {
      o.annots.check <- merge(o.annots.check, full, by = "id") #I think that it does work, it's just a pain to verify
      o.annots.check$region <- "intergenic"
      o.annots.check <- dplyr::rename(o.annots.check, feature_biotype = type)
      o.annots.check <- dplyr::rename(o.annots.check, feature_name = id)
      o.annots.check <- dplyr::rename(o.annots.check, counts = i)
      o.annots.check <- o.annots.check %>% dplyr::relocate(feature_name, .after = tx_id)
      o.annots.check <- subset(o.annots.check, counts > 0)
      cat(nrow(o.annots.check), "reads map to intergenic annotations.\n") #verify if this is sensical.
      row.names(o.annots.check) <- NULL
      
    } else
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(o.annots.check)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  updateIntergenicAsGRange <- function(intergenic) {
    tic("updateIntergenicAsGRange")
    ## Updates the intergenic annotations as a GRange object
    cat("*************************************************\n")
    cat("_____________Entering updateIntergenicAsGRange...\n")
    cat("*************************************************\n")
    tracker <- length(unique(full$first.READ))+ length(unique(intergenic))
    if (integrationCount == tracker){
      cat("Nothing to add.\n")
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return(NULL)
    } else if (integrationCount > tracker) {
    intergenic <- as.data.frame(intergenic)
    intergenic <- lapply(intergenic, `rownames<-`, NULL)
    intergenic <- makeGRangesFromDataFrame(intergenic, keep.extra.columns = TRUE)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(intergenic)
    }
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  createNonGeneInformation <- function(o.annots.check) {
    tic("createNonGeneInformation")
    ##################################################################### 
    #append o.annots.check and qual.check together, removing excess rows#
    #####################################################################
    
    
    cat("*************************************************\n")
    cat("_____________Entering createNonGeneInformation...\n")
    cat("*************************************************\n\n")
    
    tracker <- length(unique(full$first.READ))
    cat("There were", integrationCount, "integrations.\n")
    cat("The tracker count is", tracker, "right now. \n")
    
    if (integrationCount > tracker) {
      cat("Finding the reads that don't map to the already described genic or intergenic Ensembl annotations but instead map to nongenes. \nLikely empty, but want to ensure all data is captured.\n")
      t <- makeGRangesFromDataFrame(o.annots.check, keep.extra.columns = TRUE)
      ##View(as.data.frame(t)) for quality control purposes.
      pairs <- findOverlapPairs(query = t, subject = intergenic, type = "within", maxgap = -1, ignore.strand = TRUE)
      ##View(as.data.frame(pairs)) for quality control purposes.
      non.genes <- lapply(pairs, function(i) as.data.frame(i))
      if (nrow(as.data.frame(pairs)) > 0) {
        non.genes <- base::Reduce(full_join, non.genes) #above -- presents duplicated as an array? is that why it's being wonky?
        non.genes <- non.genes[!duplicated(as.list(non.genes))]
        names(non.genes) <- gsub(".X", "", names(non.genes), fixed = TRUE)
        cat(nrow(non.genes), "reads map to non-gene annotations.\n")
      } else {
        cat("Nothing has slipped through the cracks. No reads fit this category.\n")
        non.genes <- (as.data.frame(non.genes))
        print(non.genes)
        cat("Nothing to add!\n")
      }
    }
    if (integrationCount == tracker){
      cat("Nothing to add.\n")
      non.genes <- NULL
    }
    ## Don't actually use any of the three lines below ##
    fullnames <- base::colnames(full)
    nongenenames <- base::colnames(non.genes)
    common.cols <- intersect(fullnames, nongenenames)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    print(non.genes)
    return(non.genes)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  combineDatasets <- function(full, non.genes, outliers, excludeGeneList, runID) {
    tic("combineDatasets")
    ########################################################################
    ## Creates a master frame from the full gene and non gene information ##
    ########################################################################
    cat("*************************************************\n")
    cat("______________________Entering combineDatasets...\n")
    cat("*************************************************\n")
    
    cat("printing full.\n", nrow(full), "\n")
    cat("Number of columns in full:", ncol(full), "\n")
    print(head(full))
    cat("printing non.genes.\n", nrow(non.genes), "\n")
    cat("Number of columns in non.genes:", ncol(non.genes), "\n")
    print(head(non.genes))
    cat("printing outliers.\n", nrow(outliers), "\n")
    cat("Number of columns in outliers:", ncol(outliers), "\n")
    print(head(outliers))

    
    cat("binding dataframes.\n")
    
    
    master.frame <- rbind(full, non.genes, outliers)
    colnames(master.frame) <- sub("^[^.]+\\.csv\\.", "", colnames(master.frame))
    print(colnames(master.frame))
    print(nrow(master.frame))
    #head(master.frame)
    #print(master.frame)
    #print(colnames(master.frame))
    #print(excludeGeneList)
    
    master.frame <- master.frame[!grepl(excludeGeneList, master.frame$second.symbol),]
    print(master.frame)
     # Check if 'second.symbol' column exists
    #if ("second.symbol" %in% colnames(master.frame)) {
    #    master.frame <- master.frame[!grepl(excludeGeneList, master.frame$second.symbol),]
    #} else {
    #    cat("Warning: 'second.symbol' column not found in master.frame. Skipping filtering step.\n")
    #}
    ##View(master.frame) #uncomment for quality control. 
    
    
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(master.frame)
    #Error: vector memory exhausted (limit reached?) fix with R_MAX_VSIZE=100Gb but (maybe make it)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  writeMasterFrame <- function(master.frame) {
    tic("writeMasterFrame")
    ## Writes out the master frame to a csv file
    
    cat("Creating deliverable xlsx.\n")
    cat("Creating masterframe.xlsx.\n Sheet 1 will contain reads that passed quality control and corresponded to integrations that are sensical.\nSheet 2 are reads that failed quality control for any reason.\n")
    #openxlsx::write.xlsx(as.data.frame(master.frame), file="masterframe.xlsx", sheetName="High Confidence Reads", row.names=FALSE, showNA = TRUE)
    #openxlsx::write.xlsx(as.data.frame(invalidReads), file="masterframe.xlsx", sheetName="Low Confidence Reads", append=TRUE, row.names=FALSE, showNA = TRUE)
    master.frame <- subset(master.frame, select = -c(cloneCount, pcrDuplicateCount, countOfReadsMapped))
    invalidReads <- subset(invalidReads, select = -c(Clones, PCR.Dups)) 
    
    
    wb <- openxlsx::createWorkbook()
    
    openxlsx::addWorksheet(wb, sheetName = "High Confidence Reads")
    openxlsx::writeData(wb, sheet = "High Confidence Reads", as.data.frame(master.frame), rowNames = FALSE)
    
    # Add the second sheet with "Low Confidence Reads"
    openxlsx::addWorksheet(wb, sheetName = "Low Confidence Reads")
    openxlsx::writeData(wb, sheet = "Low Confidence Reads", as.data.frame(invalidReads), rowNames = FALSE)
    
    # Save the workbook
    openxlsx::saveWorkbook(wb, file = "masterframe.xlsx", overwrite = TRUE)
    
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  getFigColors <- function(master.frame) {
    tic("writeMasterFrame")
    ## gets the figure colors based on the number of gene symbols
    
    color.count <- length(unique(master.frame$first.symbol))
    AnnotColors <- randomColor(count = color.count)
    FigColors   <- tibble(symbol = unique(master.frame$first.symbol), color= AnnotColors)
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    return(FigColors)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  findOutliers <- function(intergenic, other.annots) {
    tic("findOutliers")
    ## finds outliers that are in intergenic but not in other annotations
    cat("*************************************************\n")
    cat("_________________________Entering findOutliers...\n")
    cat("*************************************************\n\n")
    cat("Finding anything else that might remain and didn't map to anything on the reference database.\n")
    #finding outliers
    cat("intergenic is", class(intergenic), "\n")
    cat("other.annots is", class(other.annots), "\n")
    cat("if one of these is NULL, assess why setdiff.Vector is now failing.\n")
    
    if (is.null(other.annots)) {
    	cat("other.annots is NULL. Returning intergenic as outliers.\n")
      elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
      elapsed_time
      return()
	} else if (is.null(intergenic)) {
    	cat("intergenic is NULL. Returning intergenic as outliers.\n")
  	  elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
  	  elapsed_time
	    return()
	} else {
		cat("Both other.annots and intergenic are non-empty. Discrepencies can be identified. \n")
		outliers <- setdiff.Vector(intergenic, other.annots) #returns everything that is in intergenic that isn't in other.annots
		outliers <- as.data.frame(outliers)
		if (nrow(outliers) > 0) {
		  outliers <- bind_rows(Map(cbind, outliers, name=names(outliers))) 
		  cat(nrow(outliers), "read(s) is/are outliers worth exploring individually.\n")
		  print(outliers)
		  elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
		  elapsed_time
		  return(outliers)
		} else {
		  cat("No reads are outliers. \n")
		  elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
		  elapsed_time
		  return(outliers)
		}
	}
    #cat("Printing outliers.\n")
    
    #outliers <- as.data.frame(outliers)
    #if (nrow(outliers) > 0) {
    #  outliers <- bind_rows(Map(cbind, outliers, name=names(outliers))) 
    #  cat(nrow(outliers), "read(s) is/are outliers worth exploring individually.\n")
    #} else {
    #  cat("No reads are outliers. \n")
    #}
    #elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    #elapsed_time
    #print(outliers)
    #return(outliers)
  }
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  getOthers <- function(combined) {
    tic("getOthers")
    cat("*************************************************\n")
    cat("____________________________Entering getOthers...\n")
    cat("*************************************************\n")
    ###############################################
    ##Mapping CpGs, Repeats, ENCODE annotations. ##
    ###############################################
    
    combinedLength <- nrow(as.data.frame(combined))
    cat("Number of rows in combinedLength:", combinedLength, "\n")
    
    if(USE_ENSEMBL == TRUE) {
      master.frame <- combined
    }
    
    if(USE_ENSEMBL == FALSE) {
      cat("Building other annotations now.\n")
      cat("Using Table Browser data.\n")
      other.annots <- read.delim(paste(coordinateDir, "externalFixed.bed", sep=""), header = FALSE)
      
      other.annots$V1 <- gsub("chr", "", as.character(other.annots$V1))
      other.annots$strand <- "*"
      other.annots$width <- abs(other.annots$V2 - other.annots$V3)
      
      col_order <- c("V1", "V2", "V3", "width", "strand", "V4", "V5")
      other.annots <- other.annots[, col_order]
      other.annots <- other.annots %>% dplyr::rename(
        "seqnames" = "V1" ,
        "start" = "V2",
        "end" = "V3",
        "category" = "V4",
        "subtype" = "V5")
      other.annots <- other.annots %>% mutate(
        category = ifelse(grepl("^CpG: \\d+$", category), "CpG Island", category),
        category = ifelse(grepl("LTR?", category), "LTR", category),
        category = ifelse(grepl("DNA?", category), "DNA", category),
        category = ifelse(grepl("RC?", category), "RC", category),
        category = ifelse(grepl("SINE?", category), "SINE", category),
      )
      rownames(other.annots) <- NULL
      #cat("other.annots with column name changes:\n")
      #print(other.annots)
      other.annotsGR <- makeGRangesFromDataFrame(unique(other.annots), keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)#, 
                                                 #seqnames.field = "V1", start.field = "V2", end.field =  "V3", strand.field = "V5")
      cat("Printing other.annotsGR:\n")
      print(head(other.annotsGR))
      names(combined) <- gsub("first.", "", names(combined), fixed = TRUE)
      #names(combined) <- gsub("second.", "", names(combined), fixed = TRUE)
      if (nrow(combined) > 0) {
        cat("Reads identified.\n")
        print(colnames(combined))
        combined <- combined %>% relocate(READ, .after = Flanking) 
        rownames(combined) <- NULL}
      if (nrow(combined) == 0) {
        cat("No reads identified.\n")
        rownames(combined) <- NULL
      }
      
      cat("Taking processed data and re-making it as a GenomicRanges object.\n")
      allGRanges <- makeGRange(combined)
      cat("This will be re-converted as a dataframe for use later, but it's being used to find overlaps now.\n")
      
      pairs <- findOverlapPairs(query = allGRanges, subject = other.annotsGR, type = "any") 
      cat("printing pairs for ENCODE elements:\n")
      print(pairs)
      ##View(as.data.frame(pairs)) for quality control purposes.
      if (length(pairs) == 0) {
        pairs <- as.data.frame(allGRanges)
        pairs$second.category <- NA
      } #else if (length(pairs) > 0) {
      	#pairs <- as.data.frame(allGRanges)
        #pairs$second.category <- NA
      #} 
      
      master.frame <- as.data.frame(pairs, row.names = NULL) 
      cat("print(master.frame)\n")
      print(master.frame)
      master.frame <- master.frame %>%
        mutate(
          isRepetitive = case_when(
            second.category == "Simple_repeat" ~ "Simple Repeat",
            second.category == "Satellite" ~ "Satellite Repeat",
            second.category == "LINE" ~ "LINE",
            second.category == "DNA" ~ "DNA Repeat Element",
            second.category == "SINE" ~ "SINE",
            second.category == "RC" ~ "Rolling Circle Repeat Element",
            second.category == "Low_complexity" ~ "Low Complexity Repeat Element",
            second.category == "Unknown" ~ "Unknown",
            TRUE ~ NA_character_
          ),
          isENCODE = case_when(
            second.category == "enhP" ~ "Proximal Enhancer-Like Signature",
            second.category == "enhD" ~ "Distal Enhancer-Like Signature",
            second.category == "prom" ~ "Promoter-like Signature",
            second.category == "CTFC" ~ "CTFC",
            second.category == "K4m3" ~ "DNase-H3K4me3",
            TRUE ~ NA_character_
          ),
          isCpG = case_when(
            second.category == "CpG Island" ~ "CpG Island",
            TRUE ~ NA_character_
          )
        )
      master.frame <- master.frame[-c(30:(ncol(master.frame)-3))] #historically 30:62
      master.frame <- as.data.frame(master.frame) ###20241216
      cat(class(master.frame))
      #print(master.frame)
      print(master.frame)
      names(master.frame) <- gsub("second.X.", "second.", names(master.frame))
      colcheck <- ncol(master.frame)
      names(master.frame) <- gsub("first.X.", "", names(master.frame))
      
      #classList <- lapply(master.frame[1:(ncol(master.frame)-3)], class)
      #for (i in 2:ncol(combined)) {
      #  class(combined[[i]]) <- classList[i]
      #}
      cat("Appending metadata.\nThis is being converted from the GenomicRanges because of the concatNonMatched function's propensity to modify class structure.\n")
      combined <- as.data.frame(allGRanges)
      
      if (nrow(master.frame) !=0) {
        master.frame <- full_join(combined,master.frame)}## is this adding more?
      #master.frame <- master.frame[-c(30:(ncol(master.frame)-3))] #this wasn't historically in here.
      if ((nrow(master.frame) == 0) && nrow(combined) == 0) {
        elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
        elapsed_time
        cat("The count of unique reads at the start was:", length(unique(reads$READ)), "\n")
        cat("The number of high confidence reads is:", length(unique(master.frame$readID)), "\n")
        cat("The number of low confidence reads is:", length(unique(invalidReads$READ)), "\n")
        cat("The number of unique reads accounted for now is:", 
        (length(unique(master.frame$readID)) + length(unique(invalidReads$READ))), "\n") #deleted this at the end let's hope this fucking works)
        cat("The difference is:", (nrow(reads)-(length(unique(master.frame$readID)) + length(unique(invalidReads$READ)))), "\n")
        return(master.frame)
      }
      
      
      length(unique(full$first.READ))+ length(unique(intergenic))
      
      master.frame <- master.frame %>%
        group_by_at(vars(-isRepetitive, -isENCODE, -isCpG)) %>%  # Group by all columns except those just added
        summarize(
          isRepetitive = ifelse(all(is.na(isRepetitive)), NA, dplyr::first(na.omit(isRepetitive))),    # Retain the first non-NA value isRepetitive
          isENCODE = ifelse(all(is.na(isENCODE)), NA, dplyr::first(na.omit(isENCODE))),            # Retain the first non-NA value isENCODE
          isCpG = ifelse(all(is.na(isCpG)), NA, dplyr::first(na.omit(isCpG))),                  # Retain the first non-NA value isCpG  
        ) %>%
        ungroup()
    }
    
    cat("Renaming columns with more human-friendly names. \n")
    
    if (nrow(master.frame) != 0) {
      master.frame <- master.frame %>% dplyr::rename(
        chromosome = seqnames,
        fivePrimeInsertion = start,
        threePrimeInsertion = end,
        width = width, #bioconductor counting -- adds one. Removing downstream because slippage is already present.
        strand = strand,
        fivePrimeShear = fivePrimeShear,
        threePrimeShear = threePrimeShear,
        flanking = Flanking,
        slippage = slippage,
        passQualityCheck = passQualityCheck, #removing because passes and fails will be stored in separate tabs.
        failedStep = failedStep, #removing because passes and fails are stored separately.
        readID = READ,
        cloneCount = Clones,
        cloneID = Clone.ID,
        pcrDuplicateCount = PCR.Dups,
        pcrID = PCR.ID,
        insertLength = INSERT_LEN,
        orientation = Human.Gene.Directionality,
        is5PIna = is.5PI.na,
        is3PIna = is.3PI.na,
        annotationFivePrimeCoordinate = second.start,
        annotationThreePrimeCoordinate = second.end,
        annotationWidth = second.width,
        annotationStrand = second.strand,
        annotationID = second.id,
        annotationType = second.type,
        annotationRegion = second.region,
        annotationFeature = second.type,
        symbol = second.symbol,
        countOfReadsMapped = second.counts #total number of reads that map to that annotation. 
      )
      master.frame$Clones <- NULL
      master.frame$PCR.Dups <- NULL
    }  
    if (nrow(master.frame) == 0) {
      cat("There is no input. Creating the correctly formatted, empty dataframe. \nThis is just to ensure consistency of formating as a deliverable.\n\n")
      master.frame$chromosome <- NULL
      master.frame$subject <- NULL
      master.frame$first.chromosome <- NULL
      master.frame$fivePrimeInsertion <- NULL
      master.frame$threePrimeInsertion <- NULL
      master.frame$second.strand  <- NA
      master.frame$second.id <- NA
      master.frame$second.type <- NA
      master.frame$second.region <- NA
      master.frame$second.symbol <- NA
      master.frame$second.counts <- NA
      
      master.frame <- master.frame %>% dplyr::rename(
        chromosome = seqnames,
        fivePrimeInsertion = start,
        threePrimeInsertion = end,
        width = width, #bioconductor counting -- adds one. Removing downstream because slippage is already present.
        strand = strand,
        fivePrimeShear = fivePrimeShear,
        threePrimeShear = threePrimeShear,
        flanking = Flanking,
        slippage = slippage,
        passQualityCheck = passQualityCheck, #removing because passes and fails will be stored in separate tabs.
        failedStep = failedStep, #removing because passes and fails are stored separately.
        readID = READ,
        cloneCount = Clones,
        cloneID = Clone.ID,
        pcrDuplicateCount = PCR.Dups,
        pcrID = PCR.ID,
        insertLength = INSERT_LEN,
        orientation = Human.Gene.Directionality,
        is5PIna = is.5PI.na,
        is3PIna = is.3PI.na,
        annotationFivePrimeCoordinate = second.start,
        annotationThreePrimeCoordinate = second.end,
        annotationWidth = second.width,
        annotationStrand = second.strand,
        annotationID = second.id,
        annotationType = second.type,
        annotationRegion = second.region,
        annotationFeature = second.type,
        symbol = second.symbol,
        countOfReadsMapped = second.counts #total number of reads that map to that annotation. 
      )
    }  
    
    if (USE_ENSEMBL == TRUE) {
      master.frame <- master.frame %>% dplyr::rename(
        annotationStableID = second.tx_id,  ### only in Ensembl?
        geneID = second.gene_id, ### remove in GenomeBrowser Version?
      )
    }
    
    cat("Columns renamed.\n")
    cat("Removing uni-flanked *assumed* columns generated to allow Bioconductor analysis to proceed.\n")
    
    master.frame$threePrimeInsertion[master.frame$flanking == "5'-Flanked"] <- NA
    master.frame$fivePrimeInsertion[master.frame$flanking == "3'-Flanked"] <- NA
    master.frame <- subset(master.frame, select = -c(annotationID))
    
    columns_vector <- c("fivePrimeShear", "threePrimeShear", "cloneCount", "pcrDuplicateCount", 
                        "insertLength", "annotationFivePrimeCoordinate", "annotationThreePrimeCoordinate", 
                        "annotationWidth", "countOfReadsMapped")
    
    for (col in columns_vector) {
      master.frame[col] <- as.numeric(unlist(master.frame[col]))
    }
    
    
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
    
    #uniqueReads <- (nrow(master.frame) + length(unique(invalidReads$READ)))
    cat("The count of unique reads at the start was", nrow(unique(reads)))
    cat("\nThe number of unique reads accounted and mapped for now is", 
    (nrow(master.frame) + nrow(invalidReads)))
    cat("\nThis breaks down to", nrow(master.frame), "and", length(unique(invalidReads$READ)), "for high confidence and low confidence reads, respectively.\n")
    
    return(master.frame)


  }
 
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  getExternalGeneAnnotations <- function(refGenome, RDS_DIR) {
    tic("getExternalGeneAnnotations")
    if(USE_ENSEMBL == TRUE) {
      ## Gets the external gene annotations from annovar or from the ones that ##
      ## we created                                                            ##
      cat("*************************************************\n")
      cat("_________________In getExternalGeneAnnotations...\n")
      cat("*************************************************\n")
      if(refGenome == "Hs") {
        gene.annot.list <- c("hg38_basicgenes")
        cat("Using Ensembl.\n")
        annots <- build_annotations(genome = "hg38", annotations = gene.annot.list)
      }
      if(refGenome == "RheMac") {
      annots <- readRDS(paste(RDS_DIR, "RheMacGeneRanges.rds", sep="\\"))     
      }
      return(annots)
    }
    if(USE_ENSEMBL == FALSE){
      #USES GENERATED ANNOTATIONS
      if(refGenome == "Hs") {
        cat("Using table browser. \n")
        annots <- read.delim(paste(coordinateDir, "genicFixed.bed", sep=""), header = FALSE)
        annots$V1 <- gsub("chr", "", as.character(annots$V1))
        annots$strand <- "*"
        annots$width <- abs(annots$V2 - annots$V3)
        
        col_order <- c("V1", "V2", "V3", "width", "strand", "V4", "V5")
        annots <- annots[, col_order]
        annots <- annots %>% dplyr::rename(
          "seqnames" = "V1",
          "start" = "V2",
          "end" = "V3",
          "symbol" = "V4",
          "type" = "V5"
        )
        if (USE_ENSEMBL == TRUE) {
        annots <- makeGRangesFromDataFrame(annots, keep.extra.columns = TRUE)
        }
        if (USE_ENSEMBL == FALSE) {
          annots <- makeGRangesFromDataFrame(annots, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
        }
        return(annots)
      }
    }
    elapsed_time <- toc(log = TRUE)  # Stop timer and log timing
    elapsed_time
  }
  
  
  #--------------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------------
  
  qcTimePlot <- function(log.lst) {
    #########################################################################
    ## Makes a 'survival' curve of how long the function ran at each step. ##
    #########################################################################
    
    cat("*************************************************\n")
    cat("___________________________Entering qcTimePlot...\n")
    cat("*************************************************\n\n")
    cat("Making QC plot.\n")
    
    times <- base::as.data.frame(unlist(log.lst))
    times <- separate(times, col = 'unlist(log.lst)', into = c("prefix", "number"), sep = ":")
    times$number <- trimws(times$number, whitespace = " sec elapsed")
    times$number <- as.numeric(times$number)
    times$cumulativeNumber <- cumsum(times$number)
    times$prefix <- factor(times$prefix, levels = unique(times$prefix))
    times$propTime <- round(times$number/sum(times$number), digits = 2)
    
    
    
    timePlot <- ggplot(times, aes(x = prefix, y = cumulativeNumber)) +
      geom_step(aes(group = 1)) + geom_point() +
      labs(x = "Process Step", y = "Cumulative Time (sec)", subtitle = "time per step (% of total time)") +
      ggtitle("Cumulative Time Taken for Each Process Step") +
      geom_text(aes(label = paste(number, " (", propTime, ")", sep = ""), y = max(cumulativeNumber)),  angle = 90, hjust = 1) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = -0.01), plot.margin = unit(c(1,1,1,1), "cm"))
      
    pdf("qcTimePlot.pdf", width = 10, height = 10, pointsize = 12, bg = "transparent") 
    print(timePlot)
    dev.off()
    
  }
  
  #--------------------------------------------------------------------------------------------------------
  
  
  #========================================================================================================
  #########################################################################################################
  ##                                              MAIN PROGRAM                                           ##
  #########################################################################################################
  #========================================================================================================



cat("Starting script!\n\n")
#InstallBioconductorPackages(biocPackages)
#InstallRPackages(RPackages)
cat("loading libraries...\n\n")
#loadLibraries(biocPackages, RPackages)
#loadLibraries(biocPackages, biocPackageVersions)
#loadLibraries(RPackages, RPackageVersions)

tic.clearlog()
#tic("Script start after loading packages.\n")
size.Reference <- getReferenceSizeFile(REF_SIZE_DIR, refGenome)
insert.refSize <- getReferenceInsertSize(insertType)  ## Set reference insertion size
chr0.order     <- getChromosomeOrder(refGenome)       ## Order of Chromosome names
reads          <- getReadsFromCSVFiles()

Cleaned.reads <- cleanReads(reads, tolerance, exciseFlag) 
Cleaned.reads <- findMappingErrors(Cleaned.reads)
Cleaned.reads <- findGapDifferences(Cleaned.reads, insert.refSize)
Cleaned.reads <- markPCRDuplicates(Cleaned.reads)
Cleaned.reads <- findOddFails(Cleaned.reads)

writePhyloMasterFrame(Cleaned.reads)

Cleaned.reads <- removeBadReads(Cleaned.reads, hostTolerance, minReadLength) 
Cleaned.reads <- markClonalExpansion(Cleaned.reads)  
integrationCount <- dplyr::filter(Cleaned.reads, passQualityCheck == TRUE) %>% nrow()

write.table(Cleaned.reads, file="Cleaned.reads.PCRDuplicatesandClonalExpansionMarked.txt", sep="\t", row.names=FALSE)
write.table(Cleaned.reads, file="Cleaned.reads.PCRDuplicatesClonalExpansionPCRSlipsMarked.txt", sep="\t", row.names=FALSE)

bed <- createBedObject(Cleaned.reads, ADD_CHR_FLAG)
invalidReads <- storeInvalids(bed)
bedgaps <- createBedGaps(bed)
annots <- getExternalGeneAnnotations(refGenome, RDS_DIR)

qual_check <- addGenicAnnotations(bedgaps, annots, chr0.order)
temp       <- createTemp(bedgaps)
s <- makeGRange(qual_check)
write.table(s, file="s.bed", sep="\t", row.names=FALSE, quote=FALSE)
write.table(bedgaps, file="bedgaps.bed", sep="\t", row.names=FALSE, quote=FALSE)

## Get all of the annotation information ##
full                 <- getFullGeneInformation(bedgaps, s, reads$Run_ID[[1]])
tracker              <- length(unique(full[[1]]$first.READ))
intergenic           <- getIntergenicInformation(bedgaps, s, temp)
other.annots         <- addAdditionalAnnotations(refGenome, RDS_DIR)
o.annots.check       <- updateOtherAnnotations(other.annots, chr0.order)
other.annots         <- makeOACaGR(o.annots.check) 
full                 <- updateIntergenicCounts(other.annots, intergenic) 
#o.annots.check      <- updateOtherAnnotationsWithIntergenicCounts(o.annots.check, full)
#intergenic          <- updateIntergenicAsGRange(intergenic) 
non.genes            <- createNonGeneInformation(o.annots.check)
outliers             <- findOutliers(intergenic, other.annots)
combined             <- combineDatasets(full, non.genes, outliers, excludeGeneList) #how to combine cleaned.reads and master.frame
master.frame         <- getOthers(combined)   

writeMasterFrame(master.frame)
#toc("script done!")
log.lst <- tic.log(format = TRUE)
write.table(x= log.lst, file = "executionTiming.txt", col.name = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
qcTimePlot(log.lst)

sessionInfo()

cat("The count of unique reads at the start was", nrow(unique(reads)))
cat("\nThe number of unique reads accounted and mapped for now is", (nrow(master.frame) + nrow(invalidReads)))
cat("\nThis breaks down to", nrow(master.frame), "and", length(unique(invalidReads$READ)), "for high confidence and low confidence reads, respectively.\n")
    

#========================================================================================================
#########################################################################################################
##                                               END PROGRAM                                           ##
#########################################################################################################
#========================================================================================================
