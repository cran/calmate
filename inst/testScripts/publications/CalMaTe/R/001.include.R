###########################################################################
# 001.include.R
#
# Author: Henrik Bengtsson
# 
# Description:
# This script defines functions and other objects that are needed by
# all scripts in this directory.
###########################################################################
library("aroma.cn");
library("calmate");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

library("R.archive");
setArchiveOption("devEval", TRUE);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
calTags <- c("<none>", "CMTN,v1", "CMTN,v1,refs=N", "CMTN,v2", "CMTN,v2,refs=N")[-2];
anTags <- c("dens", "hets", "2Mb");

hetCol <- c("#999999", "gray", "red", "blue", "orange")[1];



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
toPNG <- function(name, width=840, ...) {
  devEval("png", name=name, ..., width=width)$fullname;
} # toPNG()


loadSets <- function(..., verbose=FALSE) {
  dsT <- AromaUnitTotalCnBinarySet$byName(..., verbose=verbose);
  dsB <- AromaUnitFracBCnBinarySet$byName(..., verbose=verbose);
  list(total=dsT, fracB=dsB);
} # loadSets()



extractSignals <- function(dsList, sampleName, ..., reference=c("auto", "none", "median"), verbose=FALSE) {
  # Argument 'reference':
  reference <- match.arg(reference);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting total CN signals");
  
  verbose && enter(verbose, "Extracting sample of interest");
  verbose && cat(verbose, "Sample name: ", sampleName);
  
  idxT <- indexOf(dsList$total, sampleName);
  dfT <- getFile(dsList$total, idxT);
  idxB <- indexOf(dsList$fracB, sampleName);
  dfB <- getFile(dsList$fracB, idxB);
  
  verbose && print(verbose, list(tumor=dfT, fracB=dfB));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Extracting TCNs");
  tcn <- extractRawCopyNumbers(dfT, logBase=NULL, ...);
  verbose && print(verbose, tcn);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting BAFs");
  baf <- extractRawAlleleBFractions(dfB, ...);
  verbose && print(verbose, baf);
  verbose && exit(verbose);

  if (reference == "auto") {
    verbose && enter(verbose, "Inferring from data set name if CN ratios needs to be calculated");
    hasRatios <- hasTag(dsList$total, "CMTN");
    verbose && cat(verbose, "Has 'CMTN' tag: ", hasRatios);
    reference <- ifelse(hasRatios, "none", "median");
    verbose && cat(verbose, "Inferred argument 'reference': ", reference);
    verbose && exit(verbose);
  }
  
  if (reference == "median") {
    dfTR <- getAverageFile(dsList$total, verbose=less(verbose,5));
    
    verbose && enter(verbose, "Extracting TCNs for reference pool");
    tcnR <- extractRawCopyNumbers(dfTR, logBase=NULL, ...);
    verbose && print(verbose, tcnR);
    verbose && exit(verbose);
    
    tcn <- divideBy(tcn, tcnR);
    tcn$y <- 2*tcn$y;
  }
  
  res <- list(tcn=tcn, baf=baf);
  verbose && exit(verbose);
  
  res;
} # extractSignals()


extractCACB <- function(dsList, units, ..., reference=c("auto", "none", "median"), samplesLast=FALSE) {
  # Argument 'reference':
  reference <- match.arg(reference);

  C <- extractMatrix(dsList$total, units=units, ...);
  B <- extractMatrix(dsList$fracB, units=units, ...);
  stopifnot(identical(dim(B), dim(C)));

  if (reference == "auto") {
    verbose && enter(verbose, "Inferring from data set name if CN ratios needs to be calculated");
    hasRatios <- hasTag(dsList$total, "CMTN");
    verbose && cat(verbose, "Has 'CMTN' tag: ", hasRatios);
    reference <- ifelse(hasRatios, "none", "median");
    verbose && cat(verbose, "Inferred argument 'reference': ", reference);
    verbose && exit(verbose);
  }

  if (reference == "median") {
    muC <- matrixStats::rowMedians(C, na.rm=TRUE);
    C <- 2 * C / muC;
  }

  CB <- B*C;
  CA <- C - CB;
  data <- array(c(CA, CB), dim=c(dim(C), 2));
  dimnames(data)[[3]] <- c("A", "B");

  if (samplesLast) {
    data <- aperm(data, perm=c(1,3,2));
  }

  data;
} # extractCACB()


##########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Added option "auto" for argument 'reference' of extractSignals().
# o Created.
###########################################################################
