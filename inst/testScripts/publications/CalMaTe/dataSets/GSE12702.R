# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# GSE12702
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataSet <- "GSE12702";
tags <- "ACC,-XY,BPN,-XY,RMA,FLN,-XY";
chipTypes <- "Mapping250K_Nsp";


sampleName <- "GSM318736";
  
segments <- list(
  "normal (1,1)" = c(8,  0.0, 16.0),
  "loss (0,1)"   = c(8, 16.5, 45.0), 
  "gain (1,2)"   = c(8, 45.0,150.0)
);

getPairs <- function(dsList, ...) {
  # For the GSE12702 data set, it happens to be that samples, when sorted
  # lexicographically, are ordered in pairs of tumors and normals.
  sampleNames <- getNames(dsList$total);
  sampleNames <- matrix(sampleNames, nrow=20, ncol=2, byrow=TRUE);
  colnames(sampleNames) <- c("tumor", "normal");

  pairs <- array(sampleNames, dim=dim(sampleNames));
  colnames(pairs) <- c("tumor", "normal");
  rownames(pairs) <- sampleNames[,"tumor"];

  ## When could add patient ID annotation too.
  ## rownames(pairs) <- c(24, 25, 27, 31, 45, 52, 58, 60, 75, 110, 115, 122,
  ##                      128, 137, 138, 140, 154, 167, 80, 96);

  pairs;
} # getPairs()
