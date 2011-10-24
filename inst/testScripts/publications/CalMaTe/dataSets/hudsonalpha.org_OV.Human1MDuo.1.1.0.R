dataSet <- "hudsonalpha.org_OV.Human1MDuo.1.1.0";
tags <- "XY";
chipTypes <- "Human1M-Duo";

sampleName <- "TCGA-23-1027";

segments <- list(
  "normal (1,1)" = c(2, 105.0, 124.0),
  "gain (1,2)"   = c(2, 125.0, 140.0), 
  "CN-LOH (0,2)" = c(2, 150.0, 170.0)
);

getPairs <- function(dsList, ...) {
    fullnames <- getFullNames(dsList$total);
    sampleNames <- gsub(",(total|fracB).*", "", fullnames);
    sampleNames <- gsub("-[0-9]{2}[A-Z]-[0-9]{4}-[0-9]{2}.*", "", sampleNames);

    names <- gsub("-[0-9]{2}[A-Z]", "", sampleNames);
    uNames <- unique(names);
    pairs <- matrix(NA, nrow=length(uNames), ncol=2);
    rownames(pairs) <- uNames;
    colnames(pairs) <- c("tumor", "normal");
    for (name in uNames) {
      patternT <- sprintf("%s-(01)[A-Z]$", name);
      patternN <- sprintf("%s-(10|11)[A-Z]$", name);
      idxT <- grep(patternT, sampleNames);
      idxN <- grep(patternN, sampleNames);
      pair <- c(idxT[1], idxN[1]);
      pair <- sampleNames[pair];
      pairs[name,] <- pair;
    }
    pairs;
} # getPairs()
