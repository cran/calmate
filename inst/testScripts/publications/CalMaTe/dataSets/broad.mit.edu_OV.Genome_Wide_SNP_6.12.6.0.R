dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.12.6.0";
tags <- "ASCRMAv2";
chipTypes <- "GenomeWideSNP_6";


sampleName <- "TCGA-23-1027";

segments <- list(
  "normal (1,1)" = c(2, 108.0, 123.5),
  "gain (1,2)"   = c(2, 125.0, 140.5), 
  "CN-LOH (0,2)" = c(2, 141.5, 157.0),
  "mix"          = c(2, 108.0, 157.0)
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
