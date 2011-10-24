# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all predefined data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathT <- file.path(path, "dataSets");
filenames <- list.files(path=pathT, pattern="[.]R$");
dataSetList <- lapply(filenames, FUN=function(filename) {
  env <- new.env();
  sourceTo(filename, path=pathT, envir=env);
  env;
});
names(dataSetList) <- gsub("[.]R$", "", filenames);
dataSetList <- lapply(dataSetList, FUN=function(env) as.list(env));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Choose data set to study
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSets <- names(dataSetList);
if (interactive() && require("R.menu")) {
  dataSet <- textMenu(dataSets, title="Select data set:", value=TRUE);
} else {
  dataSet <- dataSet[1];
}

attachLocally(dataSetList[[dataSet]]);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Choose chip type to study, if more than one
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (length(chipTypes) > 1 && interactive() && require("R.menu")) {
  chipType <- textMenu(chipTypes, title="Select chip type:", value=TRUE);
} else {
  chipType <- chipTypes[1];
}

printf("Data set: %s\n", dataSet);
printf("Chip type: %s\n", chipType);
