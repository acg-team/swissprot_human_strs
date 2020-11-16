#!/usr/bin/env R
source_all <- function(directory) {
  for(file in list.files(directory)) {source(paste0(directory, file))}
}
