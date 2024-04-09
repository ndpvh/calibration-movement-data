# First, get all libraries read in
library(tidyverse)
library(data.table)
library(ggpubr)
library(modeest)

# Get all filenames within the "utility" directory, except, of course, 
# utility.R itself
filenames <- list.files(file.path("utility"),
                        pattern = ".R",
                        full.names = TRUE,
                        recursive = TRUE)
filenames <- filenames[filenames != "utility/utility.R"]

# Loop over these filenames and load in all functions that they define
lapply(filenames, source) %>% 
    invisible()
