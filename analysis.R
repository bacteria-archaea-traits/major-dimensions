#################
# Load packages #
#################

library("ggpubr")
library("png")
library("tidyverse")
library("Hmisc")
library("agricolae") #For New Tukeys method
library("devtools") #For pgls stats (Table S2)
library("phytools") #For pgls stats (Table S2)
library("phylolm") #For pgls stats (Table S2)
library("caTools") #For regression analysis (Table S2)

options(stringsAsFactors = FALSE)

###########
# Process #
###########

# Retrieve large "taxonomy_names.csv" file not in the GitHub repo
if (!file.exists("data/madin_et_al/taxonomy_names.csv")) {
  download.file(url="https://ndownloader.figshare.com/files/14875220?private_link=ab40d2a35266d729698c", destfile = "data/madin_et_al/taxonomy_names.csv")
}

# Prepare data frames
source("R/prep.R")

# Create main figures
source("R/figures.R")

# Create supplementary figures
source("R/figures_supplementary.R")

# Generate general stats output
source("R/stats.R")

# Generate table outputs
source("R/table1.R")
source("R/table2.R")
source("R/tableS2.R")
source("R/tableS3.R")
source("R/tableS4.R")
