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
