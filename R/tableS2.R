# Phylogenetic generalised least squares
###

# Prepare data 

#1: Select only columns with phylogenetic information
#2: Remove rows with any NAs
#3: Convert all to factor
taxa <- df %>% select(superkingdom,phylum,class,order,family,genus,species) %>%
  drop_na %>% 
  mutate_all(factor)

#Load modified as.phylo.formula function
#This script overrides the original as.phylo.formula function
#The new version adds branch lengths to each node which is 
#required to run phylolm below.
source("R/as.phylo.formula.R")

#Create phylo tree using formula
frm <- ~superkingdom/phylum/class/order/family/genus/species

#Create phylogenetic tree
print("Creating phylo tree. This may take a while...")
tree <- as.phylo.formula(frm, data = taxa)
print("Done")

#Use phylolm to calculate phylogenetic generalised least squares for each variable
#This has to be done on each variable independently, as the function leaves out any 
#record with any missing data


#Start output table
tb <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tb) <- c("model","r2","tips","nodes")
tb$model <- as.character()
tb$r2 <- as.numeric()
tb$tips <- as.integer()
tb$nodes <- as.integer()

####
# temp_adjusted_maxgrowth vs log10(genome_size)
####

#Describe model (for table)
model <- "Tmp. adj. max growth rate vs. log10(genome size)"

#Define traits
trait1 <- "temp_adjusted_maxgrowth"
# VS
trait2 <- "genome_size"

#1: Get trait columns with species
#2: Remove any NAs
#3: Rename columns to fit code below
tmp <- df %>% select("species",.data[[trait1]],.data[[trait2]]) %>% 
  drop_na %>%
  rename("trait1" = .data[[trait1]], "trait2" = .data[[trait2]])

#Reduce main tree to contain only species with this trait data (listed in tmp)
sub_tree <- drop.tip(phy = tree, tree$tip.label[!(tree$tip.label %in% tmp$species)])

#Reduce trait table to only contain organisms that are present in this sub tree
tmp <- tmp %>% filter(tmp$species %in% sub_tree$tip.label)

#Add species name to rownames
rownames(tmp) <- tmp$species

#Run phylogenetic generalised least squares 
#NOTE: CHECK what needs to be log10 transformed and adjust
phylm <- phylolm(tmp$trait1~log10(tmp$trait2), data = tmp, phy = sub_tree)

#Add data to table
dat <- data.frame(model = model, r2 = as.numeric(phylm$r.squared), tips = as.integer(phylm$n), nodes = as.integer(sub_tree$Nnode))
tb <- tb %>% bind_rows(dat)

####
# log10(d1_mid) vs log10(genome_size)
####

#Describe model (for table)
model <- "log10(d1_mid) vs. log10(genome_size)"

#Define traits
trait1 <- "d1_mid"
# VS
trait2 <- "genome_size"

#1: Get trait columns with species
#2: Remove any NAs
#3: Rename columns to fit code below
tmp <- df %>% select("species",.data[[trait1]],.data[[trait2]]) %>% 
  drop_na %>%
  rename("trait1" = .data[[trait1]], "trait2" = .data[[trait2]])

#Reduce main tree to contain only species with this trait data (listed in tmp)
sub_tree <- drop.tip(phy = tree, tree$tip.label[!(tree$tip.label %in% tmp$species)])

#Reduce trait table to only contain organisms that are present in this sub tree
tmp <- tmp %>% filter(tmp$species %in% sub_tree$tip.label)

#Add species name to rownames
rownames(tmp) <- tmp$species

#Run phylogenetic generalised least squares 
#NOTE: CHECK what needs to be log10 transformed and adjust
phylm <- phylolm(log10(tmp$trait1)~log10(tmp$trait2), data = tmp, phy = sub_tree)

#Add data to table
dat <- data.frame(model = model, r2 = as.numeric(phylm$r.squared), tips = as.integer(phylm$n), nodes = as.integer(sub_tree$Nnode))
tb <- tb %>% bind_rows(dat)

####
# log10(d1_mid) vs temp_adjusted_maxgrowth
####

#Describe model (for table)
model <- "log10(d1_mid) vs. temp_adjusted_maxgrowth"

#Define traits
trait1 <- "d1_mid"
# VS
trait2 <- "temp_adjusted_maxgrowth"

#1: Get trait columns with species
#2: Remove any NAs
#3: Rename columns to fit code below
tmp <- df %>% select("species",.data[[trait1]],.data[[trait2]]) %>% 
  drop_na %>%
  rename("trait1" = .data[[trait1]], "trait2" = .data[[trait2]])

#Reduce main tree to contain only species with this trait data (listed in tmp)
sub_tree <- drop.tip(phy = tree, tree$tip.label[!(tree$tip.label %in% tmp$species)])

#Reduce trait table to only contain organisms that are present in this sub tree
tmp <- tmp %>% filter(tmp$species %in% sub_tree$tip.label)

#Add species name to rownames
rownames(tmp) <- tmp$species

#Run phylogenetic generalised least squares 
#NOTE: CHECK what needs to be log10 transformed and adjust
phylm <- phylolm(log10(tmp$trait1)~tmp$trait2, data = tmp, phy = sub_tree)

#Add data to table
dat <- data.frame(model = model, r2 = as.numeric(phylm$r.squared), tips = as.integer(phylm$n), nodes = as.integer(sub_tree$Nnode))
tb <- tb %>% bind_rows(dat)

rm(dat)

#Adjust column names for output
names(tb) <- c("Model","PGLS R_squared","Number of tips (species)","Number of nodes")

#Save table S4
write.table(tb, file = "output/stats/tableS2.csv", sep = ",", quote = FALSE, row.names = F)
