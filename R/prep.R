# prep.R

###################
# Load data files #
###################

# Taxonomy files (see Madin et al. 2020)
nam <- read.csv("data/Madin et al.2020/taxonomy_names.csv", as.is=TRUE)
nam <- nam[, c("tax_id","name_txt")] #Removes columns that are not needed

tax <- read.csv("data/Madin et al.2020/ncbi_taxmap.csv", as.is=TRUE)
tax <- unique(tax[, names(tax)]) #Removes duplicate entries

# Main data frame (see Madin et al. 2020)
df <- read.csv("data/Madin et al.2020/condensed_species_GTDB[NCBI_fill]_16102020.csv", as.is=TRUE)
df <- df[!is.na(df$species),]

# Environment table (see Madin et al. 2020)
environments <- read.csv("data/Madin et al.2020/environments.csv", as.is=TRUE)

# Intracellular organismsm and mycoplasma (see Madin et al. 2020)
intr <- read.csv("data/Madin et al.2020/intracellular_organisms.csv")
myco <- read.csv("data/Madin et al.2020/mycoplasma.csv")

# Strategy table
strategies <- read.csv("data/microorganism_strategies.csv", as.is = TRUE)
strategies[strategies == ""] <- NA


#####################
# Layout parameters #
#####################

save_path <- "output/figures"

# Define fixed formats for all figures
text_size <- 9
text_color <- "black"
font_family <- "sans"
plot_line_color <- "black"
vector_colour <- "black"
plot_line_width <- 0.5 #mm

# Define theme standards 
basic_layout <- 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = plot_line_color, size = plot_line_width),
    panel.grid = element_blank(),
    axis.title = element_text(size = text_size, family=font_family, color=text_color),
    axis.text = element_text(size = text_size, family=font_family, color=text_color),
    axis.line = element_blank(),
    axis.ticks = element_line(color = plot_line_color, size = plot_line_width)
  )


#Define colour sets

colours_raw <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
#Colour of intracellular organisms
intra_col <- "yellow2"
#Rearrange colours
colours <- c(colours_raw[1],colours_raw[2],colours_raw[4],colours_raw[6],colours_raw[5],colours_raw[8],"#BDBDBD")
colours_with_intra <- c(colours_raw[1],colours_raw[2],colours_raw[4],colours_raw[6],colours_raw[5],colours_raw[8],intra_col,"#BDBDBD")
two_factor <- c("#FF8C19","#2471a3")
presence_absence <- c("gray","red","green")

# Define special axes
growth_rate_axis_text <- expression(Max~growth~rate~(hr^{-1}))
growth_rate_norm_axis_text <- "Temp-adj max growth rate (log_10 resid)"
volume_axis_text <- expression(paste("Cell volume (",mu,m^3,")", sep = ""))


#######################
# Fix any data errors #
#######################

#Fix habitat type for "Halonatronum saccharophilum". Listed as 'soil' but should be 'sediment_hypersaline'
df$isolation_source[df$species == "Halonatronum saccharophilum"] <- "sediment_hypersaline"


######################
# Create data groups #
######################

# Order data in output by setting levels manually
df$metabolism <- factor(df$metabolism, levels = c("obligate aerobic", "aerobic",  "facultative", "microaerophilic", "anaerobic","obligate anaerobic"))

#Create shapeagg column
df$shapeagg <- df$cell_shape
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("coccus")] <- "spheroid"
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("bacillus","coccobacillus","vibrio")] <- "rod"
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("pleomorphic","filament", "star", "spiral", "irregular", "flask", "spindle", "fusiform", "disc", "disc ", "square", "branced", "triangular")] <- NA


######################################
# Calculate mid diameter and volumes #
######################################

df$d1_mid <- ifelse(!is.na(df$d1_up), (df$d1_lo + df$d1_up)/2, df$d1_lo)
df$d2_mid <- ifelse(!is.na(df$d2_up), (df$d2_lo + df$d2_up)/2, df$d2_lo)

tmp <- df[!is.na(df$d1_mid) & !is.na(df$shapeagg) & df$shapeagg %in% c("rod","spheroid"), c("species_tax_id","d1_mid","d2_mid","shapeagg","cell_shape")]
tmp$volume <- NA

for(i in 1:nrow(tmp)) {
  if(!is.na(tmp$d1_mid[i])) {
    if(tmp$shapeagg[i] == "spheroid") {
      
      d <- NA
      
      #If there is a descrepancy bewteen width and lenght, use mean
      if(!is.na(tmp$d2_mid[i])) {
        #Calculate mean of the two values
        d <- (tmp$d1_mid[i]+tmp$d2_mid[i])/2
      } else {
        d <- tmp$d1_mid[i]
      }
      
      tmp$volume[i] <- 4/3*pi*(d/2)^3
      
    } else {
      if(!is.na(tmp$d2_mid[i])) {
        
        #Calculate as rod with hemispherical ends 
        
        #end volume:
        ends <- 4/3*pi*(tmp$d1_mid[i]/2)^3
        #body: length minus diameter used in ends
        body <- pi*(tmp$d1_mid[i]/2)^2*(tmp$d2_mid[i]-tmp$d1_mid[i])
        
        if(body>0) {
          tmp$volume[i] <- ends+body
        } else {
          tmp$volume[i] <- ends
        }
        
      }
    }
  }
}
tmp <- tmp[!is.na(tmp$volume),]

#Move to main data frame
df <- df %>% left_join(tmp[,c("species_tax_id","volume")], by = "species_tax_id")
rm(tmp)


###################
# Add growth rate #
###################

df$growth_rate <- log(2)/df$doubling_h


######
# Merge with Mist transporter data
######

# ms <- read.csv("../general data files/MIST/clean_mist01022019.csv")
# #Add tcp count
# df <- ms %>% group_by(species_tax_id) %>%
#   mutate(tcp_tot = sum(tcp.hk,tcp.hhk,tcp.rr,tcp.hrr,tcp.other, na.rm = TRUE)) %>%
#   mutate(tcp_hk_hkk_other_chemotaxis = sum(tcp.hk,tcp.hhk,tcp.other,tcp.chemotaxis, na.rm = TRUE)) %>%
#   mutate(tcp_hk_hkk = sum(tcp.hk,tcp.hhk, na.rm = TRUE)) %>%
#   select(species_tax_id, ocp, tcp_hk_hkk_other_chemotaxis, tcp_hk_hkk, tcp.chemotaxis, tcp_tot, majormodes_total) %>%
#   right_join(df, by = "species_tax_id")
# rm(ms)


###########################
# Add Cobo Simon habitats #
###########################

df <- df %>% left_join(environments[,c("Type","Cobo.Simon.habitat")], by = c("isolation_source"="Type"))
rm(environments)


###############################
# Create new habitat category #
###############################

#First letter upper case for environment
df$Cobo.Simon.habitat <- str_to_title(df$Cobo.Simon.habitat, locale = "en")

#Use Cobo Simon habitat scheme
df$habitat <- as.character(df$Cobo.Simon.habitat)

#Set any isolation_source that does not fit into the Cobo Simon scheme to "Other"
df$habitat[!is.na(df$isolation_source) & is.na(df$habitat)] <- "Other"

#Fix particular types of Cobo Simon habitats
df$habitat[df$habitat == "Rock"] <- "Other"
df$habitat[df$habitat == "Therm"] <- "Thermal"
df$habitat[grepl("plant|fungus|algae",df$isolation_source)] <- "Other"
df$habitat[grepl("feces|endotherm_surface|ectotherm_surface",df$isolation_source)] <- "Other"
df$habitat[df$isolation_source %in% c("host","host_animal")] <- "Other"

#Where habitat is not "Other" split host into 'endo' and 'ecto'
df$habitat[grepl("endotherm",df$isolation_source) & !(df$habitat == "Other")] <- "Endotherm"
df$habitat[grepl("ectotherm",df$isolation_source) & !(df$habitat == "Other")] <- "Ectotherm"

#Set all species with no isolation_source/habitat information to habitat = "Other"
df$habitat[is.na(df$habitat)] <- "Other"


####################################################
# Temperature adjusted growth rate using residuals #
####################################################

# With arrhenius temperature

df <- df %>% mutate(arrhenius_tmp = ifelse(!is.na(growth_tmp), 1/(growth_tmp + 273), NA))

#Get data set for creating regression model on growth rate against arrhenius growth temperature
sub <- df %>% filter(!is.na(growth_rate) & !is.na(growth_tmp))

#Get random 75% subset for training model (never use all data for model)
#set.seed(125) 
#sample = sample.split(sub$growth_rate, SplitRatio = .75)
#train = subset(sub, sample == TRUE)

#Ordinary least squares using lm() 
model <- lm(log10(growth_rate)~arrhenius_tmp, data = sub)

#Use model to calculate residuals in original data frame
#measured growth rate minus estimated growth rate = residual
df <- df %>% mutate(
  temp_adjusted_maxgrowth_arrhenius = ifelse(!is.na(growth_rate), log10(growth_rate) - (model$coefficients['arrhenius_tmp']*arrhenius_tmp+model$coefficients['(Intercept)']), NA)) 

# With growth temperature:

#Ordinary least squares using lm() 
model <- lm(log10(growth_rate)~growth_tmp, data = sub)

#Use model to calculate residuals in original data frame
#measured growth rate minus estimated growth rate = residual
df <- df %>% mutate(
  temp_adjusted_maxgrowth = ifelse(!is.na(growth_rate), log10(growth_rate) - (model$coefficients['growth_tmp']*growth_tmp+model$coefficients['(Intercept)']), NA)) 

# Clean up
rm(model,sub)


#####################
# Add strategy data #
#####################

# This requires a bit of preparation:

#Fix organism names in original data frame
strategies$organism[strategies$organism == "Planctomyces maris DSM8797"] <- "Gimesia maris DSM 8797"
strategies$organism[strategies$organism == "Vibrio sp strain Ant-300"] <- "Vibrio sp. ANT-300"
strategies$organism[strategies$organism == "Agromonas oligotrophica"] <- "Bradyrhizobium oligotrophicum"
strategies$organism[strategies$organism == "Caulobacter crescentus CB15"] <- "Caulobacter vibrioides"

#Add tax id to species at strain level (id organism)
st2 <- strategies %>% left_join(nam, by = c("organism"="name_txt"))
#Get species tax id
st3 <- st2 %>% inner_join(tax[,c("tax_id","species_tax_id")], by = "tax_id")

#Keep only relevant information
st4 <- st3 %>% select(organism,type,reference,species_tax_id) %>% 
  rename(strategy_org_name = organism, strategy = type, strategy_ref = reference) %>% 
  distinct(species_tax_id, .keep_all = TRUE)
  
#Attach to main
df <- df %>% left_join(st4, by = "species_tax_id")

#Clean up
rm(strategies,st2,st3,st4)

#Add "Unknown" strategy to rest
df <- df %>% mutate(strategy = ifelse(is.na(strategy),"Unknown",strategy))

#Convert strategy names to first letter upper case
df$strategy <- str_to_title(df$strategy, locale = "en")

#Add short names for labels
# df$short_name <- NA
# for(i in 1:nrow(df)) {
#   df$short_name[i] <- paste(c(substring(word(df$species[i],1),1,1),substring(word(df$species[i],2),1,3)), collapse = ".")
# }


##############################
# Add intracellular grouping #
##############################

df$intracellular <- NA
for(i in 1:nrow(intr)) {
  
  if(intr$phyl_level[i] == "order") {
    df$intracellular[df$order == intr$name[i]] <- "1"
  } else if (intr$phyl_level[i] == "genus") {
    df$intracellular[df$genus == intr$name[i]] <- "1"
  } else if (intr$phyl_level[i] == "species") {
    df$intracellular[df$species == intr$name[i]] <- "1"
  }
}
rm(intr)

df$mycoplasma <- NA
for(i in 1:nrow(myco)) {
    
  if(myco$phyl_level[i] == "order") {
    df$mycoplasma[df$order == myco$name[i]] <- "1"
  } else if (myco$phyl_level[i] == "genus") {
    df$mycoplasma[df$genus == myco$name[i]] <- "1"
  } else if (myco$phyl_level[i] == "species") {
    df$mycoplasma[df$species == myco$name[i]] <- "1"
  }
}
rm(myco)

df <- ungroup(df)

#Add intracellular as habitat
df$habitat[!is.na(df$intracellular)] <- "Intracellular"

#Add myco to intracellular habitat
df$habitat[!is.na(df$mycoplasma)] <- "Intracellular"


############
# Finalise #
############

# Change factor levels of habitat for plots

#Change factor levels of habitat
df$habitat <- factor(df$habitat, levels = c("Fresh","Marine","Soil","Thermal","Endotherm","Ectotherm","Intracellular","Other"))

#Convert genome size to Mbp
df <- df %>% mutate(genome_size = genome_size/1000000) 

#Save full data (with intracellular and myco)
full <- df

#Remove intracellular from main data frame
df <- df %>% filter(is.na(intracellular)) %>% 
  filter(is.na(mycoplasma))

#Remove intracellular as level from habitat in df
df$habitat <- droplevels(df$habitat)

# Clean up
#rm(nam,tax)