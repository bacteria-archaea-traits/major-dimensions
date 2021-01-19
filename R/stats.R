#Major dimensions Stats


# Function to simplify code and make output more clear
run_tests <- function(sub, x, y, file, figure, note = "") {
  
  if(!(file == "")) {
    sink(sprintf("output/stats/%s",file), append = TRUE)
  }
  
  print("############################")
  print(figure, quotes = FALSE)
  print("", quotes = FALSE)
  
  #Always log transform x variables
  sub[[x]] <- log10(sub[[x]])
  
  #Don't log transform y-axis residuals
  if(y == "temp_adjusted_maxgrowth" | y == "temp_adjusted_maxgrowth_arrehnius") {
    print(sprintf("y = %s vs x = log10( %s )",y,x), quotes = FALSE)
  } else {
  #Log transform all other y-axis variables
    sub[[y]] <- log10(sub[[y]])
    print(sprintf("y = log10( %s ) vs x = log10( %s )",y,x), quotes = FALSE)
  }
    
  if(!(note == "")) {
    print("", quotes = FALSE)
    print(sprintf("------%s------",note), quotes = FALSE)
  }
  
  # NORMALITY CHECKS
  # par(mfrow=c(1, 2))
  # plot(density(sub[[x]]), main = "Density plot log10(growth_rate)")
  # plot(density(sub[[y]]), main = "Density plot log10(genome_size)")
  # 
  # print("------NORMALITY TEST------")
  # print(sprintf("NORMALITY : %s", x))
  # print(shapiro.test(sub[[x]]))
  # 
  # print(sprintf("NORMALITY : %s", y))
  # print(shapiro.test(sub[[y]]))

  #LINEAR REGRESSIONS
  print("------LINEAR REGRESSION------")
  print(summary(lm(sub[[y]] ~ sub[[x]])))
  
  if(!(file == "")) {
    sink()
  }
}
 

########
#Fig 1a #
########

#No intra
sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(genome_size))

x <- "genome_size"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure 1a"
comment <- "NO INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

########
#Fig 1b #
########

#No intra
sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) 

x <- "d1_mid"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure 1b"
comment <- "NO INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)


########
#Fig 1c #
########

#No intra
sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid)) 

x <- "d1_mid"
y <- "genome_size"
figure <- "Figure 1c"
comment <- "NO INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)


########
#Fig 2 #
########

#Stats already produced 

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) %>%
  mutate(d1_mid_log10 = log10(d1_mid))

mod <- lm(temp_adjusted_maxgrowth ~ poly(d1_mid_log10, 2) + habitat, sub)

sink("output/stats/stats_all_figs.txt", append = TRUE)
print("############################")
print("Figure 2")
print("------NO INTRACELLULAR------")
summary(mod)
sink()


######################
# SUPPLEMENTARY FIGS #
######################

##########
# Fig S1 #
##########

sub <- full %>% filter(!is.na(habitat) & !is.na(coding_genes) & !is.na(genome_size))

x <- "genome_size"
y <- "coding_genes"
figure <- "Figure S1"
comment <- "WITH INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

##########
# Fig S2 #
##########

#With intra
sub <- full %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(genome_size)) 

x <- "genome_size"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure S2"
comment <- "WITH INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

##########
# Fig S3 #
##########

#With intra
sub <- full %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid))

x <- "d1_mid"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure S3"
comment <- "WITH INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)


#########
#Fig S4 #
#########

#With intra
sub <- full %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid)) 

x <- "d1_mid"
y <- "genome_size"
figure <- "Figure S4"
comment <- "WITH INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

##########
# Fig S5 #
##########

sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(temp_adjusted_maxgrowth))

#Bacteria
sub2 <- sub %>% filter(superkingdom == "Bacteria")
x <- "genome_size"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure S5"
comment <- "KINGDOM: BACTERIA"
save_to <- "stats_all_figs.txt"
run_tests(sub2,x,y,save_to,figure,comment)

#Archaea
sub2 <- sub %>% filter(superkingdom == "Archaea")
comment <- "KINGDOM: ARCHAEA"
run_tests(sub2,x,y,save_to,figure,comment)

##########
# Fig S6 #
##########

sub <- df %>% filter(!is.na(habitat) & !is.na(d1_mid) & !is.na(temp_adjusted_maxgrowth))

#Bacteria
sub2 <- sub %>% filter(superkingdom == "Bacteria")
x <- "d1_mid"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure S6"
comment <- "KINGDOM: BACTERIA"
save_to <- "stats_all_figs.txt"
run_tests(sub2,x,y,save_to,figure,comment)

#Archaea
sub2 <- sub %>% filter(superkingdom == "Archaea")
comment <- "KINGDOM: ARCHAEA"
run_tests(sub2,x,y,save_to,figure,comment)

##########
# Fig S7 #
##########

sub <- df %>% filter(!is.na(habitat) & !is.na(d1_mid) & !is.na(genome_size))

#Bacteria
sub2 <- sub %>% filter(superkingdom == "Bacteria")
x <- "d1_mid"
y <- "genome_size"
figure <- "Figure S7"
comment <- "KINGDOM: BACTERIA"
save_to <- "stats_all_figs.txt"
run_tests(sub2,x,y,save_to,figure,comment)

#Archaea
sub2 <- sub %>% filter(superkingdom == "Archaea")
comment <- "KINGDOM: ARCHAEA"
run_tests(sub2,x,y,save_to,figure,comment)

###########
# Fig S8 #
###########

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(rRNA16S_genes))

#With high tmp organisms
x <- "rRNA16S_genes"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure S8"
comment <- "ALL GROWTH TEMPERATURES"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

#No high tmp organisms
sub <- sub %>% filter(growth_tmp <= 50)
comment <- "GROWTH TEMPERATURE <= 50C"
run_tests(sub,x,y,save_to,figure,comment)

##########
# Fig S9 #
##########

sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(rRNA16S_genes)) 

#With high tmp organisms
x <- "rRNA16S_genes"
y <- "genome_size"
figure <- "Figure S9"
comment <- "ALL GROWTH TEMPERATURES"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

#No high tmp organisms
sub <- sub %>% filter(growth_tmp <= 50)
comment <- "GROWTH TEMPERATURE <= 50C"
run_tests(sub,x,y,save_to,figure,comment)

###########
# Fig S10 #
###########

sub <- df %>% filter(!is.na(habitat) & !is.na(d1_mid) & !is.na(rRNA16S_genes)) 

#With high tmp organisms
x <- "rRNA16S_genes"
y <- "d1_mid"
figure <- "Figure S10"
comment <- ""
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

##########
# Fig S11 #
##########

sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid)) 

#Rods
sub2 <- sub %>% filter(shapeagg == "rod")
x <- "d1_mid"
y <- "genome_size"
figure <- "Figure S11"
comment <- "SHAPE: Rods"
save_to <- "stats_all_figs.txt"
run_tests(sub2,x,y,save_to,figure,comment)

#Spheroids
sub2 <- sub %>% filter(shapeagg == "spheroid")
comment <- "SHAPE: Spheroids"
run_tests(sub2,x,y,save_to,figure,comment)

###########
# Fig S12 #
###########

#With intra
sub <- full %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(volume))

x <- "volume"
y <- "temp_adjusted_maxgrowth"
figure <- "Figure S12"
comment <- "WITH INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

#No intra
sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(volume))

comment <- "NO INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

###########
# Fig S13 #
###########

#With intra
sub <- full %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(volume))

x <- "volume"
y <- "genome_size"
figure <- "Figure S13"
comment <- "WITH INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)


#No intra
sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(volume))

comment <- "NO INTRACELLULAR"
save_to <- "stats_all_figs.txt"
run_tests(sub,x,y,save_to,figure,comment)

# END