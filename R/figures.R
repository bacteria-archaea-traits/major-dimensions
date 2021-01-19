######################################
# Figures for major dimensions paper #
######################################

#######
# Fig 1a "Species max growth rate in relation to genome size - with copios and olligos"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(genome_size))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = genome_size, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = genome_size, y = temp_adjusted_maxgrowth, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours) + 
  scale_size_manual(values = c(2,2,1.5), guide = "none") +
  scale_x_log10() + 
  annotation_logticks(sides="lb", short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(1,"mm")) +
  basic_layout + 
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "Genome size (Mbp)", y = growth_rate_norm_axis_text, fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "Fig1a.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig 1b "Species max growth rate in relation to cell radial diameter"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) 

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = d1_mid, y = temp_adjusted_maxgrowth, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours) + 
  scale_size_manual(values = c(2,2,1.5), guide = "none") +
  scale_x_log10() + 
  annotation_logticks(sides="lb", short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(1,"mm")) +
  basic_layout + 
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "Mean diameter (µm)", y = growth_rate_norm_axis_text, fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "Fig1b.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig 1c "Species genome size in relation to cell radial diameter"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid)) 

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = genome_size)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,0.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = d1_mid, y = genome_size, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours) + 
  scale_size_manual(values = c(2,2,1.5), guide = "none") +
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks(sides="lb", short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(1,"mm")) +
  basic_layout + 
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "Mean diameter (µm)", y = "Genome size (Mbp)", fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "Fig1c.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig 2 "Maximum growth rate in relation to cell diameter, with polynomial fits separated by habitat"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) %>%
  mutate(d1_mid_log10 = log10(d1_mid), growth_rate_log10 = temp_adjusted_maxgrowth)

#NOTE IN ABOVE: temp_adjusted_maxgrowth is derived from log10 scale of residuals, so should NOT be transformed again

sub$habitat <- factor(sub$habitat)

#To be able to directly calculate values for Y we need to output RAW polynomials
#NOT orthogonal polynomials (which are needed for the stats, hence raw = TRUE)
mod <- lm(growth_rate_log10 ~ poly(d1_mid_log10, 2, raw = TRUE) + habitat, data = sub)
summary(mod)

#Get variables
intercept <- mod$coefficients["(Intercept)"]
poly1 <- mod$coefficients["poly(d1_mid_log10, 2, raw = TRUE)1"]
poly2 <- mod$coefficients["poly(d1_mid_log10, 2, raw = TRUE)2"]

#y = intercept + habitat + poly1 * x + poly2 * x^2

#Create data frame with valuse for each of these regressions

#Get names of coefficients in model
coefs <- mod$coefficients[grepl("habitat",names(mod$coefficients))]

#Add Intercept to vector of names - this resembles the missing habitat
coefs <- c(coefs,mod$coefficients["(Intercept)"])

#Create data frame to hold results
polys <- data.frame("habitat" = as.character(), x = as.numeric(), y = as.numeric())

#Create data values for each polynomial for plotting

#Loop through habitats
for(i in 1:length(coefs)) {
  
  #Get name of habitat from coefficient name
  #(i.e. remove 'habitat' from coefficient name)
  name <- sub("habitat","",names(coefs[i]))
  
  #Get coefficient value for this habitat
  # + deal with special case (intercept)
  if(!(name == "(Intercept)")) {
    coef <- as.numeric(coefs[i])
  } else {
    #If this is the intercept (=Fresh), don't double up by adding coef = intercept
    coef <- 0
  }
  
  #Calculate y based on a range of x values that covers the d1_mid range
  for(a in 1:60) {
    
    #Since the model was created on log10 transformed values, we need to use 
    #similar values here - transform x using log10 before calculating
    x <- log10(a/10)
    
    #Calculate y
    y <- poly2*x^2 + poly1*x + intercept + coef
    
    #Add data point to data frame (this is a bit messy..)
    new <- as.data.frame(t(c(name,x,y)), stringsAsFactors = FALSE)
    names(new) <- c("habitat","x","y")
    new$habitat <- as.character(new$habitat)
    new$x <- as.numeric(new$x)
    new$y <- as.numeric(new$y)
    polys <- polys %>% bind_rows(new)
    
    #Loop to next x value for this habitat
  }
  #Loop to next habitat
}

#Change "(Intercept)" to "Fresh"
polys$habitat[polys$habitat == "(Intercept)"] <- "Fresh"
polys$habitat <- as.factor(polys$habitat)

#Convert x and y values back to real numbers (from log10) so
#they can be plotted on a normal scale as original data points
polys$x <- 10^polys$x
#polys$y <- 10^polys$y #When using temp_adjusted_maxgrowth the original data was log10 so no need to transform back

#Check order of peaks in each habitat
polys %>% group_by(habitat) %>% summarise(max = max(y)) %>% arrange(max)

#Arrange factor levels so colours match in plot
polys$habitat <- as.factor(polys$habitat)

#Fitted with polynomials
sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_point(data = sub2, aes(x = d1_mid, y = temp_adjusted_maxgrowth, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  
  geom_smooth(data = polys[polys$habitat == "Fresh",], aes(x = x, y = y), method = "loess", colour = "#A6CEE3") +
  geom_smooth(data = polys[polys$habitat == "Marine",], aes(x = x, y = y), method = "loess", colour = "#1F78B4") +
  geom_smooth(data = polys[polys$habitat == "Soil",], aes(x = x, y = y), method = "loess", colour = "#33A02C") +
  geom_smooth(data = polys[polys$habitat == "Thermal",], aes(x = x, y = y), method = "loess", colour = "#E31A1C") +
  geom_smooth(data = polys[polys$habitat == "Endotherm",], aes(x = x, y = y), method = "loess", colour = "#FB9A99") +
  geom_smooth(data = polys[polys$habitat == "Ectotherm",], aes(x = x, y = y), method = "loess", colour = "#FF7F00") +
  geom_smooth(data = polys[polys$habitat == "Other",], aes(x = x, y = y), method = "loess", colour = "#BDBDBD") +
  
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours) + 
  scale_colour_manual(values = colours) + 
  scale_size_manual(values = c(2,2,1.5), guide = "none") +
  scale_x_log10(limits = c(0.1,6)) + 
  annotation_logticks(sides="lb", short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(1,"mm")) +
  basic_layout + 
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )+
  guides(fill = guide_legend(override.aes = list(shape = 21)), colour = FALSE) +
  labs(x = "Mean diameter (µm)", y = growth_rate_norm_axis_text, fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "Fig2.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig 3 "Species genome size in relation to cell radial diameter - SPLIT"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid))

#Here, split Marine and Fresh into water and sediment as well
sub$habitat <- as.character(sub$habitat)
sub$habitat[grepl("sediment", sub$isolation_source) & sub$habitat == "Marine"] <- "Marine sediment"
sub$habitat[grepl("sediment", sub$isolation_source) & sub$habitat == "Fresh"] <- "Fresh sediment"
sub$habitat[sub$habitat == "Marine"] <- "Marine water"
sub$habitat[sub$habitat == "Fresh"] <- "Fresh water"

#Reorder factors
sub$habitat <- factor(sub$habitat, levels = c("Fresh water","Fresh sediment","Marine water","Marine sediment","Soil","Thermal","Endotherm","Ectotherm","Other"))
#Create colour list with 9 elements
colours_sub <- c(colours_raw[1],"Darkgrey",colours_raw[2],"Darkgrey",colours_raw[4],colours_raw[6],colours_raw[5],colours_raw[8],"#BDBDBD")


sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = genome_size)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,0.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = d1_mid, y = genome_size, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours_sub) + 
  scale_size_manual(values = c(2,2,1.5), guide = "none") +
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks(sides="lb", short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(1,"mm")) +
  basic_layout + 
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(.~habitat) +
  labs(x = "Mean diameter (µm)", y = "Genome size (Mbp)", fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "Fig3.png", plot = p, device = "png", path = save_path, units = "cm", width = 18, height = 18, dpi = 600, limitsize = TRUE)
dev.off()
