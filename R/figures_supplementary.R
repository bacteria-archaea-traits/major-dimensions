###########################################
# Supp figures for major dimensions paper #
###########################################


#######
# Fig S1 "Coding genese vs genome size"
#######

sub <- full %>% filter(!is.na(habitat) & !is.na(coding_genes) & !is.na(genome_size))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = genome_size, y = coding_genes)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_point(data = sub2, aes(x = genome_size, y = coding_genes, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours_with_intra) + 
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
  labs(x = "Genome size (Mbp)", y = "Coding genes", fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "FigS1.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S2 "Growth rate vs genome size including intracellular"
#######

sub <- full %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(genome_size))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = genome_size, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = genome_size, y = temp_adjusted_maxgrowth, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours_with_intra) + 
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

ggsave(filename = "FigS2.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S3 "Growth rate vs mean diameter including intracellular"
#######

sub <- full %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) 

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = d1_mid, y = temp_adjusted_maxgrowth, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours_with_intra) + 
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

ggsave(filename = "FigS3.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S4 "Genome size vs mean diameter including intracellular"
#######

sub <- full %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = genome_size)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,0.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = d1_mid, y = genome_size, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours_with_intra) + 
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

ggsave(filename = "FigS4.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S5 "Growth rate vs genome size factored by kingdom"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(genome_size))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = genome_size, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = superkingdom, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = genome_size, y = temp_adjusted_maxgrowth, shape = strategy, fill = superkingdom, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = two_factor) + 
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
  labs(x = "Genome size (Mbp)", y = growth_rate_norm_axis_text, fill = "Kingdom", shape = "Strategy")
p

ggsave(filename = "FigS5.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S6 "Growth rate vs mean diameter factored by kingdom"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) 

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = superkingdom, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = d1_mid, y = temp_adjusted_maxgrowth, shape = strategy, fill = superkingdom, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = two_factor) + 
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
  labs(x = "Mean diameter (µm)", y = growth_rate_norm_axis_text, fill = "Kingdom", shape = "Strategy")
p

ggsave(filename = "FigS6.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S7 "Genome size vs mean diameter factored by kingdom"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = d1_mid, y = genome_size)) + 
  geom_point(data = sub[sub$superkingdom == "Bacteria",], aes(fill = superkingdom, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_point(data = sub[sub$superkingdom == "Archaea",], aes(fill = superkingdom, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(data = sub, show.legend = FALSE, h = c(0.5,0.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = d1_mid, y = genome_size, shape = strategy, fill = superkingdom, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = two_factor) + 
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
  #facet_wrap(.~habitat) +
  labs(x = "Mean diameter (µm)", y = "Genome size (Mbp)", fill = "Kingdom", shape = "Strategy")
p

ggsave(filename = "FigS7.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S8 "Growth rate vs rRNA operons"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(rRNA16S_genes))

sub2 <- sub %>% filter(!(strategy == "Unknown"))
sub <- sub %>% filter(strategy == "Unknown")

p <- ggplot(sub, aes(x = rRNA16S_genes, y = temp_adjusted_maxgrowth)) + 
  geom_jitter(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565", width = 0.03) +
  geom_density_2d(show.legend = FALSE, h = c(1,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = rRNA16S_genes, y = temp_adjusted_maxgrowth, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
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
  labs(x = "16S rRNA operons", y = growth_rate_norm_axis_text, fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "FigS8.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S9 "Genome size vs rRNA operons"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(rRNA16S_genes))

sub2 <- sub %>% filter(!(strategy == "Unknown"))
sub <- sub %>% filter(strategy == "Unknown")

p <- ggplot(sub, aes(x = rRNA16S_genes, y = genome_size)) + 
  geom_jitter(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565", width = 0.04) + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,0.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = rRNA16S_genes, y = genome_size, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours) + 
  scale_size_manual(values = c(2,2,0.7), guide = "none") +
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
  labs(x = "16S rRNA operons", y = "Genome size (Mbp)", fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "FigS9.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S10 "Mean diameter vs rRNA operons"
#######

sub <- df %>% filter(!is.na(habitat) & !is.na(d1_mid) & !is.na(rRNA16S_genes))

sub2 <- sub %>% filter(!(strategy == "Unknown"))
sub <- sub %>% filter(strategy == "Unknown")

p <- ggplot(sub, aes(x = rRNA16S_genes, y = d1_mid)) + 
  geom_jitter(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565", width = 0.04) +
  geom_density_2d(show.legend = FALSE, h = c(0.5,0.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = rRNA16S_genes, y = d1_mid, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours) + 
  scale_size_manual(values = c(2,2,0.7), guide = "none") +
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
  labs(x = "16S rRNA operons", y = "Mean diameter (µm)", fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "FigS10.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S11 "Genome size vs mean diameter factored by shape"
#######

sub <- df %>% filter(!is.na(genome_size) & !is.na(d1_mid) & !is.na(shapeagg)) %>%
  arrange(shapeagg)

sub2 <- sub %>% filter(!(strategy == "Unknown"))

alpha <- c(0.4,0.8)

p <- ggplot(sub, aes(x = d1_mid, y = genome_size)) + 
  geom_point(aes(fill = shapeagg, shape = strategy, size = strategy, alpha = shapeagg), colour = "#656565") + 
  geom_density_2d(aes(colour = shapeagg), h = c(0.5,0.5), linetype = 2, size = 0.5, show.legend = FALSE) +
  geom_point(data = sub2, aes(x = d1_mid, y = genome_size, shape = strategy, fill = shapeagg, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = rev(two_factor)) + 
  scale_colour_manual(values = rev(two_factor), guide = "none") +
  scale_size_manual(values = c(2,2,1.5), guide = "none") +
  scale_alpha_manual(values = alpha, guide = "none") +
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
  labs(x = "Mean diameter (µm)", y = "Genome size (Mbp)", fill = "Shape", shape = "Strategy")
p

ggsave(filename = "FigS11.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S12 "Growth rate vs volume with intracellular"
#######

#NOTE: below selection seem to remove all Oligotrophs

sub <- full %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(volume))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = volume, y = temp_adjusted_maxgrowth)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,1.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = volume, y = temp_adjusted_maxgrowth, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,21)) +
  scale_fill_manual(values = colours_with_intra) + 
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
  labs(x = volume_axis_text, y = growth_rate_norm_axis_text, fill = "Habitat", shape = "Strategy")
p


ggsave(filename = "FigS12.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()


#######
# Fig S13 "Genome size vs cell volume with intracellular"
#######

sub <- full %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(volume))

sub2 <- sub %>% filter(!(strategy == "Unknown"))

p <- ggplot(sub, aes(x = volume, y = genome_size)) + 
  geom_point(aes(fill = habitat, shape = strategy, size = strategy), alpha = 0.8, colour = "#656565") + 
  geom_density_2d(show.legend = FALSE, h = c(0.5,0.5), linetype = 2, colour = "black", size = 0.3) +
  geom_point(data = sub2, aes(x = volume, y = genome_size, shape = strategy, fill = habitat, size = strategy), show.legend = TRUE) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = colours_with_intra) + 
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
  labs(x = volume_axis_text, y = "Genome size (Mbp)", fill = "Habitat", shape = "Strategy")
p

ggsave(filename = "FigS13.png", plot = p, device = "png", path = save_path, units = "cm", width = 12, height = 8, dpi = 600, limitsize = TRUE)
dev.off()