#Table 1

store <- data.frame("Trait" = character(), "Genome size" = character(), "Cell radial diameter" = character(), "rRNA operon copy number" = character())

#Growth rate against genome size, diameter and rRNA 
s1 <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(genome_size))
m1 <- lm(temp_adjusted_maxgrowth~log10(genome_size), data = s1)

s2 <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) 
m2 <- lm(temp_adjusted_maxgrowth~log10(d1_mid), data = s2)

s3 <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(rRNA16S_genes))
m3 <- lm(temp_adjusted_maxgrowth~log10(rRNA16S_genes), data = s3)

store <- rbind(store, data.frame(trait = "Normalised max growth rate", 
                                 Genome.size = as.character(sprintf("%s (df=%s)",formatC(summary(m1)$r.squared, digits = 3, format = "g"), summary(m1)$df[2])), 
                                 Cell.radial.diameter = as.character(sprintf("%s (df=%s)", formatC(summary(m2)$r.squared, digits = 3, format = "g"), summary(m2)$df[2])),
                                 rRNA.operon.copy.number = as.character(sprintf("%s (df=%s)", formatC(summary(m3)$r.squared, digits = 3, format = "g"), summary(m3)$df[2]))
))

#Genome size against diameter and rRNA 
s2 <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(d1_mid)) 
m2 <- lm(log10(genome_size)~log10(d1_mid), data = s2)

s3 <- df %>% filter(!is.na(habitat) & !is.na(genome_size) & !is.na(rRNA16S_genes))
m3 <- lm(log10(genome_size)~log10(rRNA16S_genes), data = s3)

store <- rbind(store, data.frame(trait = "Genome size", 
                                 Genome.size = NA, 
                                 Cell.radial.diameter = as.character(sprintf("%s (df=%s)", formatC(summary(m2)$r.squared, digits = 3, format = "g"), summary(m2)$df[2])),
                                 rRNA.operon.copy.number = as.character(sprintf("%s (df=%s)", formatC(summary(m3)$r.squared, digits = 3, format = "g"), summary(m3)$df[2]))
))


#Cell radial diameter against rRNA 
s3 <- df %>% filter(!is.na(habitat) & !is.na(d1_mid) & !is.na(rRNA16S_genes))
m3 <- lm(log10(d1_mid)~log10(rRNA16S_genes), data = s3)

store <- rbind(store, data.frame(trait = "Genome size", 
                                 Genome.size = NA, 
                                 Cell.radial.diameter = NA,
                                 rRNA.operon.copy.number = as.character(sprintf("%s (df=%s)", formatC(summary(m3)$r.squared, digits = 3, format = "g"), summary(m3)$df[2]))
))


write.csv(store, "output/stats/table1.csv", row.names = FALSE)

