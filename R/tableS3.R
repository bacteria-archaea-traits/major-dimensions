#Table S3

sub <- df %>% filter(!is.na(habitat) & !is.na(temp_adjusted_maxgrowth) & !is.na(d1_mid)) %>%
  mutate(d1_mid_log10 = log10(d1_mid))

m1 <- lm(temp_adjusted_maxgrowth ~ poly(d1_mid_log10, 2) + habitat, sub)

store <- data.frame("Model" = character(), "Coefficient-SE" = character(), "Probability P" = character())

for(i in 1:length(m1$coefficients)) {
  
  if(coef(summary(m1))[i,"Pr(>|t|)"] < 0.001) {
    stars = "***"
  } else if(coef(summary(m1))[i,"Pr(>|t|)"] < 0.01) {
    stars = "**"
  } else if(coef(summary(m1))[i,"Pr(>|t|)"] < 0.05) {
    stars = "*"
  } else {
    stars = ""
  }
  
  store <- rbind(store,data.frame(
    Model = names(m1$coefficients[i]), 
    Coefficient.SE = sprintf("%s +/- %s", round(coef(summary(m1))[i,"Estimate"], digits = 3), round(coef(summary(m1))[i,"Std. Error"], digits = 3)), 
    Probability.P = sprintf("%s %s",coef(summary(m1))[i,"Pr(>|t|)"], stars)))
}

write.csv(store, "output/stats/tableS3.csv", row.names = FALSE)


