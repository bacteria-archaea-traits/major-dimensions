# Data prep

dat <- df

dat$genome_size_log10 <- log10(dat$genome_size)
dat$d1_mid_log10 <- log10(dat$d1_mid)

# Models
store <- data.frame()

mod1 <- lm(temp_adjusted_maxgrowth ~ d1_mid_log10, data=dat)
temp <- summary(mod1)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod1), r2=temp$r.squared, df=mod1$df))

mod2 <- lm(temp_adjusted_maxgrowth ~ genome_size_log10, data=dat)
temp <- summary(mod2)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod2), r2=temp$r.squared, df=mod2$df))

mod3 <- lm(temp_adjusted_maxgrowth ~ d1_mid_log10 * genome_size_log10, data=dat)
temp <- summary(mod3)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod3), r2=temp$r.squared, df=mod3$df))

mod4 <- lm(temp_adjusted_maxgrowth ~ poly(d1_mid_log10, 2), data=dat[!is.na(dat$d1_mid_log10),])
temp <- summary(mod4)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod4), r2=temp$r.squared, df=mod4$df))

mod5 <- lm(temp_adjusted_maxgrowth ~ poly(genome_size_log10, 2), data=dat[!is.na(dat$genome_size_log10),])
temp <- summary(mod5)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod5), r2=temp$r.squared, df=mod5$df))

mod6 <- lm(temp_adjusted_maxgrowth ~ habitat, data=dat)
temp <- summary(mod6)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod6), r2=temp$r.squared, df=mod6$df))

mod7 <- lm(temp_adjusted_maxgrowth ~ poly(d1_mid_log10, 2) + habitat, data=dat[!is.na(dat$d1_mid_log10),])
temp <- summary(mod7)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod7), r2=temp$r.squared, df=mod7$df))

mod8 <- lm(temp_adjusted_maxgrowth ~ poly(d1_mid_log10, 2) * habitat, data=dat[!is.na(dat$d1_mid_log10),])
temp <- summary(mod8)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod8), r2=temp$r.squared, df=mod8$df))

mod9 <- lm(temp_adjusted_maxgrowth ~ poly(d1_mid_log10, 2) * genome_size_log10 * habitat, data=dat[!is.na(dat$d1_mid_log10),])
temp <- summary(mod9)
store <- rbind(store, data.frame(model=as.character(temp$call)[2], AIC=AIC(mod9), r2=temp$r.squared, df=mod9$df))

write.csv(store, "output/stats/tableS4.csv", row.names = FALSE)
