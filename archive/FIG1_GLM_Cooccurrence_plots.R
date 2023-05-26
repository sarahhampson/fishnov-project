# Title: GLM Co-occurence plots
# Author: Sarah Hampson

# 0.  Load necessary dataframes and modify --------------------------------

# Read results csvs
full_results_df <- read.csv("./detection results/full_results_df.csv")
cooccur_results_df_site <- read.csv("./inputs/model_data/cooccur_glm_site_df.csv")

# Add time series length, years passed and year position to results df
full_results_temp <- lapply(split(full_results_df, f=full_results_df$sqID),
                            function(x){
                              x$TS_Length <- nrow(x)
                              x$Yrs_passed <- x$bins - min(x$bins)
                              x$Yr_pos <- 1:nrow(x)
                              return(x)
                            })
full_results_df <- do.call("rbind", full_results_temp)

# Use scale function to be able to compare effect sizes
full_results_df$bin.lagS <- scale(full_results_df$bin.lag)[,1]
full_results_df$Yrs_passedS <- scale(full_results_df$Yrs_passed)[,1]
full_results_df$logYr_posS <- scale(log(full_results_df$Yr_pos))[,1]

# 1. Create co-occurrence models ------------------------------------------

# Time point model of novelty cooccurrence using full results df
cooccur.glm.tp <- glmer(funcnovel ~ taxnovel + bin.lagS + logYr_posS
                        + (1|Site/Quarter) + (1|Country/HydroBasin),
                        data = full_results_df,
                        family=binomial)

# Site level model of novelty cooccurrence using full results df
cooccur.glm.site <- glmer(funcnovel ~ taxnovel + bincountS + mean_ts_fillS
                          + (1|Country/HydroBasin), 
                          data = cooccur_results_df_site,
                          family = binomial)


# 2. Plot effect sizes ----------------------------------------------------

# Extract effect size data
tp_effects <- as.data.frame(effects::effect(term= "taxnovel", mod= cooccur.glm.tp))
site_effects <- as.data.frame(effects::effect(term= "taxnovel", mod= cooccur.glm.site))

# Look at these df
tp_effects
site_effects

# Install sjplot package
install.packages("sjPlot")
library(sjPlot)

# Plot the effects 
tp_effects_plot <- plot_model(cooccur.glm.tp, title="", show.values=T) + theme_sjplot()
site_effects_plot <- plot_model(cooccur.glm.site, title="", show.values=T, colors="steelblue") + theme_sjplot()

# View and save plots
tp_effects_plot
ggsave(
  filename="cooccur_tp_effectsize_plot", plot = last_plot(),
  device = "pdf", path = "./figures/GLM",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)
site_effects_plot
ggsave(
  filename="cooccur_site_effectsize_plot", plot = last_plot(),
  device = "pdf", path = "./figures/GLM",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 3. Plot binary probabilities --------------------------------------------

# Plot time point probabilities
tp_glm_cooccur_plot <- ggplot() + 
  #2
  #geom_point(data=L.site.cooccur, aes(taxnov, log(funcnov))) + 
  #3
  geom_point(data=tp_effects, aes(x=taxnovel, y=fit), color="black") +
  #4
  #geom_line(data=tp_effects, aes(taxnovel, y=fit), color="red") +
  #5
  #geom_ribbon(data= x_taxnov, aes(x=taxnov, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  geom_errorbar(data=tp_effects, aes(taxnovel, ymin=lower, ymax=upper), width=0.08, linetype=1, colour="black") +
  theme_classic() +
  #scale_x_continuous(breaks=c(0,1), limits=c(0,1)) +
  #6
  labs(x="Taxonomic novelty?", y="Probability of Size Novelty Occurence")

# Plot site level probabilities
site_glm_cooccur_plot <- ggplot() + 
  #2
  #geom_point(data=L.site.cooccur, aes(taxnov, log(funcnov))) + 
  #3
  geom_point(data=site_effects, aes(x=taxnovel, y=fit), color="black") +
  #4
  #geom_line(data=tp_effects, aes(taxnovel, y=fit), color="red") +
  #5
  #geom_ribbon(data= x_taxnov, aes(x=taxnov, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  geom_errorbar(data=site_effects, aes(taxnovel, ymin=lower, ymax=upper), width=0.08, linetype=1, colour="black") +
  theme_classic() +
  #scale_x_continuous(breaks=c(0,1), limits=c(0,1)) +
  #6
  labs(x="Taxonomic novelty?", y="Probability of Functional Novelty Occurence")

# View and save plots
tp_glm_cooccur_plot
ggsave(
  filename="cooccur_tp_glm_plot", plot = last_plot(),
  device = "pdf", path = "./figures/GLM",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)
site_glm_cooccur_plot
ggsave(
  filename="cooccur_site_glm_plot", plot = last_plot(),
  device = "pdf", path = "./figures/GLM",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 4.  Cleanup -------------------------------------------------------------

rm(full_results_df, cooccur_results_df_site, full_results_temp, cooccur.glm.tp,
   cooccur.glm.site, tp_effects, site_effects, tp_effects_plot, site_effects_plot,
   tp_glm_cooccur_plot, site_glm_cooccur_plot)
