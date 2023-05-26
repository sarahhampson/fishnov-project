# Title: Plots for trophic level changes
# Author: Sarah Hampson


# 0.  Load necessary dataframes -------------------------------------------

# All community troph
comms_troph_df <- read.csv("./inputs/model_data/trait_models/troph_df.csv")

# Taxonomic troph change model summary df
tax_troph_plot_df <- read.csv("./outputs/troph/tax_lm_troph.csv")
tax_Abtroph_plot_df <- read.csv("./outputs/troph/tax_lm_Abtroph.csv")
tax_sdtroph_plot_df <- read.csv("./outputs/troph/tax_lm_sdtroph.csv")
troph_troph_plot_df <- read.csv("./outputs/troph/troph_lm_troph.csv")
troph_Abtroph_plot_df <- read.csv("./outputs/troph/troph_lm_Abtroph.csv")
troph_sdtroph_plot_df <- read.csv("./outputs/troph/troph_lm_sdtroph.csv")

# 1.  Subset and clean dataframes ---------------------------------------------------

# Take only intercept and novel=True estimates
tax_troph_plot_df <- tax_troph_plot_df[1:2,]
tax_Abtroph_plot_df <- tax_Abtroph_plot_df[1:2,]
tax_sdtroph_plot_df <- tax_sdtroph_plot_df[1:2,]
troph_troph_plot_df <- troph_troph_plot_df[1:2,]
troph_Abtroph_plot_df <- troph_Abtroph_plot_df[1:2,]
troph_sdtroph_plot_df <- troph_sdtroph_plot_df[1:2,]

# Absolute change dataframes need to be exponentially transformed
# but only transform the mean estimate and intervals
tax_Abtroph_plot_df <- exp(tax_Abtroph_plot_df[, 3:7])
troph_Abtroph_plot_df <- exp(troph_Abtroph_plot_df[, 3:7])

# Put in a new "cat" variable for "Non novel" and "Novel" rows
tax_troph_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_Abtroph_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_sdtroph_plot_df$taxnovel <- c("FALSE", "TRUE")
troph_troph_plot_df$trophnovel <- c("FALSE", "TRUE")
troph_Abtroph_plot_df$trophnovel <- c("FALSE", "TRUE")
troph_sdtroph_plot_df$trophnovel <- c("FALSE", "TRUE")

# 2.  Absolute troph change plots ------------------------------------------

# Taxonomic
ggplot(data = tax_Abtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, Abtraitchange), trim=T) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression("|" *Delta*" CWM trophic level |"))
# Save plot
ggsave(
  filename="taxAbtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Trophic novelty
ggplot(data = troph_Abtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=trophnovel, Abtraitchange), trim=T) +
  geom_point(aes(trophnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(trophnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(trophnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=trophnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Trophic novelty")) +
  labs(x="Community category", y=expression("|" *Delta*" CWM trophic level |"))
# Save plot
ggsave(
  filename="trophAbtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 3.  Raw troph change plots -----------------------------------------------

# Taxonomic
ggplot(data = tax_troph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, traitchange), trim=F) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression(Delta*" CWM trophic level"))
# Save plot
ggsave(
  filename="taxrawtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Trophic level
ggplot(data = troph_troph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=trophnovel, traitchange), trim=F) +
  geom_point(aes(trophnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(trophnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(trophnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=trophnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Trophic novelty")) +
  labs(x="Community category", y=expression(Delta*" community mean troph"))
# Save plot
ggsave(
  filename="trophrawtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 4.  Standard deviation change plots -------------------------------------

# Taxonomic
ggplot(data = tax_sdtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, sdchange), trim=F) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression(Delta*" "*sigma*" trophic level"))
# Save plot
ggsave(
  filename="taxsdtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# Trophic novelty
ggplot(data = troph_sdtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=trophnovel, sdchange), trim=F) +
  geom_point(aes(trophnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(trophnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(trophnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=trophnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Trophic novelty")) +
  labs(x="Community category", y=expression(Delta*" "*sigma*" trophic level"))
# Save plot
ggsave(
  filename="trophsdtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# 5.  Cleanup -------------------------------------------------------------

rm(comms_troph_df, troph_Abtroph_plot_df, troph_glm_site_df, troph_sdtroph_plot_df,
   troph_troph_plot_df, tax_Abtroph_plot_df, tax_glm_site_df, tax_sdtroph_plot_df,
   tax_troph_plot_df)
