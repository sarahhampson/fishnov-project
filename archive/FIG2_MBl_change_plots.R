# Title: Plots for MBl trait changes
# Author: Sarah Hampson


# 0.  Load necessary dataframes -------------------------------------------

# All community size change data
comms_size_df <- read.csv("./inputs/model_data/trait_models/MBl_df.csv")

# Taxonomic size change model summary df
tax_size_plot_df <- read.csv("./outputs/MBl/tax_lm_size.csv")
tax_Absize_plot_df <- read.csv("./outputs/MBl/tax_lm_Absize.csv")
tax_sdsize_plot_df <- read.csv("./outputs/MBl/tax_lm_sdsize.csv")
func_size_plot_df <- read.csv("./outputs/MBl/func_lm_size.csv")
func_Absize_plot_df <- read.csv("./outputs/MBl/func_lm_Absize.csv")
func_sdsize_plot_df <- read.csv("./outputs/MBl/func_lm_sdsize.csv")

# 1.  Subset and clean dataframes ---------------------------------------------------

# Take only intercept and novel=True estimates
tax_size_plot_df <- tax_size_plot_df[1:2,]
tax_Absize_plot_df <- tax_Absize_plot_df[1:2,]
tax_sdsize_plot_df <- tax_sdsize_plot_df[1:2,]
func_size_plot_df <- func_size_plot_df[1:2,]
func_Absize_plot_df <- func_Absize_plot_df[1:2,]
func_sdsize_plot_df <- func_sdsize_plot_df[1:2,]

# Absolute change dataframes need to be exponentially transformed
# but only transform the mean estimate and intervals
tax_Absize_plot_df <- exp(tax_Absize_plot_df[, 3:7])
func_Absize_plot_df <- exp(func_Absize_plot_df[, 3:7])

# Put in a new "cat" variable for "Non novel" and "Novel" rows
tax_size_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_Absize_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_sdsize_plot_df$taxnovel <- c("FALSE", "TRUE")
func_size_plot_df$funcnovel <- c("FALSE", "TRUE")
func_Absize_plot_df$funcnovel <- c("FALSE", "TRUE")
func_sdsize_plot_df$funcnovel <- c("FALSE", "TRUE")

# 2.  Absolute size change plots ------------------------------------------

# Taxonomic
ggplot(data = tax_Absize_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=taxnovel, Abtraitchange), trim=T) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_size_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression("| " *Delta*" CWM size |"))
# Save plot
ggsave(
  filename="taxAbsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# Functional
ggplot(data = func_Absize_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=funcnovel, Abtraitchange), trim=T) +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_size_df, aes(x=funcnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Size novelty")) +
  labs(x="Community category", y=expression("| " *Delta*" CWM size |"))
# Save plot
ggsave(
  filename="funcAbsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# 3.  Raw size change plots -----------------------------------------------

# Taxonomic
ggplot(data = tax_size_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=taxnovel, traitchange), trim=T) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_size_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression(Delta*" CWM size"))
# Save plot
ggsave(
  filename="taxrawsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# Functional
ggplot(data = func_size_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=funcnovel, traitchange), trim=T) +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_size_df, aes(x=funcnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Size novelty")) +
  labs(x="Community category", y=expression(Delta*" CWM size"))
# Save plot
ggsave(
  filename="funcrawsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# 4.  Standard deviation change plots -------------------------------------

# Taxonomic
ggplot(data = tax_sdsize_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=taxnovel, sdchange), trim=F) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_size_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression(Delta* " standard deviation in size"))
# Save plot
ggsave(
  filename="taxsdsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# Functional
ggplot(data = func_sdsize_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=funcnovel, sdchange), trim=F) +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_size_df, aes(x=funcnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Functional novelty")) +
  labs(x="Community category", y=expression(Delta*" standard deviation in size"))
# Save plot
ggsave(
  filename="funcsdsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# 5.  Cleanup -------------------------------------------------------------

rm(comms_size_df, func_Absize_plot_df, func_glm_site_df, func_sdsize_plot_df,
   func_size_plot_df, tax_Absize_plot_df, tax_glm_site_df, tax_sdsize_plot_df,
   tax_size_plot_df)



