# Title: Final plots for trait changes
# Author: Sarah Hampson

# 0.  Load necessary dataframes -------------------------------------------

install.packages("patchwork")
library(patchwork)

# All community data
comms_PC1_df <- read.csv("./inputs/model_data/trait_models/PC1_df.csv")
comms_PC2_df <- read.csv("./inputs/model_data/trait_models/PC2_df.csv")
comms_size_df <- read.csv("./inputs/model_data/trait_models/MBl_df.csv")
comms_troph_df <- read.csv("./inputs/model_data/trait_models/troph_df.csv")

# Taxonomic trait change model summary dfs
tax_PC1_plot_df <- read.csv("./outputs/multitrait/tax_lm_PC1.csv")
tax_AbPC1_plot_df <- read.csv("./outputs/multitrait/tax_lm_AbPC1.csv")
tax_PC2_plot_df <- read.csv("./outputs/multitrait/tax_lm_PC2.csv")
tax_AbPC2_plot_df <- read.csv("./outputs/multitrait/tax_lm_AbPC2.csv")
multitrait_PC1_plot_df <- read.csv("./outputs/multitrait/multitrait_lm_PC1.csv")
multitrait_AbPC1_plot_df <- read.csv("./outputs/multitrait/multitrait_lm_AbPC1.csv")
multitrait_PC2_plot_df <- read.csv("./outputs/multitrait/multitrait_lm_PC2.csv")
multitrait_AbPC2_plot_df <- read.csv("./outputs/multitrait/multitrait_lm_AbPC2.csv")
tax_size_plot_df <- read.csv("./archive/tax_lm_size.csv")
tax_Absize_plot_df <- read.csv("./archive/tax_lm_Absize.csv")
size_size_plot_df <- read.csv("./archive/func_lm_size.csv")
size_Absize_plot_df <- read.csv("./archive/func_lm_Absize.csv")
tax_troph_plot_df <- read.csv("./outputs/troph/tax_lm_troph.csv")
tax_Abtroph_plot_df <- read.csv("./outputs/troph/tax_lm_Abtroph.csv")
troph_troph_plot_df <- read.csv("./outputs/troph/troph_lm_troph.csv")
troph_Abtroph_plot_df <- read.csv("./outputs/troph/troph_lm_Abtroph.csv")

# 1.  Subset and clean dataframes ---------------------------------------------------

# Take only intercept and novel=True estimates
tax_PC1_plot_df <- tax_PC1_plot_df[1:2,]
tax_AbPC1_plot_df <- tax_AbPC1_plot_df[1:2,]
tax_PC2_plot_df <- tax_PC2_plot_df[1:2,]
tax_AbPC2_plot_df <- tax_AbPC2_plot_df[1:2,]
multitrait_PC1_plot_df <- multitrait_PC1_plot_df[1:2,]
multitrait_AbPC1_plot_df <- multitrait_AbPC1_plot_df[1:2,]
multitrait_PC2_plot_df <- multitrait_PC2_plot_df[1:2,]
multitrait_AbPC2_plot_df <- multitrait_AbPC2_plot_df[1:2,]
tax_size_plot_df <- tax_size_plot_df[1:2,]
tax_Absize_plot_df <- tax_Absize_plot_df[1:2,]
size_size_plot_df <- size_size_plot_df[1:2,]
size_Absize_plot_df <- size_Absize_plot_df[1:2,]
tax_troph_plot_df <- tax_troph_plot_df[1:2,]
tax_Abtroph_plot_df <- tax_Abtroph_plot_df[1:2,]
troph_troph_plot_df <- troph_troph_plot_df[1:2,]
troph_Abtroph_plot_df <- troph_Abtroph_plot_df[1:2,]

# Absolute change dataframes need to be exponentially transformed
# but only transform the mean estimate and intervals
tax_AbPC1_plot_df <- exp(tax_AbPC1_plot_df[, 3:7])
tax_AbPC2_plot_df <- exp(tax_AbPC2_plot_df[, 3:7])
multitrait_AbPC1_plot_df <- exp(multitrait_AbPC1_plot_df[, 3:7])
multitrait_AbPC2_plot_df <- exp(multitrait_AbPC2_plot_df[, 3:7])
tax_Absize_plot_df <- exp(tax_Absize_plot_df[,3:7])
size_Absize_plot_df <- exp(size_Absize_plot_df[,3:7])
tax_Abtroph_plot_df <- exp(tax_Abtroph_plot_df[,3:7])
troph_Abtroph_plot_df <- exp(troph_Abtroph_plot_df[,3:7])

# Put in a new "cat" variable for "Non novel" and "Novel" rows
tax_PC1_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_AbPC1_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_PC2_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_AbPC2_plot_df$taxnovel <- c("FALSE", "TRUE")
multitrait_PC1_plot_df$funcnovel <- c("FALSE", "TRUE")
multitrait_AbPC1_plot_df$funcnovel <- c("FALSE", "TRUE")
multitrait_PC2_plot_df$funcnovel <- c("FALSE", "TRUE")
multitrait_AbPC2_plot_df$funcnovel <- c("FALSE", "TRUE")
tax_size_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_Absize_plot_df$taxnovel <- c("FALSE", "TRUE")
size_size_plot_df$funcnovel <- c("FALSE", "TRUE")
size_Absize_plot_df$funcnovel <- c("FALSE", "TRUE")
tax_troph_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_Abtroph_plot_df$taxnovel <- c("FALSE", "TRUE")
troph_troph_plot_df$funcnovel <- c("FALSE", "TRUE")
troph_Abtroph_plot_df$funcnovel <- c("FALSE", "TRUE")

# 2.  Make Tax PC plot objects ------------------------------------------

# Taxonomic absolute PC1 change
plot1 <- ggplot(data = tax_AbPC1_plot_df) + 
  geom_violin(data=comms_PC1_df, aes(x=taxnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.9)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM PC1 |")) #x="Community category"
plot1
# Save plot
ggsave(
  filename="taxAbPC1plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Taxonomic raw PC1 change
plot2 <- ggplot(data = tax_PC1_plot_df) + 
  geom_violin(data=comms_PC1_df, aes(x=taxnovel, traitchange), scale="area", trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.8, 0.8)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM PC1")) #x="Community category"
plot2
# Save plot
ggsave(
  filename="taxPC1plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Taxonomic absolute PC2 change
plot3 <- ggplot(data = tax_AbPC2_plot_df) + 
  geom_violin(data=comms_PC2_df, aes(x=taxnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.7)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM PC2 |")) #x="Community category"
plot3
# Save plot
ggsave(
  filename="taxAbPC2plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Taxonomic raw PC2 change
plot4 <- ggplot(data = tax_PC2_plot_df) + 
  geom_violin(data=comms_PC2_df, aes(x=taxnovel, traitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.45, 0.45)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM PC2")) #x="Community category"
plot4
# Save plot
ggsave(
  filename="taxPC2plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 3.  Make Tax size plot objects -----------------------------------------------

# Taxonomic absolute size change
plot5 <- ggplot(data = tax_Absize_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=taxnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.7)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM size |")) #x="Community category"
plot5
# Save plot
ggsave(
  filename="taxAbsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Taxonomic raw size change
plot6 <- ggplot(data = tax_size_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=taxnovel, traitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.2, 0.2)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM size")) #x="Community category"
plot6
# Save plot
ggsave(
  filename="taxsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 4.  Make tax troph plot objects -------------------------------------

# Taxonomic absolute troph change
plot7 <- ggplot(data = tax_Abtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.3)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM trophic level |")) #x="Community category"
# Save plot
ggsave(
  filename="taxAbtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Taxonomic raw size change
plot8 <- ggplot(data = tax_troph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, traitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.2, 0.2)) +
  scale_x_discrete(labels = c("Background", "Novel SpC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM trophic level")) #x="Community category"
# Save plot
ggsave(
  filename="taxtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)


# 5.  Make Multivariate PC plot objects ------------------------------------------

# Multitrait absolute PC1 change
plot9 <- ggplot(data = multitrait_AbPC1_plot_df) + 
  geom_violin(data=comms_PC1_df, aes(x=funcnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.9)) +
  scale_x_discrete(labels = c("Background", "Novel MC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM PC1 |")) #x="Community category"
plot9
# Save plot
ggsave(
  filename="multitraitAbPC1plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Multitrait raw PC1 change
plot10 <- ggplot(data = multitrait_PC1_plot_df) + 
  geom_violin(data=comms_PC1_df, aes(x=funcnovel, traitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.8, 0.8)) +
  scale_x_discrete(labels = c("Background", "Novel MC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM PC1")) #x="Community category"
plot10
# Save plot
ggsave(
  filename="multitraitPC1plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Multitrait absolute PC2 change
plot11 <- ggplot(data = multitrait_AbPC2_plot_df) + 
  geom_violin(data=comms_PC2_df, aes(x=funcnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.7)) +
  scale_x_discrete(labels = c("Background", "Novel MC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM PC2 |")) #x="Community category"
plot11
# Save plot
ggsave(
  filename="multitraitAbPC2plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Multitrait raw PC2 change
plot12 <- ggplot(data = multitrait_PC2_plot_df) + 
  geom_violin(data=comms_PC2_df, aes(x=funcnovel, traitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.45, 0.45)) +
  scale_x_discrete(labels = c("Background", "Novel MC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM PC2")) #x="Community category"
plot12
# Save plot
ggsave(
  filename="multitraitPC2plot", plot = last_plot(),
  device = "pdf", path = "./figures/PC1",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 6.  Make Size size plot objects -----------------------------------------------

# Taxonomic absolute size change
plot13 <- ggplot(data = size_Absize_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=funcnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.7)) +
  scale_x_discrete(labels = c("Background", "Novel SizeC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM size |")) #x="Community category"
# Save plot
ggsave(
  filename="sizeAbsizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# Taxonomic raw size change
plot14 <- ggplot(data = size_size_plot_df) + 
  geom_violin(data=comms_size_df, aes(x=funcnovel, traitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.2, 0.2)) +
  scale_x_discrete(labels = c("Background", "Novel SizeC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM size")) #x="Community category"
# Save plot
ggsave(
  filename="sizesizeplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

# 7.  Make Troph troph plot objects -------------------------------------

plot15 <- ggplot(data = troph_Abtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=trophnovel, Abtraitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 0.3)) +
  scale_x_discrete(labels = c("Background", "Novel TC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression("|" *Delta*" CWM trophic level |")) #x="Community category"
# Save plot
ggsave(
  filename="trophAbtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)

plot16 <- ggplot(data = troph_troph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=trophnovel, traitchange), trim=T, colour="white", fill="lightgrey") +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.05, linetype=1) +
  #geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  theme_classic() + geom_hline(yintercept=0,linetype=2) +
  theme(axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(-0.2, 0.2)) +
  scale_x_discrete(labels = c("Background", "Novel TC")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL, y=expression(Delta*" CWM trophic level")) #x="Community category"
# Save plot
ggsave(
  filename="trophtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 10,
  units = "cm", dpi = 300)


# 8.  Make final whole plot for tax trait change --------------------------

plot_tax <- (plot1 + plot2)/
  (plot3 + plot4)/
  (plot5 + plot6)/
  (plot7 + plot8)

plot_tax + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(face = 'bold'))

ggsave(
  filename="taxtraitsplot", plot = last_plot(),
  device = "pdf", path = "./figures",
  scale = 1, width = 16, height = 25,
  units = "cm", dpi = 300)

# 9.  Make final whole plot for func trait change --------------------------

plot_func <- (plot9 + plot10)/
  (plot11 + plot12)/
  (plot13 + plot14)/
  (plot15 + plot16)

plot_func + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(face = 'bold'))

ggsave(
  filename="functraitsplot", plot = last_plot(),
  device = "pdf", path = "./figures",
  scale = 1, width = 16, height = 25,
  units = "cm", dpi = 300)


# 10.  Make final emergence probability plot ------------------------------

# Get means, CIs and category
means_plot_df1 <- make.glm.df(tax.glm.site, "Species", "plot")[1, 4:7]
means_plot_df2 <- make.glm.df(func.glm.site, "Morphology", "plot")[1, 4:7]
means_plot_df3 <- make.glm.df(size.glm.site, "Size", "plot")[1, 4:7]
means_plot_df4 <- make.glm.df(troph.glm.site, "Trophic Level", "plot")[1, 4:7]

means_plot_df <- rbind(means_plot_df1, means_plot_df2, means_plot_df3, means_plot_df4)
means_plot_df$Category <- factor(c("Species", "Morphology", "Size", "Trophic Level"), 
                                 levels=c("Species", "Morphology", "Size", "Trophic Level"))

# Make plot of means
means_plot <- ggplot(data=means_plot_df, aes(x = Category, y = `Mean Estimate`)) +
  geom_point(aes(x = Category, y = `Mean Estimate`), size = 2) +
  geom_errorbar(aes(ymin = `CI Lower`, ymax = `CI Upper`), width = 0.1, color = "black", size = 0.5) +
  labs(x = "Novel composition", y = "Emergence Probability") + 
  coord_cartesian(ylim=c(0, 0.2)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()
means_plot +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave(
  filename="meanprobsplot", plot = last_plot(),
  device = "pdf", path = "./figures",
  scale = 1, width = 12, height = 12,
  units = "cm", dpi = 300)

#Table 1: Results of the Generalized Linear Model predicting Y from X1 and X2

#| Coefficients | Standard Errors | p-values |
#  | ------------ | --------------- | -------- |
#  | Intercept    | 1.24            | 0.003    |
#  | X1           | 0.12            | 0.05     |
##  | X2           | -0.20           | 0.01     |
  
#Note: The model was fit using the logit link function and the binomial distribution. 
#The degree of freedom was 198. The deviance was 225.5. The AIC was 235.5. The McFadden's pseudo-R2 value was 0.32.


# 10.  Cleanup -------------------------------------------------------------

rm(comms_PC1_df, comms_size_df, comms_troph_df, multitrait_AbPC1_plot_df, 
   multitrait_PC1_plot_df, plot_func, plot_tax, plot1, plot2, plot3, plot4,
   plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, size_Absize_plot_df,
   size_size_plot_df, tax_Absize_plot_df, tax_troph_plot_df, tax_Abtroph_plot_df,
   tax_size_plot_df, troph_Abtroph_plot_df, troph_troph_plot_df)
