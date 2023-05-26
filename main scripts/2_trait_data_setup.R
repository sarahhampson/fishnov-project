# Title: Trait data setup
# Author: Sarah Hampson


# 1.  Identify species ----------------------------------------------------

# Import survey_data and fishmorph data if not already in environment
survey_data <- read.csv("./inputs/data/survey_data.csv")
fishmorph_data <- read.csv("./inputs/data/fishmorph_data.csv")

# Make list of unique species
species_list <- survey_data %>% dplyr::summarise(Species = unique(Species))

# 2.  Fishmorph traits ----------------------------------------------------

# Make a list of traits that we are interested in from fishmorph
traitlist_FM <- as.list(c("MBl", "BEl", "VEp", "REs", "OGp", "RMl", "BLs", "PFv", "PFs", "CPt"))

# Use fishmorph_trait_finder to select all data from fishmorph
trait_data_FM <- fishmorph.trait.finder(traitlist_FM, fishmorph_data, species_list)

# 3.  Trophic level -------------------------------------------------------

# Retrieve trophic level data from rfishbase
trophic_data_FB <- ecology(species_list$Species, 
                           fields = c("Species", "FoodTroph", "DietTroph"))

# Make list of species missing from troph data
missing_troph_sp <- species_list %>% 
  filter(!Species %in% trophic_data_FB$Species)

# Show more missing species from above troph list and filter from df
missing_troph_sp2 <- trophic_data_FB %>% filter(is.na(FoodTroph), is.na(DietTroph))
trophic_data_FB <- trophic_data_FB %>% filter(!Species %in% missing_troph_sp2$Species)

# Make sure no double ups in troph data from FB by taking mean diet and food troph
trophic_data_FB <- trophic_data_FB %>% 
  group_by(Species) %>% 
  dplyr::summarise(FoodTroph = mean(FoodTroph, na.rm=T),
                   DietTroph = mean(DietTroph, na.rm=T))

# There are 151+37 so 188 species missing, which I found information on via online fishbase

# Read in csv of model based fishbase species without specific trophic information
troph_Est_data <- read_csv('./inputs/raw_data/trophic_level_FBModels.csv', col_names = T)
troph_Est_data$Species <- gsub("_", " ", troph_Est_data$Species)
troph_Est_data <- troph_Est_data %>% 
  rename(troph_level = Trophic_level) %>% 
  mutate(troph_cat = "Est") %>% 
  dplyr::select(-Cat)

# Filter troph diet data first
troph_data_diet <- trophic_data_FB %>% 
  filter(!is.na(DietTroph)) %>% 
  rename(troph_level = DietTroph) %>% 
  dplyr::select(-FoodTroph) %>% 
  mutate(troph_cat = "Diet")

# Filter troph food data
troph_data_food <- trophic_data_FB %>% 
  filter(!Species %in% troph_data_diet$Species) %>% 
  rename(troph_level = FoodTroph) %>% 
  dplyr::select(-DietTroph) %>% 
  mutate(troph_cat = ifelse(troph_level>0, "Food", NA))

# Combine the data into one trophic df
troph_data <- rbind(troph_Est_data, troph_data_diet, troph_data_food)

# 4.  Make full trait dataframe -------------------------------------------

# Put all trait data together and remove unnecessary columns
full_trait_data <- left_join(trait_data_FM, troph_data, by="Species")
full_trait_data <- full_trait_data %>% dplyr::select(-troph_cat)
full_trait_data$troph_level <- as.numeric(full_trait_data$troph_level)

# 5.  Add missing MBl -----------------------------------------------------

# Determine  NA trait gaps (7 species)
NA_trait_species <- full_trait_data %>%
  filter(is.na(MBl))

# Sub in found values for missing trait data
full_trait_data <- full_trait_data %>% 
  mutate(MBl = ifelse(Species == "Hypophthalmus oremaculatus", 44,
                ifelse(Species == "Moxostoma collapsum", 45,
                 ifelse(Species == "Percina nevisense", 77.7,
                  ifelse(Species == "Pterygoplichthys ambrosettii", 50, MBl)))))

# Check NA species again and then remove unecessary df
NA_trait_species <- full_trait_data %>%
  filter(is.na(MBl))

# 6.  Collinearity and coverage -------------------------------------------

# Checking to see whether co-linearity exists between traits
# Using VIF (variance inflation factor)
# If VIF = 1 there is no correlation
# If 1<VIF<5 there is moderate correlation
# If VIF> 5 there is strong correlation

# Install VIF package and define multiple linear regression model
install.packages("car")
library(car)
# Trial all versions of this where each trait is the response
traits <- full_trait_data %>% dplyr::select(-X, -Species)
# Create correlation matrix and save it
cor_mat <- cor(traits, use="pairwise.complete.obs")
write.csv(cor_mat, "./archive/trait_correlation_matrix.csv")
# Create vif values and save it
vif_values <- apply(traits, 2, function(x) vif(lm(x ~ ., data = traits)))
write.csv(vif_values, "./archive/vif_trait_values.csv")

# Define a function to calculate proportion of NA values for traits
trait_na_table <- sapply(trait_data[2:ncol(trait_data)], function(x) {
  prop.table(table(is.na(x)))
})
trait_na_table <- do.call("rbind", trait_na_table)

# Define a function to calculate proportion of NA values for species row
sp_na_table <- trait_data %>% dplyr::select(Species)
sp_na_table$na_count <- apply(trait_data, 1, function(x) {
  sum(is.na(x))
})

# Print the table
print(trait_na_table)
sp_na_table %>% filter(na_count>1) %>% arrange(desc(na_count))

rm(trait_na_table, sp_na_table)

# 7.  Save trait data -----------------------------------------------------

# Create filename and write csv to inputs for later use
filename_traits <- "./inputs/data/trait_data.csv"
write.csv(full_trait_data, filename_traits)

# 8.  Clean up ------------------------------------------------------------

rm(fishmorph_data, species_list, full_trait_data, func_dissims_matrix, 
   troph_data_diet, troph_data_food, survey_data, trait_data_FM, 
   traitlist_FM, troph_Est_data, trophic_data_FB, missing_troph_sp, 
   missing_troph_sp2, colin.model, NA_trait_species, traits, vif_values,
   filename_dissims, filename_traits, troph_data, trait_na_table,
   sp_na_table, cor_mat, vif_values)
