# Title: Species information
# Author: Sarah Hampson

# 0.  Load necessary dataframes -------------------------------------------

# Load full results csv
trait_data <- read.csv("./inputs/data/trait_data.csv")

# 1.  Make species list ---------------------------------------------------

species_list <- trait_data$Species

# 2.  Find genus and family level info ------------------------------------

# Create an empty data frame to store the results
species_info_df <- data.frame(Species = character(),
                        Genus = character(),
                        Family = character(),
                        stringsAsFactors = FALSE)

# Loop through the species list and extract genus and family information
for (i in seq_along(species_list)) {
  species_info <- species(species_list[i], fetch = "Genus")
  genus_name <- species_info$Genus[1]
  family_name <- species_info$Family[1]
  
  # add the results to the data frame
  species_info_df[i, "Species"] <- species_list[i]
  species_info_df[i, "Genus"] <- genus_name
  species_info_df[i, "Family"] <- family_name
}

# print the resulting data frame
print(species__df)

