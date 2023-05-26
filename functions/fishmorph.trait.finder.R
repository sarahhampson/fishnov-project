# A function which subset trait data from fishmorph based on a species list

# Define function
fishmorph.trait.finder <- function(traits, fishmorphdata, specieslist){
  
  # Verification variables
  #fishmorphdata <- fishmorph_data
  #specieslist <- species_list
  #traits <- traitlist.FM
  
  # Convert traits to a vector
  traits <- as.character(traits)
  
  # Subset fishmorph data to only species occurding in species list
  fishmorphdata <- subset(fishmorphdata, fishmorphdata$Species %in% specieslist$Species)
  
  # Subset fishmorph trait data to only traits of interest and species name
  trait_data <- fishmorphdata %>% group_by(Species) %>%
    dplyr::select(all_of(traits))
  
  # Return trait values and species
  return(trait_data)
}
