# A function which subset trait data from fishbase based on a species list

# Define function
fishbase.trait.finder <- function(specieslist, traits){
  
  # Verification variables
  #fishmorphdata <- fishmorph_data
  #specieslist <- species_list
  traits = "Length"
  
  # Change traits to a character string
  traits = as.character(traits)
    
  # If trait is body length, 
  if("Length" %in% traits){
    # Search length
    trait_data <- fb_tbl("species") %>% 
      mutate(Species = paste(Genus, Species)) %>%
      filter(Species %in% specieslist$Species) %>% 
      dplyr::select(Species, "Length")
    
    }
    
  # If trait is Trophic level (food or diet)
  if("TrophicLevel" %in% traits){
    # Retrieve trophic level data
    trophic.data <- ecology(as.character(unlist(specieslist)), fields = c("Species", "FoodTroph", "DietTroph"))
    
    # Check for NAs
    nrow(trophic.data[is.na(trophic.data$FoodTroph),]) # 130 species without data
    nrow(trophic.data[is.na(trophic.data$DietTroph),]) # 294 species without data
    
    # Replace NAs with the genus-level or family-level average of both food and diet trophic levels
    trophic.data <- trophic_level_FB(as.character(unlist(specieslist)), type = "food items") 
  }
  
  # Create a list of species that DO NOT appear in fishmorph
  #fishbase_species <- specieslist %>%
    #summarise(Species = Species) %>%
    #filter(!Species %in% fishmorphdata$Species)
  
  return(trait_data)
}
