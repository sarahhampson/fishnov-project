# A function which checks and standardises all species names in RivFishTime
# or other databases according to the FishBase database

standardise.species.names <- function(database){ # Give database with a species column
  
  # Verification variables
  #database <- read.csv("inputs/raw_data/FISHMORPH_Database.csv", sep = ";")
  #database <- database %>% dplyr::rename("Species" = "Genus.species")
  #database <- read.csv("inputs/raw_data/1873_2_RivFishTIME_SurveyTable.csv")
  
  # Required packages
  require(rfishbase)
  require(dplyr)
  require(tidyr)
  
  # Create list of all species occurring in database
  specieslist <- database %>% dplyr::summarise(Species = unique(Species))

  # Filter through fishbase for synonyms for database
  synonymslist <- specieslist %>% pull(Species) %>% synonyms() %>%
    dplyr::select(provided_name = synonym, valid_name = Species, Comment = Status) 
  
  # Replace tricky case issues
  #synonymslist$Comment[synonymslist$Comment == "Accepted name"] <- "accepted name"
  #for(i in 1:length(synonymslist)){
   # if((synonymslist$provided_name[i] != synonymslist$valid_name[i]) &
      # (synonymslist$Comment[i] == "accepted name")){
       # synonymslist$valid_name[i] = synonymslist$provided_name[i]
     # }
   # if((is.na(synonymslist$valid_name[i])) &
      # (synonymslist$Comment[i] == "accepted name")){
      #synonymslist$valid_name[i] <- synonymslist$provided_name[i]
      #}
  #}
    
  # Remove weird case where valid_name is "Genus Species"
  synonymslist <- synonymslist %>% filter(valid_name != "Genus Species")
  
  # Create df for the synonyms list which summarises counts for each comment category
  synonymslist.wide <- synonymslist %>% 
    dplyr::select(-valid_name) %>% 
    group_by(provided_name, Comment) %>%
    dplyr::summarise(count = n()) %>% 
    tidyr::pivot_wider(names_from = Comment, values_from = count)
  
  # Create list of species to take synonyms for and then create a clean synonyms list
  if("provisionally accepted name" %in% colnames(synonymslist.wide)){
    take_synonym_list <- synonymslist.wide %>% 
      dplyr::filter((is.na(`accepted name`)) & 
               is.na(`provisionally accepted name`) & 
               synonym == 1) %>% 
      dplyr::select(provided_name)
    synonymslist.clean <- synonymslist %>% 
      dplyr::filter((provided_name == valid_name & Comment %in% c("accepted name", "provisionally accepted name"))|
               (provided_name != valid_name & Comment %in% c("accepted name"))|
               (provided_name %in% take_synonym_list$provided_name & Comment == "synonym")|
               (is.na(valid_name) & is.na(Comment)))
  }else{
    take_synonym_list <- synonymslist.wide %>% 
      dplyr::filter((is.na(`accepted name`)) & 
               synonym == 1) %>% 
      dplyr::select(provided_name)
    synonymslist.clean <- synonymslist %>% 
      dplyr::filter((provided_name == valid_name & Comment %in% c("accepted name"))|
               (provided_name != valid_name & Comment %in% c("accepted name")) |
               (provided_name %in% take_synonym_list$provided_name & Comment == "synonym")|
               (is.na(valid_name) & is.na(Comment)))
  }

  # Note for fishmorph the database entries goes down by one species 
  # but this species doesn't occur in survey data
  
  # Add species reference ID to synonyms list df
  synonymslist.clean <- tibble::rowid_to_column(synonymslist.clean, "sp.ID")
  
  # Add species reference ID to original database
  database.ID <- database %>% group_by(Species) %>%
    mutate(sp.ID = match(unique(Species), synonymslist.clean$provided_name))
  
  # Merge dataframes
  database.new <- merge(database.ID, synonymslist.clean, by = "sp.ID")
  
  # Replace column "Species" with "Accepted Name"
  database.new <- database.new %>% dplyr::select(-Species, -provided_name, -Comment, -sp.ID)
  
  # Change valid name column to Species
  database.new <- database.new %>% dplyr::rename("Species" = "valid_name")
  
  # Return the new database with fishbase name harmonisation
  return (database.new)

}


