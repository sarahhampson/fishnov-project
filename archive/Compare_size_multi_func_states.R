# Just a random script to look at comparing func novelty for size only vs multi trait approach

# First import results for both novelty detection methods
multi_results <- full_results_df # Already in environment
size_results <- read.csv("./archive/size_results_df.csv")

# Subset only sqID, bins, taxnovel and funcnovel
multi_results_tax <- multi_results %>% dplyr::select(sqID, bins, taxnovel) %>% 
  filter(taxnovel==TRUE) %>% mutate(ID = paste0(sqID, " ", bins))
size_results_tax <- size_results %>% dplyr::select(sqID, bins, taxnovel) %>% 
  filter(taxnovel==TRUE) %>% mutate(ID = paste0(sqID, " ", bins))
multi_results_func <- multi_results %>% dplyr::select(sqID, bins, funcnovel) %>% 
  filter(funcnovel==TRUE) %>% mutate(ID = paste0(sqID, " ", bins))
size_results_func <- size_results %>% dplyr::select(sqID, bins, funcnovel) %>% 
  filter(funcnovel==TRUE) %>% mutate(ID = paste0(sqID, " ", bins))

# Find intersect
common_sites_tax <- as.data.frame(intersect(multi_results_tax$ID, size_results_tax$ID)) # 251 --> they match woo!
common_sites_func <- as.data.frame(intersect(multi_results_func$ID, size_results_func$ID)) # 37 --> not all match

# Find IDs which are now novel which weren't before
diff_sites_multi_func <- multi_results_func %>% filter(!ID %in% common_sites_func[,1])

# Save these in archive
write.csv(diff_sites_multi_func, "./archive/diff_sites_multi_func.csv")

# Cleanup
rm(multi_results, size_results, multi_results_func, size_results_func, 
   multi_results_tax, size_results_tax, common_sites_func, common_sites_tax,
   diff_sites_multi_func)
  