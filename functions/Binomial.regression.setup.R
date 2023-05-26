# This function sets up the binomial linear regressiona models for looking at
# co-occurence of different novelty categories between func and tax

######################### STEP 1: SETUP MODEL DF ###########################

# For now I am going to "pretend" my tax results are functional results also 
#(i.e. they should always co-occur)
#func.results.TEMP.df <- pilot.tax.results.df
#tax.results.df <- pilot.tax.results.df

# Rename cat column to func and tax cat
#func.results.TEMP.df <- func.results.TEMP.df %>% rename(func.cat = cat)
#tax.results.df <- tax.results.df %>% rename(tax.cat = cat)

# Merge two novelty results dataframes
##model.df <- merge(tax.results.df, func.results.TEMP.df, by = c("site", "bins", "bin.lag")) 

# Modify dataframe to only include information relevant for binomial regression
#model.df <- model.df %>%
 # group_by(site) %>% dplyr::summarise(site=site, bins, bin.lag,
                #  TaxNov=ifelse(tax.cat=="novel", 1, 0),
                #  TaxInst=ifelse(tax.cat=="instant", 1, 0),
                #  TaxCumul=ifelse(tax.cat=="cumul", 1, 0),
                 # FuncNov=ifelse(func.cat=="novel", 1, 0),
                  #FuncInst=ifelse(func.cat=="instant", 1, 0),
                  #FuncCumul=ifelse(func.cat=="cumul", 1, 0))

# Convert into table of co-occurences between categories
#model.frequency.df <- model.df %>%
 # group_by(site) %>%
 # dplyr::summarise()
                  
# Create binomial regression models for each condition
# Need a way to make this simpler using tapply or lapply or other

# Function to create a list of glmm outputs based on model.df
#create.glmm.df <- function(model.df)
  
  # Create empty list for glmm results
  #glmm.list <- ()

  # 
  

  
  
#glmer(model.df$TaxNov ~ model.df$FuncNov + model.df$FuncInst +
#                      model.df$FuncCumul + (1|model.df$site), family = binomial)
#taxInst.glmm <- glmer(model.df$TaxInst ~ model.df$FuncNov + model.df$FuncInst +
#                        model.df$FuncCumul + (1|model.df$site), family = binomial)
#taxCumul.glmm <- glmer(model.df$TaxCumul ~ model.df$FuncNov + model.df$FuncInst +
#                        model.df$FuncCumul + (1|model.df$site), family = binomial)
#funcNov.glmm <- glmer(model.df$FuncNov ~ model.df$TaxNov + model.df$TaxInst +
#                        model.df$TaxCumul + (1|model.df$site), family = binomial)
#funcInst.glmm <- glmer(model.df$FuncInst ~ model.df$TaxNov + model.df$TaxInst +
#                        model.df$TaxCumul + (1|model.df$site), family = binomial)
#funcCumul.glmm <- glmer(model.df$FuncCumul ~ model.df$TaxNov + model.df$TaxInst +
#                        model.df$TaxCumul + (1|model.df$site), family = binomial)

# Create plots of 
#plot(model.df$TaxNov, model.df$FuncNov)


#func_prob_results <- lapply(func_comm_state, function(x){
 # glmer(cbind(success, failure) ~ 1 + (1|site), #+ (1|season)
        #data=x, family=binomial)

# For now I am going to "pretend" my tax results are functional results also 
#(i.e. they should always co-occur)



# I need to to two analyses: 
#1. conditional probability of tax nov if func nov
#2. condition probability of func nov if tax nov
  
# I will try this first on an example result for func and tax nov
# Dummy data is 5 time series

##### Pseudo code below ####

# Create 

# Tax novelty results
#tax_nov_results

# In these two df the novelty categories are stored as T/F 
# So I need to coerce these into a df which looks like this

# 



  

  
  
  
  
  
  
  
  