rm(list = ls()) #Clear working memory

# Required packages
my.packs <- c(
  
  #Accessing data
  'RPostgreSQL', 
  
  #Data manipulation
  'dplyr', 'reshape2', 'raster',
  
  #Data analysis
  'jagsUI',
  
  #Plotting
  'ggplot2', 'viridis', 'ggrepel',"RColorBrewer")

# if any of them are not installed, install them
if (any(!my.packs %in% installed.packages()[, 'Package'])) {install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}

# load packages
lapply(my.packs, require, character.only = TRUE)

####################################################################
# Study parameters
####################################################################
# Minimum number of years of counts per site to include in analysis
min_obs = 10

# Range of years to consider in analysis
min_year = 1984
max_year = 2017

################################################################
# PART 2: LOAD DOWNLOADED COUNT DATA AND FORMAT FOR ANALYSIS
# NOTE THAT PENGUIN COUNTS AND ASSOCIATED PRECISION SCORES WERE DOWNLOADED FROM
# MAPPPD USING CODE FROM CHE-CASTALDO ET AL. 2017 (SEE THEIR APPENDIX)
################################################################

# Read in data
siteloc = read.csv("../../data/SiteLocations.csv") # colony locations
load(file = "../generated_files/yframe2.rda")      # counts and associated precision from MAPPPD
load(file = "../generated_files/SiteList2.rda")    # names and ccamlr regions for each colony

# Merge sitelist info
yframe = merge(yframe, SiteList[,c("site_id","ccamlr_id")], by = "site_id")

# Number of years with data for each site
number_of_counts = aggregate(season~site_id, data = yframe, FUN = function(x) length(unique(x)))
number_of_counts = merge(number_of_counts, siteloc[,c("site_id","ccamlr_id","latitude","longitude")], all.x = TRUE)

# Only keep sites with enough data
sites_to_keep = unique(subset(number_of_counts, season >= min_obs)$site_id)
rows_to_keep = which(SiteList$site_id %in% sites_to_keep)
SiteList_keep = subset(SiteList, site_id %in% sites_to_keep)
sites_ordered = SiteList_keep$site_id
regions_ordered = SiteList_keep$ccamlr_id

# spatial information for each site
siteloc = subset(siteloc, site_id %in% sites_ordered)
siteloc$region_factor = as.factor(siteloc$ccamlr_id)

# Generate data matrices for JAGS
season = subset(yframe, site_id %in% sites_ordered)$season
n_sites = length(sites_to_keep)
n_seasons = length(min(season):max(season))
first_year = min(season)
year_seq = min(season):max(season)

nest_counts = matrix(NA,nrow = n_sites, ncol = n_seasons)
rownames(nest_counts) = sites_ordered
colnames(nest_counts) = year_seq
chick_counts = adult_counts = nest_precision = chick_precision = adult_precision = nest_counts

# Populate matrices
for (s in sites_ordered){
  site_dat = subset(yframe, site_id == s)
  for (y in year_seq){
    nests_y = subset(site_dat, count_type == "nests" & season == y)
    chicks_y = subset(site_dat, count_type == "chicks" & season == y)
    adults_y = subset(site_dat, count_type == "adults" & season == y)
    
    if (nrow(nests_y)>1) print ("More than 1 nest count for site",s,"in year",y)
    if (nrow(chicks_y)>1) print ("More than 1 chick count for site",s,"in year",y)
    if (nrow(adults_y)>1) print ("More than 1 adult count for site",s,"in year",y)
    
    if (nrow(nests_y)>0){
      nest_counts[which(sites_ordered == s),which(year_seq == y)] = nests_y$count
      nest_precision[which(sites_ordered == s),which(year_seq == y)] = nests_y$precision
    }
    
    if (nrow(chicks_y)>0){
      chick_counts[which(sites_ordered == s),which(year_seq == y)] = chicks_y$count
      chick_precision[which(sites_ordered == s),which(year_seq == y)] = chicks_y$precision
    }
    
    if (nrow(adults_y)>0){
      adult_counts[which(sites_ordered == s),which(year_seq == y)] = adults_y$count
      adult_precision[which(sites_ordered == s),which(year_seq == y)] = adults_y$precision
    }
    
  }
}

#Rearrange count matrices into dataframe (for plotting with ggplot)
nest_df = melt(nest_counts); chick_df = melt(chick_counts); adult_df = melt(adult_counts)
colnames(nest_df) = colnames(chick_df) = colnames(adult_df) = c("Site","Year","Count")
nest_df[,"Type"] = "Nests";chick_df[,"Type"] = "Chicks";adult_df[,"Type"] = "Adults"
count_df = rbind(nest_df, chick_df, adult_df)
count_df = merge(count_df, siteloc, by.x = "Site", by.y = "site_id", all.x = TRUE)

#Reorganize dataframe by region
count_df <- count_df[order(count_df$region_factor, count_df$Site, count_df$Year),] 
count_df$Site = factor(count_df$Site, levels = unique(count_df$Site))

countplot = ggplot() + 
  geom_point(data = count_df,aes(x = Year, y = Count, col = region_factor, shape = Type), alpha = 0.75)+
  theme_bw()+
  scale_color_manual(values = viridis(length(unique(count_df$region_factor)), option = "D", end = 0.8))+
  scale_shape_manual(values = c(1,22,19))+
  facet_grid(Site~., scales = "free_y")

ADPE_data = list(nest_counts = nest_counts,
                 chick_counts = chick_counts,
                 adult_counts = adult_counts,
                 
                 nest_precision = nest_precision,
                 chick_precision = chick_precision,
                 adult_precision = adult_precision,
                 
                 n_sites = n_sites,
                 n_seasons = n_seasons,
                 year_seq = year_seq,
                 sites_ordered = sites_ordered,
                 regions_ordered = regions_ordered,
                 
                 siteloc = siteloc,
                 count_df = count_df
)

# # Add data from Terre Adelie (first count is in 1984)
# NOTE: THIS DATA IS NOT IN MAPPPD, BUT IS AVAILABLE UPON REQUEST FROM C. BARBRAUD
#TA_data = read.csv("../../data/ADPE_counts_Terre_Adelie_Petrel_Island.csv")

# Create a blank placeholder
TA_data = data.frame(YEAR = seq(1952,2017),
                     SPECIES = "ADPE",
                     ABUN = NA,
                     CHICKS = NA,
                     BS = NA,
                     PHEN = NA)
TA_data = subset(TA_data, YEAR %in% ADPE_data$year_seq)

for (i in 1:6) ADPE_data[[i]] = as.data.frame(ADPE_data[[i]])
ADPE_data$nest_counts["PGEO",] = TA_data$ABUN
ADPE_data$chick_counts["PGEO",] = TA_data$CHICKS
ADPE_data$adult_counts["PGEO",] = NA # no adult counts are present at TA

#Assume data at TA is collected with maximum precision
ADPE_data$nest_precision["PGEO",] = max(ADPE_data$nest_precision,na.rm=TRUE)
ADPE_data$chick_precision["PGEO",] = max(ADPE_data$nest_precision,na.rm=TRUE)
ADPE_data$adult_precision["PGEO",] = NA

ADPE_data$siteloc$site_id = as.character(ADPE_data$siteloc$site_id)
ADPE_data$siteloc$site_name = as.character(ADPE_data$siteloc$site_name)
ADPE_data$siteloc$region = as.character(ADPE_data$siteloc$region)

ADPE_data$siteloc["PGEO",] = c(site_id = "PGEO",
                               site_name = "Point Geologie",
                               region = "East Antarctica",
                               ccamlr_id = "58.4.1",
                               latitude = -66.67,
                               latitude = 140.01,
                               input_srid = 4326,
                               region_factor = "58.4.1")

ADPE_data$siteloc$site_id = as.factor(ADPE_data$siteloc$site_id)
ADPE_data$siteloc$site_name = as.factor(ADPE_data$siteloc$site_name)
ADPE_data$siteloc$region = as.factor(ADPE_data$siteloc$region)
ADPE_data$siteloc$latitude = as.numeric(ADPE_data$siteloc$latitude)
ADPE_data$siteloc$longitude = as.numeric(ADPE_data$siteloc$longitude)

ADPE_data$n_sites = nrow(ADPE_data$siteloc)
ADPE_data$n_seasons = ncol(ADPE_data$nest_counts)
ADPE_data$year_seq = ADPE_data$year_seq
ADPE_data$sites_ordered = as.factor(as.character(ADPE_data$siteloc$site_id))
ADPE_data$regions_ordered = as.factor(as.character(ADPE_data$siteloc$region_factor))

# Map of Study
world <- map_data("world") 
worldmap <- ggplot(world, aes(x=long, y=lat, group=group)) + 
  geom_polygon(aes(fill = group), col = "gray50") + 
  geom_point(data = ADPE_data$siteloc, aes(y = latitude, x = longitude, group = 1, col = latitude), alpha = 1, size = 2)+
  geom_label_repel(data = ADPE_data$siteloc, aes(y = latitude, x = longitude, group = 1, label = site_id, col = latitude), alpha = 1, size = 3)+
  
  theme_bw()+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  
  scale_color_gradientn(colors = viridis(10)[1:9], guide = F)+
  scale_fill_gradientn(colors = c("gray90"), guide = F)+
  
  scale_y_continuous(breaks=(-2:2) * 30, limits = c(-90,-60)) + 
  scale_x_continuous(breaks=(-4:4) * 45)
worldmap = worldmap + coord_map("ortho", orientation=c(-90, 0, 0)) 

#pdf("../figures/worldmap.pdf", width = 8, height = 8)
print(worldmap)
#dev.off()

print(length(ADPE_data$sites_ordered))

# Export site data
#write.csv(ADPE_data$siteloc, file = paste0("../generated_files/siteloc_analysis_",min_obs,"_yrs.csv"), row.names = FALSE)
#save(ADPE_data, file = paste0("../generated_files/ADPE_data_",min_obs,"_yrs.rda"))

################################################################
# PART 2: GENERATE SEA ICE COVARIATES
################################################################

sites = ADPE_data$sites_ordered #unique sites for use in Adelie model
regions = ADPE_data$regions_ordered
years = ADPE_data$year_seq #should match the range specified at the top of script

# # Monthly SIC at each site in 500 km radius around breeding colony
# sic500 = read.csv("../../data/ADPE_SIC_500.csv")
# 
# sic500_melt = melt(sic500, id = c("site_id","year")); sic500_melt$month = as.numeric(sic500_melt$variable)
# 
# covar_df = data.frame()
# for (s in sites){
#   start = Sys.time()
#   for (y in years){
# 
#     #Months used in analysis
#     relevant_mos = c(3,4,5,6)
# 
#     #1 winter prior to current breeding season (0 years prior)
#     temp = subset(sic500_melt, site_id == s & year %in% y & month %in% relevant_mos)
#     SIC_1 = NA
#     if (nrow(temp) == length(relevant_mos)) SIC_1 = mean(temp$value)
# 
#     #2 winters prior...
#     temp = subset(sic500_melt, site_id == s & year %in% (y-1) & month %in% relevant_mos)
#     SIC_2 = NA
#     if (nrow(temp) == length(relevant_mos)) SIC_2 = mean(temp$value)
# 
#     #3
#     temp = subset(sic500_melt, site_id == s & year %in% (y-2) & month %in% relevant_mos)
#     SIC_3 = NA
#     if (nrow(temp) == length(relevant_mos)) SIC_3 = mean(temp$value)
# 
#     #4
#     temp = subset(sic500_melt, site_id == s & year %in% (y-3) & month %in% relevant_mos)
#     SIC_4 = NA
#     if (nrow(temp) == length(relevant_mos)) SIC_4 = mean(temp$value)
# 
#     #5
#     temp = subset(sic500_melt, site_id == s & year %in% (y-4) & month %in% relevant_mos)
#     SIC_5 = NA
#     if (nrow(temp) == length(relevant_mos)) SIC_5 = mean(temp$value)
# 
#     #6 winters prior...
#     temp = subset(sic500_melt, site_id == s & year %in% (y-5) & month %in% relevant_mos)
#     SIC_6 = NA
#     if (nrow(temp) == length(relevant_mos)) SIC_6 = mean(temp$value)
# 
#     #########################################################################
#     #Add variables to dataframe
#     #########################################################################
#     row_dat = data.frame(site_id = s,
# 
#                          year = y,
# 
#                          region = regions[which(ADPE_data$sites_ordered == s)],
# 
#                          SIC_1 = SIC_1,
#                          SIC_2 = SIC_2,
#                          SIC_3 = SIC_3,
#                          SIC_4 = SIC_4,
#                          SIC_5 = SIC_5,
#                          SIC_6 = SIC_6
#     )
# 
#     covar_df = rbind(covar_df, row_dat)
#   }
#   stop = Sys.time()
#   print(paste(s, round(stop-start,3), "seconds"))
# }

#save covariates
#save(covar_df, file = paste0("../generated_files/ADPE_covariates_",min_obs,"_yrs.rda"))

################################################################
# PART 3: PREPARE DATA FOR BAYESIAN ANALYSIS
################################################################

load(file = paste0("../generated_files/ADPE_covariates_",min_obs,"_yrs.rda"))

# Prepare data for analysis
nest_train = as.matrix(ADPE_data$nest_counts)
chick_train = as.matrix(ADPE_data$chick_counts)
adult_train = as.matrix(ADPE_data$adult_counts)

# Long format data (for use in observation model)
nest_df = melt(nest_train); colnames(nest_df) = c("Site","Year","Nests")
chick_df = melt(chick_train); colnames(chick_df) = c("Site","Year","Chicks")
adult_df = melt(adult_train); colnames(adult_df) = c("Site","Year","Adults")
count_df = cbind(nest_df,chick_df$Chicks, adult_df$Adults)
count_df_full = count_df
count_df_full[is.na(count_df_full)] = -1
count_df_full = count_df_full[-which(apply(count_df_full[,3:5],1,sum) == -3),] #Remove rows with NAs for nests, chicks, and adults
count_df_full = count_df_full[order(count_df_full$Site,count_df_full$Year),]

# Fill in precision for counts with NA values (this doesn't affect model estimates)
ADPE_data$nest_precision[is.na(ADPE_data$nest_precision)]   = 1612.798
ADPE_data$chick_precision[is.na(ADPE_data$chick_precision)] = 1612.798
ADPE_data$adult_precision[is.na(ADPE_data$adult_precision)] = 1612.798

ADPE_data$nest_precision = as.matrix(ADPE_data$nest_precision)
ADPE_data$chick_precision = as.matrix(ADPE_data$chick_precision)
ADPE_data$adult_precision = as.matrix(ADPE_data$adult_precision)

# Package data to send to JAGS
jags.data = list()

jags.data$site = as.numeric(count_df_full$Site)
jags.data$n_seasons = ADPE_data$n_seasons
jags.data$n_sites = ADPE_data$n_sites
jags.data$season = as.numeric(count_df_full$Year) - min(as.numeric(count_df_full$Year)) + 1
jags.data$n_length = length(jags.data$site)

# Number of sites in each region
jags.data$nest_train = nest_train
jags.data$chick_train = chick_train
jags.data$adult_train = adult_train

jags.data$nest_precision = ADPE_data$nest_precision
jags.data$chick_precision = ADPE_data$chick_precision
jags.data$adult_precision = ADPE_data$adult_precision

jags.data$first_obs = aggregate(Year~Site, data = count_df_full, FUN = min)$Year - min(ADPE_data$year_seq) + 1
jags.data$last_obs = aggregate(Year~Site, data = count_df_full, FUN = max)$Year - min(ADPE_data$year_seq) + 1

#Format SIC covariates for jags analysis
#Function to create appropriate covariate matrices for each region
covar_mat_fn = function(df, variable_name){
  
  df = df[,c("year","region","site_id",variable_name)]
  
  covar_mat = matrix(NA, nrow = length(unique(df$site_id)), ncol = length(ADPE_data$year_seq))
  colnames(covar_mat) = ADPE_data$year_seq
  rownames(covar_mat) = unique(df$site_id)
  
  for (s in 1:length(unique(df$site_id))){
    for (y in 1:length(ADPE_data$year_seq)){
      covar_mat[s,y] = df[which(df$year == ADPE_data$year_seq[y] & df$site_id == unique(df$site_id)[s]),4]
    }
  }
  
  # Standardize each row (i.e., site) separately
  covar_mat_zstand = covar_mat*NA
  for (row in 1:nrow(covar_mat)){
    covar_mat_zstand[row,] = (covar_mat[row,] - mean(covar_mat[row,]))/sd(covar_mat[row,])
  }
  
  # Or: standardize entire matrix
  #covar_mat_zstand = (covar_mat - mean(covar_mat))/sd(covar_mat)
  
  return(list(cov_mat = covar_mat,
              cov_mat_zstand = covar_mat_zstand))
}

# Covariates
# mean summer sea ice over previous seasons
SIC_1yr = covar_mat_fn(df = covar_df,variable_name = "SIC_1")
SIC_2yr = covar_mat_fn(df = covar_df,variable_name = "SIC_2")
SIC_3yr = covar_mat_fn(df = covar_df,variable_name = "SIC_3")
SIC_4yr = covar_mat_fn(df = covar_df,variable_name = "SIC_4")
SIC_5yr = covar_mat_fn(df = covar_df,variable_name = "SIC_5")
SIC_6yr = covar_mat_fn(df = covar_df,variable_name = "SIC_6")

#long-term mean SIC at site
SIC_sitemean = apply(SIC_1yr$cov_mat,1,mean)

site_order = names(sort(SIC_sitemean)) #sites ordered by mean SIC

jags.data$SIC_anomaly_1yr = SIC_1yr$cov_mat_zstand
jags.data$SIC_anomaly_2yr = SIC_2yr$cov_mat_zstand
jags.data$SIC_anomaly_3yr = SIC_3yr$cov_mat_zstand
jags.data$SIC_anomaly_4yr = SIC_4yr$cov_mat_zstand
jags.data$SIC_anomaly_5yr = SIC_5yr$cov_mat_zstand
jags.data$SIC_anomaly_6yr = SIC_6yr$cov_mat_zstand

jags.data$SIC_sitemean = (SIC_sitemean - mean(SIC_sitemean))/sd(SIC_sitemean)
jags.data$SIC_sitemean_unscaled = SIC_sitemean

# Since 0 is not supported in logit transform
jags.data$chick_train[which(jags.data$chick_train == 0)] = 1
jags.data$region = as.numeric(ADPE_data$regions_ordered)
jags.data$n_regions = length(unique(jags.data$region))

save(jags.data, file = paste0("../generated_files/jags_data_multilevel_SIC_",min_obs,"_yrs.rda"))
