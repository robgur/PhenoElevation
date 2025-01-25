library(lubridate)
library(phenesse)
library(terra)
library(prism)
library(raster)
library(readr)
library(car)
library(performance)
library(lme4)
library(lmerTest)
library(MuMIn)
library(splines)
library(rtrees)
library(phyr)
library(stringr)
library(brms)
library(curl)
library(sf)

######################### load up some data files
SpeciesTraits2 <- read_csv("C:/Users/robgu/Downloads/SpeciesTraits2.csv")
Quantile_all449 <- read_csv("C:/Users/robgu/Downloads/Quantile_all449.csv")

######################### get phylogeny
sp_list <- subset(SpeciesTraits2, select=c('species','genus','family'))
sp_list$species <- str_to_title(sp_list$species)
plant_tree1 = get_tree(sp_list = sp_list, taxon = "plant", scenario = "at_basal_node", show_grafted = TRUE)

######################## get_climate
test <- rnaturalearth::ne_states(country = c("United States of America"),returnclass = "sf")
test2 <-st_union(test, by_feature = FALSE)

Quantile_all449$origin <- as.Date(paste0(Quantile_all449$year, "-01-01"),tz = "UTC") - days(1)
Quantile_all449$onset_day <- as.integer(Quantile_all449$onset)
Quantile_all449$onset_date<- as.Date(Quantile_all449$onset_day, origin = Quantile_all449$origin, tz = "UTC") 
Quantile_all449$start_date <- Quantile_all449$onset_date - 90


Quantile_all449_sf <- Quantile_all449  %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = "+proj=longlat +datum=WGS84 +no_defs")

Quantile_all4492 <- st_intersection(Quantile_all449_sf,test)

Quantile_all4492_2 <- Quantile_all4492  %>%
  dplyr::mutate(longitude = sf::st_coordinates(.)[,1],
                latitude = sf::st_coordinates(.)[,2])

Quantile_all4492_3 <- Quantile_all4492_2[,c(1:49,84:85)]

#####Accumulate_climate_for_onset
clim_slope <- data.frame()
clim_resid <-data.frame()
clim_mean <- data.frame() 
clim_precip <- data.frame()

for (i in 1:nrow(Quantile_all4492_3)){
  dm_url <- paste("https://daymet.ornl.gov/single-pixel/api/data?lat=",Quantile_all4492_3$latitude[i],"&lon=",Quantile_all4492_3$longitude[i],"&vars=prcp,tmax,tmin&start=",Quantile_all4492_3$start_date[i],"&end=",Quantile_all4492_3$onset_date[i])
  dm_url2 <- gsub(" ", "", dm_url)
  test2 <- curl(dm_url2)
  out <- readLines(test2)
  out2 <- out[7:length(out)]
  out3 <- read.table(text=paste(out2, collapse='\n'), header = TRUE, stringsAsFactors = FALSE, sep=',')
  out3$tmean <- (out3$tmax..deg.c. + out3$tmin..deg.c.)/2
  tot_mean <- mean(out3$tmean)
  model <- lm(tmean ~ yday, data = out3)
  slope <- model$coefficients[2]
  resid <- sum(resid(model)^2)
  precip <- mean(out3$prcp..mm.day.)
  clim_mean = rbind(clim_mean, tot_mean)
  clim_resid = rbind(clim_resid, resid)
  clim_slope = rbind(clim_slope, slope)
  clim_precip = rbind(clim_precip,precip)
}

Quantile_all4492_5 <- Quantile_all4492_3 %>%  st_drop_geometry()
Quantile_all4492_5 <- cbind(Quantile_all4492_5,clim_mean,clim_resid,clim_slope,clim_precip)
Quantile_all4492_5 <- Quantile_all4492_5 %>%
  rename("clim_mean_onset" = 52) %>%
  rename("clim_resid_onset" = 53) %>%
  rename("clim_slope_onset" =54) %>%
  rename("clim_precip_onset" =55)

##### accumulate_climate_offset

Quantile_all4492_5$mean_day <- as.integer(Quantile_all4492_5$mean)
Quantile_all4492_5$mean_date<- as.Date(Quantile_all4492_5$mean_day, origin = Quantile_all4492_5$origin, tz = "UTC") 
Quantile_all4492_5$end_date <- Quantile_all4492_5$mean_date + 90

clim_slope_off <- data.frame()
clim_resid_off <-data.frame()
clim_mean_off <- data.frame() 
clim_precip_off <- data.frame()

for (i in 1:nrow(Quantile_all4492_5)){
  dm_url <- paste("https://daymet.ornl.gov/single-pixel/api/data?lat=",Quantile_all4492_5$latitude[i],"&lon=",Quantile_all4492_5$longitude[i],"&vars=prcp,tmax,tmin&start=",Quantile_all4492_5$mean_date[i],"&end=",Quantile_all4492_5$end_date[i])
  dm_url2 <- gsub(" ", "", dm_url)
  test2 <- curl(dm_url2)
  out <- readLines(test2)
  out2 <- out[7:length(out)]
  out3 <- read.table(text=paste(out2, collapse='\n'), header = TRUE, stringsAsFactors = FALSE, sep=',')
  out3$tmean <- (out3$tmax..deg.c. + out3$tmin..deg.c.)/2
  tot_mean <- mean(out3$tmean)
  model <- lm(tmean ~ yday, data = out3)
  slope = model$coefficients[2]
  resid <- sum(resid(model)^2)
  precip <- mean(out3$prcp..mm.day.)
  clim_mean_off = rbind(clim_mean_off, tot_mean)
  clim_resid_off = rbind(clim_resid_off, resid)
  clim_slope_off = rbind(clim_slope_off, slope)
  clim_precip_off = rbind(clim_precip_off,precip)
}

Quantile_all4492_5 <- cbind(Quantile_all4492_5,clim_mean_off,clim_resid_off,clim_slope_off,clim_precip_off)
Quantile_all4492_5 <- Quantile_all4492_5 %>%
  rename("clim_mean_offset" = 59) %>%
  rename("clim_slope_offset" = 60) %>%
  rename("clim_resid_offset" =61) %>%
  rename("clim_precip_offset" =62)

####climate_context
Quantile_all4492_5$offset_day <- as.integer(Quantile_all4492_5$offset)
Quantile_all4492_5$offset_date<- as.Date(Quantile_all4492_5$offset_day, origin = Quantile_all4492_5$origin, tz = "UTC") 

clim_mean_dur <- data.frame()
clim_precip_dur <- data.frame()

for (i in 1:nrow(Quantile_all4492_5)){
  dm_url <- paste("https://daymet.ornl.gov/single-pixel/api/data?lat=",Quantile_all4492_5$latitude[i],"&lon=",Quantile_all4492_5$longitude[i],"&vars=prcp,tmax,tmin&start=",Quantile_all4492_5$start_date[i],"&end=",Quantile_all4492_5$offset_date[i])
  dm_url2 <- gsub(" ", "", dm_url)
  test2 <- curl(dm_url2)
  out <- readLines(test2)
  out2 <- out[7:length(out)]
  out3 <- read.table(text=paste(out2, collapse='\n'), header = TRUE, stringsAsFactors = FALSE, sep=',')
  out3$tmean <- (out3$tmax..deg.c. + out3$tmin..deg.c.)/2
  tot_mean <- mean(out3$tmean)
  precip <- mean(out3$prcp..mm.day.)
  clim_mean_dur = rbind(clim_mean_dur, tot_mean)
  clim_precip_dur = rbind(clim_precip_dur, precip) 
}

Quantile_all4492_5 <- cbind(Quantile_all4492_5,clim_mean_dur,clim_precip_dur)
Quantile_all4492_5 <- Quantile_all4492_5 %>%
  rename("clim_mean_duration" = 65) %>%
  rename("clim_precip_duration" =66)

############################################################
# Phy_r version



#############################################################
#brms version

Quantile_all4492_5$species <- str_to_title(Quantile_all4492_5$species)
plant_tree1$tip.label <- str_remove(plant_tree1$tip.label, pattern = "\\*")
plant_tree1$tip.label <- str_remove(plant_tree1$tip.label, pattern = "\\*")
A <- ape::vcv.phylo(plant_tree1)
duration_test <- brm(duration ~ (ns(elevation_sc,2)+ns(latitude_sc,2)+bio12_sc)^2  + bio4_sc + obsnum_sc + (1 + elevation_sc + latitude_sc | gr(species, cov = A)), data=Quantile_all449, data2 = list(A = A), family=Gamma(link="log"), chains = 4,cores = 4, iter=10000)
#duration_test <- brm(duration ~ ts(elevation_sc)*ts(latitude_sc)*bio13_sc  + bio4_sc + obsnum_sc + (1 + elevation_sc + latitude_sc | gr(species, cov = A)), data=Quantile_all449, data2 = list(A = A), family=Gamma(link="log"), chains = 4,cores = 4, iter=10000)

#############################################
#check_model_results_including_fit

plot(duration_test, N=2, ask=FALSE)
summary(duration_test)  #note that Rhat is 1, which is good - convergence
plot(duration_test, N=2, ask=FALSE)  #plots look good
pp_check(duration_test)  #very helful - model fits the data quite well
plot_model(duration_test, type = "pred", terms = c("elevation_sc","latitude_sc","bio12_sc"))
#plot_model works but uncertainty of error terms are not taken into account. We may want to use `rstantools::posterior_predict()`.

######################################
#check for phylogenetic_signal

hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + shape^2) = 0"
hyp <- hypothesis(duration_test,hyp,class=NULL)
plot(hyp)
#this suggests there is almost no phylogenetic signal in this model. We can likely simpify and run a model without phylogenetic autocorrelation term
