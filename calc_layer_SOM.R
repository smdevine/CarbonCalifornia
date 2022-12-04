#updated to include depth-weighted SOM for soils < 100 cm deep
#TO-DO: (1) add soil depth calc to map-units (2) convert mukey to character class
#left off at line 273, waiting for component restriction table information
#functions to work with ProfileApply
#first define top and bottom of layer of interest
bottom_depth <- 50
depth_top <- 30
layer_thickness <- bottom_depth - depth_top
subDir_name <- paste0(depth_top, '_', bottom_depth, 'cm_aggregation')
fname_depth <- paste0('_', depth_top, '_', bottom_depth, 'cm_')

wtd.mean <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzdepb_r - x$hzdept_r
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=TRUE) #modified to TRUE on 2/16/22
  m
}
kgOrgC_sum <- function(x, slice_it=FALSE, depth, rm.NAs=TRUE, om_to_c=1.72) { #removing NAs so as to capture SOC conent in soils <100 cm deep; but true NAs must be dealt with below based on comparison of "soil depth" and depth to which SOM data exists
  if (slice_it) { 
    x <- horizons(x)[1:depth, ]
    depths(x) <- cokey ~ hzdept_r + hzdepb_r
  }
  thick <- x$hzdepb_r - x$hzdept_r
  sum((thick / 10) * (x$om_r / om_to_c) * x$dbthirdbar_r * (1 - x$fragvol_r_sum / 100), na.rm = rm.NAs)
} #don't need slice_it, because that taken care of within horizon_to_comp function

data_depth <- function(x, varname) {
  thick <- x$hzdepb_r - x$hzdept_r
  sum(thick[!is.na(x[[varname]])])
}
horizon_to_comp_layer <- function(horizon_SPC, depth, comp_df, vars_of_interest = 'om_r', varnames = 'SOM', top_depth=0, layername) {
  columnames <- paste0(varnames, '_', layername, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  assign("top_depth", top_depth, envir = .GlobalEnv)
  sliced_SPC <- slice(horizon_SPC, top_depth:(depth-1) ~ .) #depth was '0:depth' in previous version
  stopifnot(unique(sliced_SPC$pedon_key)==site(sliced_SPC)$pedon_key)
  for (i in seq_along(vars_of_interest)) {
    s <- site(sliced_SPC)
    s[[columnames[i]]] <- profileApply(sliced_SPC, FUN = wtd.mean, y=vars_of_interest[i])
    site(sliced_SPC) <- s
  }
  s <- site(sliced_SPC)
  s[[paste0('kgOrg.m2_', layername, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum)
  s[['SOM_depth']] <- profileApply(sliced_SPC, FUN=data_depth, varname='om_r')
  s[['BD_depth']] <- profileApply(sliced_SPC, FUN=data_depth, varname='dbthirdbar_r')
  s[['Rock_Frag_depth']] <- profileApply(sliced_SPC, FUN = data_depth, varname='fragvol_r_sum')
  columnames <- c(columnames, paste0('kgOrg.m2_', layername, 'cm'), 'SOM_depth', 'BD_depth', 'Rock_Frag_depth')
  rm(depth, envir = .GlobalEnv) #because we had to put it there earlier
  rm(top_depth, envir = .GlobalEnv)#because we had to put it there earlier
  s$compname <- comp_df$compname[match(s$cokey, comp_df$cokey)]
  s$mukey <- comp_df$mukey[match(s$cokey, comp_df$cokey)]
  s$comppct <- comp_df$comppct_r[match(s$cokey, comp_df$cokey)]
  s <- s[,c('mukey', 'cokey', 'compname', 'comppct', columnames)]
  s
}

MUaggregate <- function(df1, varname) {
  sapply(split(x=df1, f=df1$mukey), FUN=function(x) {if(sum(!is.na(x[[varname]]))==0) {NA} 
    else{sum(x$comppct[!is.na(x[[varname]])] * x[[varname]][!is.na(x[[varname]])] / sum(x$comppct[!is.na(x[[varname]])]))}
  })
}
MUAggregate_wrapper <- function(df1, varnames) {
  x <- sapply(varnames, FUN=MUaggregate, df1=df1)
  as.data.frame(cbind(mukey=as.integer(row.names(x)), x))
}
library(aqp)
mainDir <- 'C:/Users/smdevine/Desktop/post doc/soil organic carbon'
SSURGOdir <- file.path(mainDir, 'SSURGO data')
list.files(SSURGOdir)
rockNA_to_0 <- TRUE #to convert all NA for rock frag content to 0
add_missing_horizons <- FALSE

#read-in mapunit area info
mu_data <- read.csv(file.path(SSURGOdir, 'ca_mapunit_area.csv'), stringsAsFactors = FALSE)
mu_data$hectares <- mu_data$sqm/10000
sum(mu_data$hectares) #41,477,744 hectares

#read-in component data and define major component as >15% aerial coverage
list.files(SSURGOdir)
comp_data <- read.csv(file.path(SSURGOdir, 'ca_components.csv'), na.strings = c('', ' '), stringsAsFactors = FALSE)
lapply(comp_data, summary)
length(unique(comp_data$cokey)) #91658 unique cokeys
length(unique(comp_data$mukey)) #19036 unique mukeys
if(sum(is.na(comp_data$majcompflag)) > 0) {stop(print('there are NAs in majcomp column!'))}
if(sum(is.na(comp_data$comppct_r)) > 0) {stop(print('there are NAs in the comppct column!'))}
sum(comp_data$comppct_r[comp_data$majcompflag=='Yes'] < 15) #159 are major components with <15% aereal coverage
sum(comp_data$comppct_r[comp_data$majcompflag=='No '] >= 15) #284 are minor components with >15% aereal coverage
comp_data$majcompflag[comp_data$majcompflag=='No ' & comp_data$comppct_r>=15] <- 'Yes'
comp_data$majcompflag[comp_data$majcompflag=='Yes' & comp_data$comppct_r < 15] <- 'No '

#read-in soil horizon data and convert all 0% SOM reporting to NA
horizon_data <- read.csv(file.path(SSURGOdir, 'ca_horizons.csv'), na.strings = c('', ' '), stringsAsFactors = FALSE)
colnames(horizon_data)
length(unique(horizon_data$cokey)) #42115 unique cokeys
sum(horizon_data$om_r==0, na.rm = TRUE) #4408 instances of 0% SOM
hznames_zeroOM <- unique(horizon_data$hzname[which(horizon_data$om_r==0)])
hznames_zeroOM
horizon_data$om_r[horizon_data$om_r==0] <- NA
horizon_data$majcompflag <- comp_data$majcompflag[match(horizon_data$cokey, comp_data$cokey)]
table(horizon_data$majcompflag)
horizon_data_majcomps <- horizon_data[horizon_data$majcompflag=='Yes',]
horizon_data_majcomps$mukey <- comp_data$mukey[match(horizon_data_majcomps$cokey, comp_data$cokey)]
if(add_missing_horizons)  {
  horizon_data_majcomps <- rbind(horizon_data_majcomps, data.frame(cokey=as.integer(21156280), hzname=NA, hzdept_r=61, hzdepb_r=94, om_r=NA, dbthirdbar_r=NA, fragvol_r_sum=NA, majcompflag='Yes', mukey=as.integer(3251996))) #this derived from below analysis finding missing horizons from various cokeys to add missing horizons before the analylsis
  horizon_data_majcomps <- rbind(horizon_data_majcomps, data.frame(cokey=as.integer(22009763), hzname='H2', hzdept_r=15, hzdepb_r=41, om_r=NA, dbthirdbar_r=NA, fragvol_r_sum=NA, majcompflag='Yes', mukey=as.integer(471569)))
}
lapply(horizon_data_majcomps, class)
# horizon_data_majcomps <- subset(horizon_data_majcomps, horizon_data_majcomps$cokey!=22009763) #see below for rationale

#inspect bulk density data
summary(horizon_data_majcomps$dbthirdbar_r)
sum(horizon_data_majcomps$dbthirdbar_r < 0.5, na.rm = TRUE) #2886 have BD < 0.5 g cm^3

#inspect rock fragment data
summary(horizon_data_majcomps$fragvol_r_sum)
sum(horizon_data_majcomps$fragvol_r_sum > 100, na.rm = TRUE) #12 are greater than 0
horizon_data_majcomps[which(horizon_data_majcomps$fragvol_r_sum > 100),] #12 are 
sum(horizon_data_majcomps$fragvol_r_sum == 100, na.rm = TRUE) #7 equal to zero
sum(horizon_data_majcomps$fragvol_r_sum==0, na.rm = TRUE) #4704
sum(is.na(horizon_data_majcomps$fragvol_r_sum)) #28692
sum(is.na(horizon_data_majcomps$fragvol_r_sum) & !is.na(horizon_data_majcomps$om_r)) #9602
sum(is.na(horizon_data_majcomps$dbthirdbar_r) & !is.na(horizon_data_majcomps$om_r)) #366

#convert rock fragment NA to 0% rock fragment [don't do for now]
if(rockNA_to_0) {
  horizon_data_majcomps$fragvol_r_sum[is.na(horizon_data_majcomps$fragvol_r_sum)] <- 0
}

#then convert rock fragment data >= 100% to NA
horizon_data_majcomps$fragvol_r_sum[horizon_data_majcomps$fragvol_r_sum >= 100] <- NA
summary(horizon_data_majcomps$fragvol_r_sum)

#inspect some duripans
duripans <- horizon_data_majcomps[which(horizon_data_majcomps$hzname=='Bqm'),]
summary(duripans$om_r)
hist(duripans$hzdept_r)
sum(duripans$hzdept_r < 30) #10
duripans[duripans$hzdept_r < 30,]

#inspect O horizons
unique(horizon_data_majcomps$hzname[grepl('Oa', horizon_data_majcomps$hzname)])
length(horizon_data_majcomps$hzname[grepl('Oa', horizon_data_majcomps$hzname)])
unique(horizon_data_majcomps$hzname[grepl('Oi', horizon_data_majcomps$hzname)])
length(horizon_data_majcomps$hzname[grepl('Oi', horizon_data_majcomps$hzname)])
unique(horizon_data_majcomps$hzname[grepl('Oe', horizon_data_majcomps$hzname)])
length(horizon_data_majcomps$hzname[grepl('Oe', horizon_data_majcomps$hzname)])
summary(horizon_data_majcomps$dbthirdbar_r[grepl('Oa', horizon_data_majcomps$hzname)])
summary(horizon_data_majcomps$dbthirdbar_r[grepl('Oi', horizon_data_majcomps$hzname)])
summary(horizon_data_majcomps$dbthirdbar_r[grepl('Oe', horizon_data_majcomps$hzname)])
summary(horizon_data_majcomps$om_r[grepl('Oi', horizon_data_majcomps$hzname)])
summary(horizon_data_majcomps$om_r[grepl('Oa', horizon_data_majcomps$hzname)])
summary(horizon_data_majcomps$om_r[grepl('Oe', horizon_data_majcomps$hzname)])
sum(!is.na(horizon_data_majcomps$om_r) & is.na(horizon_data_majcomps$dbthirdbar_r) & horizon_data_majcomps$hzdept_r < 100) #309
table(horizon_data_majcomps$hzname[!is.na(horizon_data_majcomps$om_r) & is.na(horizon_data_majcomps$dbthirdbar_r) & horizon_data_majcomps$hzdept_r < 100])
cokey_test <- unique(horizon_data_majcomps$cokey[!is.na(horizon_data_majcomps$om_r) & is.na(horizon_data_majcomps$dbthirdbar_r) & horizon_data_majcomps$hzdept_r < 100]) #267
sum(mu_data$hectares[mu_data$mukey %in% comp_data$mukey[comp_data$cokey %in% cokey_test]]) #affects 234523.7 hectares
rm(cokey_test)

#convert data.frame to a SoilProfileCollection object
#see https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/aqp-intro.html?root=aqp
horizon_data_spc <- horizon_data_majcomps
depths(horizon_data_spc) <- cokey ~ hzdept_r + hzdepb_r
class(horizon_data_spc)
print(horizon_data_spc)
depth_logic_result <- checkHzDepthLogic(horizon_data_spc)
lapply(depth_logic_result[,2:6], table)
if(!all(depth_logic_result$valid)) {stop(print('There are errors in SSURGO horizonation that need to be fixed.'))}
depth_logic_result[depth_logic_result$valid==FALSE,]

# #fix horizonation for Feb 2022 SSURGO data
# horizon_data_majcomps[horizon_data_majcomps$cokey==21156280,] #no longer irrelevant as missing horizon is at 61-94 cm
# comp_data[comp_data$cokey==21156280,]
# mu_data[mu_data$mukey==3251996,] #2327 hectares
# horizon_data_majcomps[horizon_data_majcomps$cokey==22009763,] #this was removed above for now
# comp_data[comp_data$cokey==22009763,]
# mu_data[mu_data$mukey==471569,] #613 hectares

#estimate soil depth
soil_comp_depths <- profileApply(horizon_data_spc, FUN = estimateSoilDepth)
length(soil_comp_depths) #29306
sum(soil_comp_depths < 30) #1973
sum(soil_comp_depths < 100) #10494
head(soil_comp_depths)
sum(horizon_data_majcomps$om_r==0, na.rm = TRUE) #0 because of fix above

#aggregate horizon data to specified layer
# test <- horizon_to_comp_layer(horizon_SPC = horizon_data_spc[1], depth = 50, comp_df = comp_data, top_depth = 30, layername='30_50')
# rm(test)
comp_SOM <- horizon_to_comp_layer(horizon_SPC = horizon_data_spc, depth = 50, comp_df = comp_data, top_depth = 30, layername='30_50')
comp_SOM$soil_depth <- soil_comp_depths[match(comp_SOM$cokey, names(soil_comp_depths))]
colnames(comp_SOM)
comp_SOM$kgOrg.m2_30_50cm[is.na(comp_SOM$SOM_30_50cm)] <- NA #because otherwise these become 0 due to SOC content calc function that removes NAs

length(unique(comp_SOM$cokey)) #29306, same as above
summary(comp_SOM$kgOrg.m2_30_50cm) #3860 NA
sum(comp_SOM$kgOrg.m2_30_50cm==0, na.rm = TRUE) #45 still reported as zero
comp_SOM[which(comp_SOM$kgOrg.m2_30_50cm==0),]
head(comp_SOM)
lapply(comp_SOM, summary)
summary(comp_SOM$SOM_30_50cm)
sum(is.na(comp_SOM$SOM_30_50cm))
hist(comp_SOM$SOM_30_50cm)
unique(comp_SOM$compname[comp_SOM$kgOrg.m2_30_50cm==0])

#number of instances where soil depth is less than bottom depth of interest and SOM depth is less than reported soil depth
sum(!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$soil_depth <= bottom_depth & comp_SOM$SOM_depth < layer_thickness) #1750 instances
head(comp_SOM[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$soil_depth <= bottom_depth & comp_SOM$SOM_depth < layer_thickness,])

#number of instances where SOC stock is non-NA and SOM depth is greater than reported BD depth
sum(!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$BD_depth) #58 instances
if(!dir.exists(file.path(SSURGOdir, 'intermediate results', subDir_name))) {
  dir.create(file.path(SSURGOdir, 'intermediate results', subDir_name))
}
write.csv(comp_SOM[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$BD_depth,], file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('SOMdepth_greater_BDdepth', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)

cokeys_SOMdepth_greater_BD_depth <- unique(comp_SOM$cokey[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$BD_depth])
test <- horizon_data_majcomps[horizon_data_majcomps$cokey %in% cokeys_SOMdepth_greater_BD_depth,]
dim(test)
test <- test[test$hzdept_r < bottom_depth,]
unique(test$hzname[is.na(test$dbthirdbar_r)])
test <- test[is.na(test$dbthirdbar_r),]
write.csv(test, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('SOMdepth_greater_BDdepth_missingBDhorizons', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
rm(test)
rm(cokeys_SOMdepth_greater_BD_depth)

#convert SOC content to NA when SOM data depth exceeds BD depth
comp_SOM$kgOrg.m2_30_50cm[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$BD_depth] <- NA
sum(comp_SOM$kgOrg.m2_30_50cm==0, na.rm = TRUE) #now 3
comp_SOM[which(comp_SOM$kgOrg.m2_30_50cm==0),]

#convert these to NA
comp_SOM$kgOrg.m2_30_50cm[which(comp_SOM$kgOrg.m2_30_50cm==0 & !is.na(comp_SOM$SOM_30_50cm))] <- NA

#how many instances where SOM data depth exceeds rock frag depth
sum(!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$Rock_Frag_depth) #3
comp_SOM[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$Rock_Frag_depth,]
horizon_data_majcomps[horizon_data_majcomps$cokey %in% comp_SOM$cokey[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$Rock_Frag_depth],]

#convert SOC content to NA when SOM data depth exceeds rock frag depth
comp_SOM$kgOrg.m2_30_50cm[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$Rock_Frag_depth] <- NA
sum(!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$SOM_depth > comp_SOM$Rock_Frag_depth) #0 now

#how many instances of shallow soils with complete SOM data
sum(!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$soil_depth <= bottom_depth & comp_SOM$SOM_depth == layer_thickness) #55

#how many instances of shallow soils with "incomplete" SOM data
sum(!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$soil_depth <= bottom_depth & comp_SOM$SOM_depth < layer_thickness) #1741, but not all of these are of concern

#investigate these cokeys further 
problematic_cokeys <- comp_SOM[!is.na(comp_SOM$kgOrg.m2_30_50cm) & comp_SOM$soil_depth <= bottom_depth & comp_SOM$SOM_depth < layer_thickness,]
dim(problematic_cokeys) #1741 rows
length(unique(problematic_cokeys$cokey)) #1741

write.csv(problematic_cokeys, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_cokeys', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)

problematic_horizons <- horizon_data_majcomps[horizon_data_majcomps$cokey %in% problematic_cokeys$cokey,]

write.csv(problematic_horizons, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_horizons', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
unique(problematic_horizons$hzname)

problematic_horizons_H123 <- problematic_horizons[problematic_horizons$hzname %in% c('H1', 'H2', 'H3', 'H4', 'H5', 'H6'),]
problematic_cokeys_H123 <- problematic_cokeys[problematic_cokeys$cokey %in% problematic_horizons_H123$cokey,]

write.csv(problematic_cokeys_H123, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_cokeys_H123', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
problematic_horizons_H123 <- problematic_horizons[problematic_horizons$cokey %in% problematic_cokeys_H123$cokey, ]
length(unique(problematic_cokeys_H123$cokey)) #359
length(unique(problematic_horizons_H123$cokey)) #359

write.csv(problematic_horizons_H123, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_horizons_H123', fname_depth, Sys.Date(),'.csv')), row.names = FALSE)

dim(problematic_cokeys_H123) #359 that should be dealt with
length(unique(problematic_cokeys_H123$mukey)) #affecting 345 mukeys
sum(mu_data$hectares[mu_data$mukey %in% problematic_cokeys_H123$mukey]) #1,039,816 ha
unique(problematic_horizons_H123$hzname)

problematic_cokeys_noH123 <- problematic_cokeys[!(problematic_cokeys$cokey %in% problematic_cokeys_H123$cokey),]

dim(problematic_cokeys_noH123) #1382 rows
sum(mu_data$hectares[mu_data$mukey %in% problematic_cokeys_noH123$mukey]) #affecting 2,648,930 hectares

problematic_horizons_noH123 <- problematic_horizons[problematic_horizons$cokey %in% problematic_cokeys_noH123$cokey,]

test_cokeys <- problematic_cokeys_noH123$cokey[problematic_cokeys_noH123$cokey %in% problematic_cokeys_H123$cokey]
length(test_cokeys) #no overlap now if 0
rm(test_cokeys)

#further refine noH123 cokeys based on SOM_depth
summary(problematic_cokeys_noH123$soil_depth - depth_top)
problematic_cokeys_noH123[problematic_cokeys_noH123$soil_depth < depth_top,]
horizon_data_majcomps[horizon_data_majcomps$cokey==22460201,] #Br horizon?
horizon_data_majcomps[horizon_data_majcomps$cokey==22460571,] #Br horizon?
horizon_data_majcomps[horizon_data_majcomps$cokey==22588383,] #two Cr horizons
horizon_data_majcomps[horizon_data_majcomps$cokey==22823730,] #two R horizons, makes no sense
summary((problematic_cokeys_noH123$soil_depth - depth_top) == problematic_cokeys_noH123$SOM_depth)
test <- problematic_cokeys_noH123[(problematic_cokeys_noH123$soil_depth - depth_top) > problematic_cokeys_noH123$SOM_depth,]
dim(test)
test
rm(test)

problematic_cokeys_noH123 <- problematic_cokeys_noH123[(problematic_cokeys_noH123$soil_depth - depth_top) > problematic_cokeys_noH123$SOM_depth,]
problematic_horizons_noH123 <- problematic_horizons_noH123[problematic_horizons_noH123$cokey %in% problematic_cokeys_noH123$cokey,]

write.csv(problematic_cokeys_noH123, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_cokeys_noH123', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
write.csv(problematic_horizons_noH123, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_horizons_noH123', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)

#make SOC and SOM NA when the SOM data was less than soil depth when soil depth is < 100
sum(is.na(comp_SOM$kgOrg.m2_30_50cm)) #3924
sum(comp_SOM$cokey %in% problematic_cokeys_noH123$cokey) #17
comp_SOM$kgOrg.m2_30_50cm[comp_SOM$cokey %in% problematic_cokeys_noH123$cokey] <- NA
sum(is.na(comp_SOM$kgOrg.m2_30_50cm)) #3941

sum(is.na(comp_SOM$SOM_30_50cm)) #3860
comp_SOM$SOM_30_50cm[comp_SOM$cokey %in% problematic_cokeys_noH123$cokey] <- NA
sum(is.na(comp_SOM$SOM_30_50cm)) #3877

#look at cases where soil depth > bottom depth but SOM data depth < layer thickness
#will use also component restriction table to make decisions about these cases
sum(!is.na(comp_SOM$kgOrg.m2_30_50cm) & (comp_SOM$soil_depth > bottom_depth) & (comp_SOM$SOM_depth < layer_thickness)) #1809 soils
problematic_cokeys_v2 <- comp_SOM[!is.na(comp_SOM$kgOrg.m2_30_50cm) & (comp_SOM$soil_depth > bottom_depth) & (comp_SOM$SOM_depth < layer_thickness),]
dim(problematic_cokeys_v2) #1809 rows
length(unique(problematic_cokeys_v2$cokey)) #1809

sum(problematic_cokeys_v2$cokey %in% problematic_cokeys$cokey) #0 overlap

problematic_horizons_v2 <- horizon_data_majcomps[horizon_data_majcomps$cokey %in% problematic_cokeys_v2$cokey,]
unique(problematic_horizons_v2$hzname)

problematic_horizons_H123_v2 <- problematic_horizons_v2[problematic_horizons_v2$hzname %in% c('H1', 'H2', 'H3', 'H4', 'H5', 'H6'),]
problematic_cokeys_H123_v2 <- problematic_cokeys_v2[problematic_cokeys_v2$cokey %in% problematic_horizons_H123_v2$cokey,]

write.csv(problematic_cokeys_H123_v2, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_cokeys_H123_v2_', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
problematic_horizons_H123_v2 <- problematic_horizons_v2[problematic_horizons_v2$cokey %in% problematic_cokeys_H123_v2$cokey, ]
length(unique(problematic_cokeys_H123_v2$cokey)) #1571
length(unique(problematic_horizons_H123_v2$cokey)) #1571
unique(problematic_horizons_H123_v2$hzname)
write.csv(problematic_horizons_H123_v2, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_horizons_H123_v2', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
dim(problematic_cokeys_H123_v2) #1571 that should be dealt with
length(unique(problematic_cokeys_H123_v2$mukey)) #affecting 1463 mukeys
sum(mu_data$hectares[mu_data$mukey %in% problematic_cokeys_H123_v2$mukey]) #3,370,817 ha

problematic_cokeys_noH123_v2 <- problematic_cokeys_v2[!(problematic_cokeys_v2$cokey %in% problematic_cokeys_H123_v2$cokey),]
write.csv(problematic_cokeys_noH123_v2, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_cokeys_noH123_v2', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
dim(problematic_cokeys_noH123_v2) #238 rows
sum(problematic_cokeys_H123_v2$cokey %in% problematic_cokeys_noH123_v2$cokey)
sum(mu_data$hectares[mu_data$mukey %in% problematic_cokeys_noH123_v2$mukey]) #affecting  516775.9 hectares

problematic_horizons_noH123_v2 <- problematic_horizons_v2[problematic_horizons_v2$cokey %in% problematic_cokeys_noH123_v2$cokey,]
write.csv(problematic_horizons_noH123_v2, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_horizons_noH123_v2', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)

#which problematic_cokeys_noH123_v2 are actually a problem
sum(!is.na(problematic_cokeys_noH123_v2$SOM_30_50cm) & problematic_cokeys_noH123_v2$SOM_depth < layer_thickness) #238
sum(!is.na(problematic_cokeys_noH123_v2$kgOrg.m2_30_50cm) & problematic_cokeys_noH123_v2$SOM_depth < layer_thickness) #238

sum(problematic_cokeys_H123$cokey %in% problematic_cokeys_H123_v2$cokey)
problematic_cokeys_H123_all <- rbind(problematic_cokeys_H123, problematic_cokeys_H123_v2)
write.csv(problematic_cokeys_H123_all, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_cokeys_H123_all', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
dim(problematic_cokeys_H123_all) #1930 rows
length(unique(problematic_cokeys_H123_all$cokey)) #1930

#read in restriction data
list.files(SSURGOdir)
reskind_df <- read.csv(file.path(SSURGOdir, 'ca_restrictions.csv'), stringsAsFactors = FALSE, na.strings = c('', ' '))
dim(reskind_df)
unique(reskind_df$reskind)
length(unique(reskind_df$cokey)) #25838
reskind_df <- reskind_df[which(reskind_df$reskind != 'Abrupt textural change'),] #this gets rid of NAs also
length(unique(reskind_df$cokey)) #25183
sum(reskind_df$reskind=='Duripan', na.rm = TRUE) #1786
sum(is.na(reskind_df$resdept_r)) #3
reskind_df[is.na(reskind_df$resdept_r),]
reskind_df <- reskind_df[!is.na(reskind_df$resdept_r),]

#split into list by cokey
reskind_list <- split(reskind_df, reskind_df$cokey)
table(sapply(reskind_list, nrow))
test <- do.call(rbind, lapply(reskind_list, function(x) {
  if(nrow(x)==1) {
    x
  } else {
    x[order(x$resdept_r)[1],]
  }}))
sum(test$resdept_r < 100) #25180
head(test)
unique(test$reskind)
length(unique(test$cokey)) #25180
test[test$reskind=='Duripan',]
rm(test)

#problematic horizons
sum(problematic_horizons_H123$cokey %in% problematic_horizons_H123_v2$cokey) #0
problematic_horizons_H123_all <- rbind(problematic_horizons_H123, problematic_horizons_H123_v2)
write.csv(problematic_horizons_H123_all, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_horizons_H123_all', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
head(problematic_horizons_H123_all)
problematic_horizons_H123_all$reskind_match <- FALSE
problematic_horizons_H123_all$reskind_match <- apply(problematic_horizons_H123_all, MARGIN=1, FUN=function(z) {
  if(z[1] %in% as.character(reskind_df$cokey)) {
    y <- reskind_df[as.character(reskind_df$cokey)==z[1],]
    if(trimws(z[3]) %in% as.character(y$resdept_r)) {
      TRUE
    } else {FALSE}
  } else {FALSE}
}
)
length(unique(problematic_horizons_H123_all$cokey)) #1930
table(problematic_horizons_H123_all$reskind_match) #FALSE 4068 TRUE 1912
head(problematic_horizons_H123_all)
length(unique(problematic_horizons_H123_all$cokey[problematic_horizons_H123_all$reskind_match==TRUE]))#1848
problematic_horizons_H123_all[problematic_horizons_H123_all$cokey=='22485463',]
reskind_df[reskind_df==22485463,]
comp_SOM[comp_SOM$cokey==22485463,]
colnames(problematic_cokeys_H123_all)

#determine whether comp data is acceptable
#modified so that SOM_depth is added to depth_top to get whether there's a match to the top of a restrictive horizon
problematic_cokeys_H123_all$SOM_QC <- apply(X=problematic_cokeys_H123_all, MARGIN=1, FUN=function(z) {
  y <- problematic_horizons_H123_all[as.character(problematic_horizons_H123_all$cokey)==z[2], ]
  print(y)
  if(any(y$reskind_match==TRUE)) {
    print('Inside if 1')
    y <- y[y$reskind_match==TRUE,]
    temp <- as.integer(z[7]) + depth_top
    if(trimws(temp) %in% as.character(y$hzdept_r)) {
      print('Inside if 2')
      TRUE
    } else {FALSE}
  } else {FALSE}
})
table(problematic_cokeys_H123_all$SOM_QC) #FALSE 165 TRUE 1765 
problematic_cokeys_H123_all[problematic_cokeys_H123_all$SOM_QC==FALSE,]
problematic_cokeys_H123_all[problematic_cokeys_H123_all$cokey==22484893,]
problematic_horizons_H123_all[problematic_horizons_H123_all$cokey==22484893,]
reskind_df[reskind_df$cokey==22484893,]
write.csv(problematic_cokeys_H123_all, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('problematic_cokeys_H123_FINAL', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)
problematic_cokeys_H123_SOM_NA <- problematic_cokeys_H123_all[problematic_cokeys_H123_all$SOM_QC==FALSE,]
problematic_cokeys_H123_SOM_ck <- problematic_cokeys_H123_all[problematic_cokeys_H123_all$SOM_QC==TRUE,]
if(all(problematic_cokeys_H123_SOM_ck$BD_depth >= problematic_cokeys_H123_SOM_ck$SOM_depth)==TRUE){print('All data with good SOM data also had data to calc SOC content')} else{print('Code needs to be modified here')}
#finalize comp level data now
#make SOC and SOM NA when the SOM data was less than soil depth
sum(is.na(comp_SOM$kgOrg.m2_30_50cm)) #3941
comp_SOM$kgOrg.m2_30_50cm[comp_SOM$cokey %in% problematic_cokeys_H123_SOM_NA$cokey] <- NA
sum(is.na(comp_SOM$kgOrg.m2_30_50cm)) #4106

sum(is.na(comp_SOM$SOM_30_50cm)) #3877
comp_SOM$SOM_30_50cm[comp_SOM$cokey %in% problematic_cokeys_H123_SOM_NA$cokey] <- NA
sum(is.na(comp_SOM$SOM_30_50cm)) #4042

#dealing with problematic_cokeys_noH123_v2 issue, where SOM_depth is less than layer thickness but soil depth is 
dim(problematic_cokeys_noH123_v2) #238
summary(problematic_cokeys_noH123_v2$SOM_30_50cm)
summary(problematic_cokeys_noH123_v2$kgOrg.m2_30_50cm)
summary(problematic_cokeys_noH123_v2$SOM_depth)
head(problematic_cokeys_noH123_v2)
problematic_horizons_noH123_v2[problematic_horizons_noH123_v2$cokey==22460526,]
reskind_df[reskind_df$cokey==22460526,]
problematic_horizons_noH123_v2[problematic_horizons_noH123_v2$cokey==22529890,]
reskind_df[reskind_df$cokey==22529890,]

#make SOC and SOM NA when the SOM data was less than soil depth
sum(is.na(comp_SOM$kgOrg.m2_30_50cm)) #4106
comp_SOM$kgOrg.m2_30_50cm[comp_SOM$cokey %in% problematic_cokeys_noH123_v2$cokey] <- NA
sum(is.na(comp_SOM$kgOrg.m2_30_50cm)) #4344

sum(is.na(comp_SOM$SOM_30_50cm)) #4042
comp_SOM$SOM_30_50cm[comp_SOM$cokey %in% problematic_cokeys_noH123_v2$cokey] <- NA
sum(is.na(comp_SOM$SOM_30_50cm)) #4280

#cokeys have been updated
#sum((thick / 10) * (x$om_r / om_to_c) * x$dbthirdbar_r * (1 - x$fragvol_r_sum / 100)
# horizon_data_majcomps[horizon_data_majcomps$cokey==21292622,]
# comp_SOM[comp_SOM$cokey==21292622,]
# 80*5/100+5*28/100+1.5*67/100 #6.405 is correct
# sum(c((5/10)*(80/1.72) * 0.05 *(1-33/100), (28/10)*(5/1.72) * 1.18 *(1-15/100), (67/10)*(1.5/1.72) * 1.40 *(1-75/100))) #10.98808 is correct
# 
# horizon_data_majcomps[horizon_data_majcomps$cokey==21158686,]
# comp_SOM[comp_SOM$cokey==21158686,]
# 
# horizon_data_majcomps[horizon_data_majcomps$cokey==21139431,]
# comp_SOM[comp_SOM$cokey==21139431,]

write.csv(comp_SOM, file.path(SSURGOdir, 'intermediate results', subDir_name, paste0('comp_SOM', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)

comppct_by_mukey_om_data <- data.frame(mukey=row.names(tapply(comp_SOM$comppct[!is.na(comp_SOM$SOM_30_50cm)], comp_SOM$mukey[!is.na(comp_SOM$SOM_30_50cm)], sum)), comppct_tot=as.numeric(tapply(comp_SOM$comppct[!is.na(comp_SOM$SOM_30_50cm)], comp_SOM$mukey[!is.na(comp_SOM$SOM_30_50cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_om_data$comppct_tot)
sum(comppct_by_mukey_om_data$comppct_tot < 70) #2672 mukeys less than 70
sum(comppct_by_mukey_om_data$comppct_tot==15) #85 have only 15% of components with data
sum(comppct_by_mukey_om_data$comppct_tot < 15) #0, because major component defined as >15% of mapunit

comppct_by_mukey_soc_data <- data.frame(mukey=row.names(tapply(comp_SOM$comppct[!is.na(comp_SOM$kgOrg.m2_30_50cm)], comp_SOM$mukey[!is.na(comp_SOM$kgOrg.m2_30_50cm)], sum)), comppct_tot=as.numeric(tapply(comp_SOM$comppct[!is.na(comp_SOM$kgOrg.m2_30_50cm)], comp_SOM$mukey[!is.na(comp_SOM$kgOrg.m2_30_50cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_soc_data$comppct_tot)

CA_30_50cm_muagg <- MUAggregate_wrapper(df1=comp_SOM, varnames = c('SOM_30_50cm', 'kgOrg.m2_30_50cm'))
CA_30_50cm_muagg$compct_SOM_30_50cm <- comppct_by_mukey_om_data$comppct_tot[match(CA_30_50cm_muagg$mukey, comppct_by_mukey_om_data$mukey)]
CA_30_50cm_muagg$compct_SOCkg_30_50cm <- comppct_by_mukey_soc_data$comppct_tot[match(CA_30_50cm_muagg$mukey, comppct_by_mukey_soc_data$mukey)]
CA_30_50cm_muagg$hectares <- mu_data$hectares[match(CA_30_50cm_muagg$mukey, mu_data$mukey)]
summary(CA_30_50cm_muagg$hectares)
sum(CA_30_50cm_muagg$hectares[which(CA_30_50cm_muagg$compct_SOM==15.00)]) #269429.9
sum(CA_30_50cm_muagg$hectares) #32,827,071 total hectares
sum(CA_30_50cm_muagg$hectares[!is.na(CA_30_50cm_muagg$SOM_30_50cm)]) #29,196,422
sum(mu_data$hectares) #41,477,744

#add result to file with all mukeys
mu_data$SOM_30_50cm <- CA_30_50cm_muagg$SOM_30_50cm[match(mu_data$mukey, CA_30_50cm_muagg$mukey)]
mu_data$compct_SOM_30_50cm <- CA_30_50cm_muagg$compct_SOM_30_50cm[match(mu_data$mukey, CA_30_50cm_muagg$mukey)]
mu_data$kgSOC.m2_30_50cm <- CA_30_50cm_muagg$kgOrg.m2_30_50cm[match(mu_data$mukey, CA_30_50cm_muagg$mukey)]
mu_data$compct_SOCkg_30_50cm <- CA_30_50cm_muagg$compct_SOCkg_30_50cm[match(mu_data$mukey, CA_30_50cm_muagg$mukey)]

sum(mu_data$hectares[!is.na(mu_data$SOM_30_50cm)])/sum(mu_data$hectares) #70.4% of AOI with SOM summaries
length(mu_data$mukey[!is.na(mu_data$SOM_30_50cm)])/length(mu_data$mukey) #89.4% of mukeys
write.csv(mu_data, file.path(SSURGOdir, 'summaries', paste0('CA_mapunit_SOM', fname_depth, Sys.Date(), '.csv')), row.names = FALSE)

#cokeys have been updated
#manual check of calculations
# mu_data[mu_data$mukey==455490,]
# comp_data[comp_data$mukey==455490,]
# horizon_data_majcomps[horizon_data_majcomps$cokey==22128903,]
# 0.5*(20/100)+0.25*(80/100) #0.3 % SOM
# comp_SOM[comp_SOM$cokey==22128903,]
# horizon_data_majcomps[horizon_data_majcomps$cokey==22128904,]
# 0.5*(15/100)+0.25*(85/100) #0.2875 % SOM
# comp_SOM[comp_SOM$cokey==22128904,]
# (0.5*(20/100)+0.25*(80/100))*(55/80) + (0.5*(15/100)+0.25*(85/100))*(25/80) #0.2960938 % SOM for mapunit major component percent weighted avg.
# mu_data$SOM_30_50cm[mu_data$mukey==455490] #0.2960937, correct
# 
# mu_data[mu_data$mukey==457102,]
# comp_data[comp_data$mukey==457102,]
# horizon_data_majcomps[horizon_data_majcomps$cokey==19720470,] #2.5%
# horizon_data_majcomps[horizon_data_majcomps$cokey==19720469,] #1.5%
# 2.5*50/95+1.5*45/95
# mu_data$SOM_30_50cm[mu_data$mukey==457102] #correct

#identify some breakpoints
SOM_wtd_data <- unlist(mapply(function(x, y) {rep(x, times=round(y/10, 0))}, x=mu_data$SOM_30_50cm, y=mu_data$hectares))

length(SOM_wtd_data) #3247044
writeClipboard(as.character(quantile(SOM_wtd_data, na.rm = TRUE, probs = c(0.25, 0.5, 0.75))))
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.2, 0.4, 0.6, 0.8), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.167, 0.334, 0.501, 0.668, 0.835), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.143, 0.286, 0.429, 0.572, 0.715, 0.858), na.rm = TRUE)))

SOC_wtd_data <- unlist(mapply(function(x, y) {rep(x, times=round(y/10, 0))}, x=mu_data$kgSOC.m2_30_50cm, y=mu_data$hectares))
head(SOC_wtd_data, 20)
writeClipboard(as.character(quantile(SOC_wtd_data, na.rm = TRUE, probs = c(0.25, 0.5, 0.75))))
writeClipboard(as.character(quantile(SOC_wtd_data, probs = c(0.2, 0.4, 0.6, 0.8), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOC_wtd_data, probs = c(0.167, 0.334, 0.501, 0.668, 0.835), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOC_wtd_data, probs = c(0.143, 0.286, 0.429, 0.572, 0.715, 0.858), na.rm = TRUE)))

#manual check of spatially weighted 40th percentile
sum(mu_data$hectares[which(mu_data$SOM_30_50cm < 1.1911765)]) #11,866,719 hectares
11866719/sum(mu_data$hectares[!is.na(mu_data$SOM_30_50cm)])
#40% of area of interest with SOM data < 1.191%

#compare with 30-cm aggregation
list.files(file.path(SSURGOdir, 'summaries'))
mu_data_30cm <- read.csv(file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_30cm_10.7.22.csv'), stringsAsFactors = FALSE)
all(mu_data$mukey==mu_data_30cm$mukey)
mu_diffs <- mu_data_30cm[!is.na(mu_data_30cm$SOM_30cm) & is.na(mu_data$SOM_30_50cm),]
dim(mu_diffs) #1557
head(mu_diffs, 20)
mu_data_30cm[mu_data_30cm$mukey==1034683,]
mu_data[mu_data$mukey==1034683,]
comp_SOM[comp_SOM$mukey==1034683,]
horizon_data_majcomps[horizon_data_majcomps$cokey==21501161,]
reskind_df[reskind_df$cokey==21501161,]
horizon_data_majcomps[horizon_data_majcomps$cokey==21501162,]
reskind_df[reskind_df$cokey==21501162,]
horizon_data_majcomps[horizon_data_majcomps$cokey==21501163,]
reskind_df[reskind_df$cokey==21501163,]
tail(mu_diffs, 20)

#om data goes to 46 cm but lithic indicated at H3 at 66 cm
mu_data_30cm[mu_data_30cm$mukey==848226,]
mu_data[mu_data$mukey==848226,]
comp_SOM[comp_SOM$mukey==848226,]
horizon_data_majcomps[horizon_data_majcomps$cokey==22011423,]
reskind_df[reskind_df$cokey==22011423,]
