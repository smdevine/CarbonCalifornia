#updated to include depth-weighted SOM for soils < 30 cm deep
#TO-DO: (1) add soil depth calc to map-units (2) convert mukey to character class
library(aqp)
mainDir <- 'C:/Users/smdevine/Desktop/post doc/soil organic carbon'
SSURGOdir <- file.path(mainDir, 'SSURGO data')
list.files(SSURGOdir)
rockNA_to_0 <- TRUE #to convert all NA for rock frag content to 0

#read-in mapunit area info
mu_data <- read.csv(file.path(SSURGOdir, 'ca_mapunit_area.csv'), stringsAsFactors = FALSE)
mu_data$hectares <- mu_data$sqm/10000
sum(mu_data$hectares) #41,477,744 hectares

#read-in component data and define major component as >15% aerial coverage
list.files(SSURGOdir)
comp_data <- read.csv(file.path(SSURGOdir, 'ca_components.csv'), na.strings = c('', ' '), stringsAsFactors = FALSE)
lapply(comp_data, summary)
length(unique(comp_data$cokey)) #91970 unique cokeys
length(unique(comp_data$mukey)) #19005 unique mukeys
if(sum(is.na(comp_data$majcompflag)) > 0) {stop(print('there are NAs in majcomp column!'))}
if(sum(is.na(comp_data$comppct_r)) > 0) {stop(print('there are NAs in the comppct column!'))}
sum(comp_data$comppct_r[comp_data$majcompflag=='Yes'] < 15) #150 are major components with <15% aereal coverage
sum(comp_data$comppct_r[comp_data$majcompflag=='No '] >= 15) #285 are minor components with >15% aereal coverage
comp_data$majcompflag[comp_data$majcompflag=='No ' & comp_data$comppct_r>=15] <- 'Yes'
comp_data$majcompflag[comp_data$majcompflag=='Yes' & comp_data$comppct_r < 15] <- 'No '

#read-in soil horizon data and convert all 0% SOM reporting to NA
horizon_data <- read.csv(file.path(SSURGOdir, 'ca_horizons.csv'), na.strings = c('', ' '), stringsAsFactors = FALSE)
colnames(horizon_data)
length(unique(horizon_data$cokey)) #41772 unique cokeys
sum(horizon_data$om_r==0, na.rm = TRUE) #4294 instances of 0% SOM
hznames_zeroOM <- unique(horizon_data$hzname[which(horizon_data$om_r==0)])
hznames_zeroOM
horizon_data$om_r[horizon_data$om_r==0] <- NA
horizon_data$majcompflag <- comp_data$majcompflag[match(horizon_data$cokey, comp_data$cokey)]
table(horizon_data$majcompflag)
horizon_data_majcomps <- horizon_data[horizon_data$majcompflag=='Yes',]
horizon_data_majcomps$mukey <- comp_data$mukey[match(horizon_data_majcomps$cokey, comp_data$cokey)]
horizon_data_majcomps <- subset(horizon_data_majcomps, horizon_data_majcomps$cokey!=22009763) #see below for rationale

#inspect bulk density data
summary(horizon_data_majcomps$dbthirdbar_r)
sum(horizon_data_majcomps$dbthirdbar_r < 0.5, na.rm = TRUE) #2866 have BD < 0.5 g cm^3

#inspect rock fragment data
summary(horizon_data_majcomps$fragvol_r_sum)
sum(horizon_data_majcomps$fragvol_r_sum > 100, na.rm = TRUE) #12 are greater than 0
horizon_data_majcomps[which(horizon_data_majcomps$fragvol_r_sum > 100),] #12 are 
sum(horizon_data_majcomps$fragvol_r_sum == 100, na.rm = TRUE) #7 equal to zero
sum(horizon_data_majcomps$fragvol_r_sum==0, na.rm = TRUE) #4671
sum(is.na(horizon_data_majcomps$fragvol_r_sum)) #28539
sum(is.na(horizon_data_majcomps$fragvol_r_sum) & !is.na(horizon_data_majcomps$om_r))
sum(is.na(horizon_data_majcomps$dbthirdbar_r) & !is.na(horizon_data_majcomps$om_r)) #361

#convert rock fragment NA to 0% rock fragment [don't do for now]
if(rockNA_to_0) {
  horizon_data_majcomps$fragvol_r_sum[is.na(horizon_data_majcomps$fragvol_r_sum)] <- 0
}

#then convert rock fragment data >= 100% to NA
horizon_data_majcomps$fragvol_r_sum[horizon_data_majcomps$fragvol_r_sum >= 100] <- NA
summary(horizon_data_majcomps$fragvol_r_sum)
duripans <- horizon_data_majcomps[which(horizon_data_majcomps$hzname=='Bqm'),]
summary(duripans$om_r)
hist(duripans$hzdept_r)
sum(duripans$hzdept_r < 30) #3
duripans[duripans$hzdept_r < 30,]
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

#fix horizonation?
horizon_data[horizon_data$cokey==21156280,] #irrelevant as missing horizon is at 61-94 cm
horizon_data[horizon_data$cokey==22009763,] #this was removed above for now

#estimate soil depth
soil_comp_depths <- profileApply(horizon_data_spc, FUN = estimateSoilDepth)
sum(soil_comp_depths < 30) #1950
head(soil_comp_depths)
sum(horizon_data_majcomps$om_r==0, na.rm = TRUE) #0 because of fix above

#functions to work with ProfileApply
wtd.mean <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzdepb_r - x$hzdept_r
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=TRUE) #modified to TRUE on 2/16/22
  m
}
kgOrgC_sum <- function(x, slice_it=FALSE, depth, rm.NAs=TRUE, om_to_c=1.72) { #removing NAs so as to capture SOC conent in soils <30 cm deep; but true NAs must be dealt with below based on comparison of "soil depth" and depth to which SOM data exists
  if (slice_it) { 
    x <- horizons(x)[1:depth, ]
    depths(x) <- cokey ~ hzdept_r + hzdepb_r
  }
  thick <- x$hzdepb_r - x$hzdept_r
  sum((thick / 10) * (x$om_r / om_to_c) * x$dbthirdbar_r * (1 - x$fragvol_r_sum / 100), na.rm = rm.NAs)
}
# kgOrgC <- function(x, om_to_c=1.72) {
#   thick <- x$hzdepb_r - x$hzdept_r
#   (thick / 10) * (x$om_r / om_to_c) * x$dbthirdbar_r * (1 - x$fragvol_r_sum / 100)
# }
data_depth <- function(x, varname) {
  thick <- x$hzdepb_r - x$hzdept_r
  sum(thick[!is.na(x[[varname]])])
}
horizon_to_comp <- function(horizon_SPC, depth, comp_df, vars_of_interest = 'om_r', varnames = 'SOM') {
  columnames <- paste0(varnames, '_', depth, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  sliced_SPC <- slice(horizon_SPC, 0:(depth-1) ~ .) #depth was '0:depth' in previous version
  stopifnot(unique(sliced_SPC$pedon_key)==site(sliced_SPC)$pedon_key)
  for (i in seq_along(vars_of_interest)) {
    s <- site(sliced_SPC)
    s[[columnames[i]]] <- profileApply(sliced_SPC, FUN = wtd.mean, y=vars_of_interest[i])
    site(sliced_SPC) <- s
  }
  s <- site(sliced_SPC)
  s[[paste0('kgOrg.m2_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum)
  s[['SOM_depth']] <- profileApply(sliced_SPC, FUN=data_depth, varname='om_r')
  s[['BD_depth']] <- profileApply(sliced_SPC, FUN=data_depth, varname='dbthirdbar_r')
  s[['Rock_Frag_depth']] <- profileApply(sliced_SPC, FUN = data_depth, varname='fragvol_r_sum')
  columnames <- c(columnames, paste0('kgOrg.m2_', depth, 'cm'), 'SOM_depth', 'BD_depth', 'Rock_Frag_depth')
  rm(depth, envir = .GlobalEnv) #because we had to put it there earlier
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

length(unique(comp_data$cokey))
comp_SOM_30cm <- horizon_to_comp(horizon_SPC = horizon_data_spc, depth = 30, comp_df = comp_data)
comp_SOM_30cm$soil_depth <- soil_comp_depths[match(comp_SOM_30cm$cokey, names(soil_comp_depths))]
comp_SOM_30cm$kgOrg.m2_30cm[is.na(comp_SOM_30cm$SOM_30cm)] <- NA #because otherwise these become 0 due to SOC content calc function that removes NAs

length(unique(comp_SOM_30cm$cokey))
summary(comp_SOM_30cm$kgOrg.m2_30cm)
sum(comp_SOM_30cm$kgOrg.m2_30cm==0, na.rm = TRUE) #45 still reported as zero
comp_SOM_30cm[which(comp_SOM_30cm$kgOrg.m2_30cm==0),]
head(comp_SOM_30cm)
lapply(comp_SOM_30cm, summary)
summary(comp_SOM_30cm$SOM_30cm)
sum(is.na(comp_SOM_30cm$SOM_30cm))
hist(comp_SOM_30cm$SOM_30cm)
unique(comp_SOM_30cm$compname[comp_SOM_30cm$kgOrg.m2_30cm==0])
sum(!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$soil_depth <= 30 & comp_SOM_30cm$SOM_depth < comp_SOM_30cm$soil_depth) #178 instances
comp_SOM_30cm[!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$soil_depth <= 30 & comp_SOM_30cm$SOM_depth < comp_SOM_30cm$soil_depth,]
sum(!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$SOM_depth > comp_SOM_30cm$BD_depth) #204 instances
comp_SOM_30cm[!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$SOM_depth > comp_SOM_30cm$BD_depth,]

#how many instances of shallow soils with complete SOM data
sum(!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$soil_depth < 30 & comp_SOM_30cm$SOM_depth == comp_SOM_30cm$soil_depth) #572
sum(comp_SOM_30cm$kgOrg.m2_30cm==0, na.rm = TRUE)

#convert SOC content to NA when SOM data depth exceeds BD depth
comp_SOM_30cm$kgOrg.m2_30cm[!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$SOM_depth > comp_SOM_30cm$BD_depth] <- NA
sum(comp_SOM_30cm$kgOrg.m2_30cm==0, na.rm = TRUE) #now only 1
sum(!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$soil_depth <= 30 & comp_SOM_30cm$SOM_depth < comp_SOM_30cm$soil_depth) #still 178 of these

#investigate these 178 cokeys further 
problematic_cokeys <- comp_SOM_30cm[!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$soil_depth <=30 & comp_SOM_30cm$SOM_depth < comp_SOM_30cm$soil_depth,]
write.csv(problematic_cokeys, file.path(SSURGOdir, 'intermediate results', 'problematic_cokeys_2_24_22.csv'), row.names = FALSE)
problematic_horizons <- horizon_data_majcomps[horizon_data_majcomps$cokey %in% problematic_cokeys$cokey,]
write.csv(problematic_horizons, file.path(SSURGOdir, 'intermediate results', 'problematic_horizons_2_24_22.csv'), row.names = FALSE)
unique(problematic_horizons$hzname)

problematic_horizons_noH123 <- problematic_horizons[!(problematic_horizons$hzname %in% c('H1', 'H2', 'H3')),]
write.csv(problematic_horizons_noH123, file.path(SSURGOdir, 'intermediate results', 'problematic_horizons_noH123_2_24_22.csv'), row.names = FALSE)
problematic_cokeys_noH123 <- problematic_cokeys[problematic_cokeys$cokey %in% problematic_horizons_noH123$cokey,]
write.csv(problematic_cokeys_noH123, file.path(SSURGOdir, 'intermediate results', 'problematic_cokeys_noH123_2_24_22.csv'), row.names = FALSE)

problematic_horizons_H123 <- problematic_horizons[problematic_horizons$hzname %in% c('H1', 'H2', 'H3'),]
write.csv(problematic_horizons_H123, file.path(SSURGOdir, 'intermediate results', 'problematic_horizons_H123_2_24_22.csv'), row.names = FALSE)
problematic_cokeys_H123 <- problematic_cokeys[problematic_cokeys$cokey %in% problematic_horizons_H123$cokey,]
write.csv(problematic_cokeys_H123, file.path(SSURGOdir, 'intermediate results', 'problematic_cokeys_H123_2_24_22.csv'), row.names = FALSE)
dim(problematic_cokeys_H123) #141 that should be dealt with
length(unique(problematic_cokeys_H123$mukey)) #affecting 137 mukeys
sum(mu_data$hectares[mu_data$mukey %in% problematic_cokeys_H123$mukey]) #444,986 ha

#make SOC and SOM NA when the SOM data was less than soil depth
sum(is.na(comp_SOM_30cm$kgOrg.m2_30cm)) #1922
comp_SOM_30cm$kgOrg.m2_30cm[comp_SOM_30cm$cokey %in% problematic_cokeys_noH123$cokey] <- NA
sum(is.na(comp_SOM_30cm$kgOrg.m2_30cm)) #1960

sum(is.na(comp_SOM_30cm$SOM_30cm)) #1718
comp_SOM_30cm$SOM_30cm[comp_SOM_30cm$cokey %in% problematic_cokeys_noH123$cokey] <- NA
sum(is.na(comp_SOM_30cm$SOM_30cm)) #1756

#what about when SOM_depth is less than 30 but estimated soil depth is greater than 30 (many of these probably still H1, H2, H3 terminology)
sum(!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$soil_depth > 30 & comp_SOM_30cm$SOM_depth < 30) #1346
problematic_cokeys_v2 <- comp_SOM_30cm[!is.na(comp_SOM_30cm$kgOrg.m2_30cm) & comp_SOM_30cm$soil_depth > 30 & comp_SOM_30cm$SOM_depth < 30,]
sum(problematic_cokeys_v2$cokey %in% problematic_cokeys$cokey)

problematic_horizons_v2 <- horizon_data_majcomps[horizon_data_majcomps$cokey %in% problematic_cokeys_v2$cokey,]
unique(problematic_horizons_v2$hzname)

problematic_horizons_noH123_v2 <- problematic_horizons_v2[!(problematic_horizons_v2$hzname %in% c('H1', 'H2', 'H3', 'H4', 'H5')),]
write.csv(problematic_horizons_noH123_v2, file.path(SSURGOdir, 'intermediate results', 'problematic_horizons_noH123_2_24_22_v2.csv'), row.names = FALSE)
problematic_cokeys_noH123_v2 <- problematic_cokeys_v2[problematic_cokeys_v2$cokey %in% problematic_horizons_noH123_v2$cokey,]
write.csv(problematic_cokeys_noH123_v2, file.path(SSURGOdir, 'intermediate results', 'problematic_cokeys_noH123_2_24_22_v2.csv'), row.names = FALSE)
dim(problematic_cokeys_noH123_v2) #534

problematic_horizons_H123_v2 <- problematic_horizons_v2[problematic_horizons_v2$hzname %in% c('H1', 'H2', 'H3', 'H4', 'H5'),]
write.csv(problematic_horizons_H123_v2, file.path(SSURGOdir, 'intermediate results', 'problematic_horizons_H123_2_24_22_v2.csv'), row.names = FALSE)
problematic_cokeys_H123_v2 <- problematic_cokeys_v2[problematic_cokeys_v2$cokey %in% problematic_horizons_H123_v2$cokey,]
write.csv(problematic_cokeys_H123_v2, file.path(SSURGOdir, 'intermediate results', 'problematic_cokeys_H123_2_24_22_v2.csv'), row.names = FALSE)

problematic_cokeys_H123_all <- rbind(problematic_cokeys_H123, problematic_cokeys_H123_v2)
write.csv(problematic_cokeys_H123_all, file.path(SSURGOdir, 'intermediate results', 'problematic_cokeys_H123_all_2_24_22.csv'), row.names = FALSE)
dim(problematic_cokeys_H123_all) #956
length(unique(problematic_cokeys_H123_all$cokey))

#read in restriction data
list.files(SSURGOdir)
reskind_df <- read.csv(file.path(SSURGOdir, 'problematic_cokeys_H123_2_24_22_v2.csv'), stringsAsFactors = FALSE, na.strings = c('', ' '))
dim(reskind_df)
unique(reskind_df$reskind)
reskind_df <- reskind_df[which(reskind_df$reskind != 'Abrupt textural change'),] #this gets rid of NAs also
length(unique(reskind_df$cokey)) #791
sum(reskind_df$reskind=='Duripan', na.rm = TRUE) #86

#split into list by cokey
reskind_list <- split(reskind_df, reskind_df$cokey)
table(sapply(reskind_list, nrow))
test <- do.call(rbind, lapply(reskind_list, function(x) {
  if(nrow(x)==1) {
    x
  } else {
    x[order(x$resdept_r)[1],]
}}))
sum(test$resdept_r < 30) #650
head(test)
unique(test$reskind)
length(unique(test$cokey)) #791
test[test$reskind=='Duripan',]

#problematic horizons
problematic_horizons_H123_all <- rbind(problematic_horizons_H123, problematic_horizons_H123_v2)
write.csv(problematic_horizons_H123_all, file.path(SSURGOdir, 'intermediate results', 'problematic_horizons_H123_2_24_22.csv'), row.names = FALSE)
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
table(problematic_horizons_H123_all$reskind_match) #832 matches
head(problematic_horizons_H123_all)
length(unique(problematic_horizons_H123_all$cokey[problematic_horizons_H123_all$reskind_match==TRUE]))
head(problematic_cokeys_H123_all)
problematic_horizons_H123_all[problematic_horizons_H123_all$cokey=='21096542',]

#determine whether comp data is acceptable
problematic_cokeys_H123_all$SOM_QC <- apply(X=problematic_cokeys_H123_all, MARGIN=1, FUN=function(z) {
  y <- problematic_horizons_H123_all[as.character(problematic_horizons_H123_all$cokey)==z[2], ]
  print(y)
  if(any(y$reskind_match==TRUE)) {
    print('Inside if 1')
    y <- y[y$reskind_match==TRUE,]
    if(trimws(z[7]) %in% as.character(y$hzdept_r)) {
      print('Inside if 2')
      TRUE
    } else {FALSE}
  } else {FALSE}
})
table(problematic_cokeys_H123_all$SOM_QC)
write.csv(problematic_cokeys_H123_all, file.path(SSURGOdir, 'intermediate results', 'problematic_cokeys_H123_all_FINAL.csv'), row.names = FALSE)
problematic_cokeys_H123_SOM_NA <- problematic_cokeys_H123_all[problematic_cokeys_H123_all$SOM_QC==FALSE,]
problematic_cokeys_H123_SOM_ck <- problematic_cokeys_H123_all[problematic_cokeys_H123_all$SOM_QC==TRUE,]
if(all(problematic_cokeys_H123_SOM_ck$BD_depth >= problematic_cokeys_H123_SOM_ck$SOM_depth)==TRUE){print('All data with good SOM data also had data to calc SOC content')} else{print('Code needs to be modified here')}
#finalize comp level data now
#make SOC and SOM NA when the SOM data was less than soil depth
sum(is.na(comp_SOM_30cm$kgOrg.m2_30cm)) #1960
comp_SOM_30cm$kgOrg.m2_30cm[comp_SOM_30cm$cokey %in% problematic_cokeys_H123_SOM_NA$cokey] <- NA
sum(is.na(comp_SOM_30cm$kgOrg.m2_30cm)) #2271

sum(is.na(comp_SOM_30cm$SOM_30cm)) #1756
comp_SOM_30cm$SOM_30cm[comp_SOM_30cm$cokey %in% problematic_cokeys_H123_SOM_NA$cokey] <- NA
sum(is.na(comp_SOM_30cm$SOM_30cm)) #2067


horizon_data_majcomps[horizon_data_majcomps$cokey==21292622,]
comp_SOM_30cm[comp_SOM_30cm$cokey==21292622,]

horizon_data_majcomps[horizon_data_majcomps$cokey==21158686,]
comp_SOM_30cm[comp_SOM_30cm$cokey==21158686,]

horizon_data_majcomps[horizon_data_majcomps$cokey==21139431,]
comp_SOM_30cm[comp_SOM_30cm$cokey==21139431,]

write.csv(comp_SOM_30cm, file.path(SSURGOdir, 'intermediate results', 'comp_SOM_30cm_2_25_22.csv'), row.names = FALSE)

comppct_by_mukey_om_data <- data.frame(mukey=row.names(tapply(comp_SOM_30cm$comppct[!is.na(comp_SOM_30cm$SOM_30cm)], comp_SOM_30cm$mukey[!is.na(comp_SOM_30cm$SOM_30cm)], sum)), comppct_tot=as.numeric(tapply(comp_SOM_30cm$comppct[!is.na(comp_SOM_30cm$SOM_30cm)], comp_SOM_30cm$mukey[!is.na(comp_SOM_30cm$SOM_30cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_om_data$comppct_tot)
sum(comppct_by_mukey_om_data$comppct_tot < 70) #2188 mukeys less than 70
sum(comppct_by_mukey_om_data$comppct_tot==15) #40 have only 15% of components with data
sum(comppct_by_mukey_om_data$comppct_tot < 15) #0, because major component defined as >15% of mapunit

comppct_by_mukey_soc_data <- data.frame(mukey=row.names(tapply(comp_SOM_30cm$comppct[!is.na(comp_SOM_30cm$kgOrg.m2_30cm)], comp_SOM_30cm$mukey[!is.na(comp_SOM_30cm$kgOrg.m2_30cm)], sum)), comppct_tot=as.numeric(tapply(comp_SOM_30cm$comppct[!is.na(comp_SOM_30cm$kgOrg.m2_30cm)], comp_SOM_30cm$mukey[!is.na(comp_SOM_30cm$kgOrg.m2_30cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_soc_data$comppct_tot)


CA_30cm_muagg <- MUAggregate_wrapper(df1=comp_SOM_30cm, varnames = c('SOM_30cm', 'kgOrg.m2_30cm'))
CA_30cm_muagg$compct_SOM <- comppct_by_mukey_om_data$comppct_tot[match(CA_30cm_muagg$mukey, comppct_by_mukey_om_data$mukey)]
CA_30cm_muagg$compct_SOCkg <- comppct_by_mukey_soc_data$comppct_tot[match(CA_30cm_muagg$mukey, comppct_by_mukey_soc_data$mukey)]
CA_30cm_muagg$hectares <- mu_data$hectares[match(CA_30cm_muagg$mukey, mu_data$mukey)]
summary(CA_30cm_muagg$hectares)
sum(CA_30cm_muagg$hectares[which(CA_30cm_muagg$compct_SOM==15.00)]) #70920.33
sum(CA_30cm_muagg$hectares) #32,811,572 total hectares
sum(CA_30cm_muagg$hectares[!is.na(CA_30cm_muagg$SOM_30cm)]) #31,460,565
sum(mu_data$hectares) #41477744

#add result to file with all mukeys
mu_data$SOM_30cm <- CA_30cm_muagg$SOM_30cm[match(mu_data$mukey, CA_30cm_muagg$mukey)]
mu_data$compct_SOM <- CA_30cm_muagg$compct_SOM[match(mu_data$mukey, CA_30cm_muagg$mukey)]
mu_data$kgSOC.m2_30cm <- CA_30cm_muagg$kgOrg.m2_30cm[match(mu_data$mukey, CA_30cm_muagg$mukey)]
mu_data$compct_SOCkg <- CA_30cm_muagg$compct_SOCkg[match(mu_data$mukey, CA_30cm_muagg$mukey)]

sum(mu_data$hectares[!is.na(mu_data$SOM_30cm)])/sum(mu_data$hectares) #75.8% of AOI with SOM summaries
length(mu_data$mukey[!is.na(mu_data$SOM_30cm)])/length(mu_data$mukey) #94.6% of mukeys
write.csv(mu_data, file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_30cm_2.25.22.csv'), row.names = FALSE)

#manual check of calculations
mu_data[mu_data$mukey==455490,]
comp_data[comp_data$mukey==455490,]
horizon_data_majcomps[horizon_data_majcomps$cokey==19535419,]
0.5*(20/30)+0.25*(10/30) #0.4166667 % SOM
comp_SOM_30cm[comp_SOM_30cm$cokey==19535419,]
horizon_data_majcomps[horizon_data_majcomps$cokey==19535420,]
0.5*(15/30)+0.25*(15/30) #0.375 % SOM
(0.5*(20/30)+0.25*(10/30))*(55/80) + (0.5*(15/30)+0.25*(15/30))*(25/80) #0.4036458 % SOM for mapunit major component percent weighted avg.
mu_data$SOM_30cm[mu_data$mukey==455490] #0.4036458, correct

mu_data[mu_data$mukey==457102,]
comp_data[comp_data$mukey==457102,]
horizon_data_majcomps[horizon_data_majcomps$cokey==19720470,] #2.5%
horizon_data_majcomps[horizon_data_majcomps$cokey==19720469,] #1.5%
2.5*50/95+1.5*45/95
mu_data$SOM_30cm[mu_data$mukey==457102] #correct

#identify some breakpoints
SOM_wtd_data <- unlist(mapply(function(x, y) {rep(x, times=round(y/10, 0))}, x=mu_data$SOM_30cm, y=mu_data$hectares))

length(SOM_wtd_data) #3247044
writeClipboard(as.character(quantile(SOM_wtd_data, na.rm = TRUE, probs = c(0.25, 0.5, 0.75))))
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.2, 0.4, 0.6, 0.8), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.167, 0.334, 0.501, 0.668, 0.835), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.143, 0.286, 0.429, 0.572, 0.715, 0.858), na.rm = TRUE)))

SOC_wtd_data <- unlist(mapply(function(x, y) {rep(x, times=round(y/10, 0))}, x=mu_data$kgSOC.m2_30cm, y=mu_data$hectares))
head(SOC_wtd_data, 20)
writeClipboard(as.character(quantile(SOC_wtd_data, na.rm = TRUE, probs = c(0.25, 0.5, 0.75))))
writeClipboard(as.character(quantile(SOC_wtd_data, probs = c(0.2, 0.4, 0.6, 0.8), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOC_wtd_data, probs = c(0.167, 0.334, 0.501, 0.668, 0.835), na.rm = TRUE)))
writeClipboard(as.character(quantile(SOC_wtd_data, probs = c(0.143, 0.286, 0.429, 0.572, 0.715, 0.858), na.rm = TRUE)))



#manual check of spatially weighted 40th percentile
sum(mu_data$hectares[which(mu_data$SOM_30cm < 1.1911765)]) #11,866,719 hectares
11866719/sum(mu_data$hectares[!is.na(mu_data$SOM_30cm)])
#40% of area of interest with SOM data < 1.191%