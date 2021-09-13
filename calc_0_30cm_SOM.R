library(aqp)
mainDir <- 'C:/Users/smdevine/Desktop/post doc/soil organic carbon'
SSURGOdir <- file.path(mainDir, 'SSURGO data')
list.files(SSURGOdir)

#read-in mapunit area info
mu_data <- read.csv(file.path(SSURGOdir, 'ca_mapunit_area.csv'), stringsAsFactors = FALSE)
mu_data$hectares <- mu_data$sqm/10000
sum(mu_data$hectares) #41,477,751 hectares

#read-in component data and define major component as >15% aerial coverage
list.files(SSURGOdir)
comp_data <- read.csv(file.path(SSURGOdir, 'ca_components.csv'), na.strings = c('', ' '), stringsAsFactors = FALSE)
length(unique(comp_data$cokey)) #90589 unique cokeys
length(unique(comp_data$mukey)) #18737 unique mukeys
if(sum(is.na(comp_data$majcompflag)) > 0) {stop(print('there are NAs in majcomp column!'))}
if(sum(is.na(comp_data$comppct_r)) > 0) {stop(print('there are NAs in the comppct column!'))}
sum(comp_data$comppct_r[comp_data$majcompflag=='Yes'] < 15) #143 are major components with <15% aereal coverage
sum(comp_data$comppct_r[comp_data$majcompflag=='No '] >= 15) #277 are minor components with >15% aereal coverage
comp_data$majcompflag[comp_data$majcompflag=='No ' & comp_data$comppct_r>=15] <- 'Yes'
comp_data$majcompflag[comp_data$majcompflag=='Yes' & comp_data$comppct_r < 15] <- 'No '

#read-in soil horizon data and convert all 0% SOM reporting to NA
horizon_data <- read.csv(file.path(SSURGOdir, 'ca_horizons.csv'), na.strings = c('', ' '), stringsAsFactors = FALSE)
colnames(horizon_data)
length(unique(horizon_data$cokey)) #40357 unique cokeys
sum(horizon_data$om_r==0, na.rm = TRUE) #4206 instances of 0% SOM
horizon_data$om_r[horizon_data$om_r==0] <- NA
horizon_data$majcompflag <- comp_data$majcompflag[match(horizon_data$cokey, comp_data$cokey)]
table(horizon_data$majcompflag)
horizon_data_majcomps <- horizon_data[horizon_data$majcompflag=='Yes',]
horizon_data_majcomps$mukey <- comp_data$mukey[match(horizon_data_majcomps$cokey, comp_data$cokey)]

#convert data.frame to a SoilProfileCollection object
#see https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/aqp-intro.html?root=aqp
horizon_data_spc <- horizon_data_majcomps
depths(horizon_data_spc) <- cokey ~ hzdept_r + hzdepb_r
class(horizon_data_spc)
print(horizon_data_spc)
depth_logic_result <- checkHzDepthLogic(horizon_data_spc)
lapply(depth_logic_result[,2:6], table)
if(!all(depth_logic_result$valid)) {stop(print('There are errors in SSURGO horizonation that need to be fixed.'))}

#functions to work with ProfileApply
wtd.mean <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzdepb_r - x$hzdept_r
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=FALSE)
  m
}
horizon_to_comp <- function(horizon_SPC, depth, comp_df, vars_of_interest = 'om_r', varnames = 'SOM', SOC_content=FALSE) {
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

comp_SOM_30cm <- horizon_to_comp(horizon_SPC = horizon_data_spc, depth = 30, comp_df = comp_data)
summary(comp_SOM_30cm$SOM_30cm)
hist(comp_SOM_30cm$SOM_30cm)
write.csv(comp_SOM_30cm, file.path(SSURGOdir, 'intermediate results', 'comp_SOM_30cm.csv'), row.names = FALSE)

comppct_by_mukey_om_data <- data.frame(mukey=row.names(tapply(comp_SOM_30cm$comppct[!is.na(comp_SOM_30cm$SOM_30cm)], comp_SOM_30cm$mukey[!is.na(comp_SOM_30cm$SOM_30cm)], sum)), comppct_tot=as.numeric(tapply(comp_SOM_30cm$comppct[!is.na(comp_SOM_30cm$SOM_30cm)], comp_SOM_30cm$mukey[!is.na(comp_SOM_30cm$SOM_30cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_om_data$comppct_tot)
sum(comppct_by_mukey_om_data$comppct_tot < 70) #2505 mukeys less than 70
sum(comppct_by_mukey_om_data$comppct_tot==15) #70 have only 15% of components with data
sum(comppct_by_mukey_om_data$comppct_tot < 15) #0, because major component defined as >15% of mapunit

CA_30cm_muagg <- MUAggregate_wrapper(df1=comp_SOM_30cm, varnames = 'SOM_30cm')
CA_30cm_muagg$compct_SOM <- comppct_by_mukey_om_data$comppct_tot[match(CA_30cm_muagg$mukey, comppct_by_mukey_om_data$mukey)]
CA_30cm_muagg$hectares <- mu_data$hectares[match(CA_30cm_muagg$mukey, mu_data$mukey)]
summary(CA_30cm_muagg$hectares)
sum(CA_30cm_muagg$hectares[which(CA_30cm_muagg$compct_SOM==15.00)])
sum(CA_30cm_muagg$hectares) #32,470,962 total hectares
sum(CA_30cm_muagg$hectares[!is.na(CA_30cm_muagg$SOM_30cm)]) #29,652,204 hectares with data

#add result to file with all mukeys
mu_data$SOM_30cm <- CA_30cm_muagg$SOM_30cm[match(mu_data$mukey, CA_30cm_muagg$mukey)]
mu_data$compct_SOM <- CA_30cm_muagg$compct_SOM[match(mu_data$mukey, CA_30cm_muagg$mukey)]
sum(mu_data$hectares[!is.na(mu_data$SOM_30cm)])/sum(mu_data$hectares) #71.5% of AOI with SOM summaries
length(mu_data$mukey[!is.na(mu_data$SOM_30cm)])/length(mu_data$mukey) #91% of mukeys
write.csv(mu_data, file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_30cm.csv'), row.names = FALSE)

#manual check of calculations
mu_data[mu_data$mukey==455490,]
comp_data[comp_data$mukey==455490,]
horizon_data_majcomps[horizon_data_majcomps$cokey==19535419,]
0.5*(20/30)+0.25*(10/30) #0.4166667 % SOM
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
head(SOM_wtd_data)
length(SOM_wtd_data) #3247044
writeClipboard(as.character(quantile(SOM_wtd_data, na.rm = TRUE, probs = c(0.25, 0.5, 0.75))))
#25%  50%  75% 
#0.75 1.50 2.50 
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.2, 0.4, 0.6, 0.8), na.rm = TRUE)))
#20%       40%       60%       80% 
#0.7433333 1.1911765 1.7600000 2.7333333
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.167, 0.334, 0.501, 0.668, 0.835), na.rm = TRUE)))
#10%       25%       50%       75%       90% 
#0.4960784 0.7500000 1.5000000 2.5000000 4.5952381 
writeClipboard(as.character(quantile(SOM_wtd_data, probs = c(0.143, 0.286, 0.429, 0.572, 0.715, 0.858), na.rm = TRUE)))
#10%       20%       40%       60%       80%       90% 
#0.4960784 0.7433333 1.1911765 1.7600000 2.7333333 4.5952381

#manual check of spatially weighted 40th percentile
sum(mu_data$hectares[which(mu_data$SOM_30cm < 1.1911765)]) #11,866,719 hectares
11866719/sum(mu_data$hectares[!is.na(mu_data$SOM_30cm)])
#40% of area of interest with SOM data < 1.191%