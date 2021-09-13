workDir <- 'C:/Users/smdevine/Desktop/PostDoc/SoilC/SSURGO'
library(raster)
list.files(workDir)
OM_to_OC <- 1.72 #OM = OC * 1.72
components_CA <- read.csv(file.path(workDir, 'ca_component_data.csv'), stringsAsFactors = FALSE)
dim(components_CA) #92,807 entries
length(unique(components_CA$cokey)) #91,193 unique cokeys
length(unique(components_CA$cokey[components_CA$majcompflag=='Yes']))
unique(components_CA$majcompflag)
unique(components_CA$reskind)
sum(components_CA$reskind=="")
components_CA$reskind[components_CA$reskind==""] <- 'None'
components_restrictions <- data.frame(cokey=unique(components_CA$cokey), reskinds=as.character(tapply(components_CA$reskind, components_CA$cokey, unique)), reskind_n=as.integer(tapply(components_CA$reskind, components_CA$cokey, function(x) length(x))))
dim(components_restrictions)
summary(components_restrictions$reskind_n)
sum(components_restrictions$reskind_n==3)
sum(components_restrictions$reskind_n==2)
sum(components_restrictions$reskind_n==0)
summary(components_restrictions$reskinds)
components_restrictions$reskinds

components_paralithic <- components_restrictions[grepl('Paralithic bedrock', components_restrictions$reskinds), ]
dim(components_paralithic)
length(unique(components_paralithic$cokey))
summary(components_paralithic$reskinds)
components_lithic <- components_restrictions[grepl('Lithic bedrock', components_restrictions$reskinds), ]
dim(components_lithic)
summary(components_lithic)
cokeys_rock <- unique(c(components_lithic$cokey, components_paralithic$cokey))
length(cokeys_rock)

horizons_CA <- read.csv(file.path(workDir, 'ca_horizon_data.csv'), stringsAsFactors = FALSE)
dim(horizons_CA)
summary(horizons_CA) #17 horizons have fragment vol > 100%
horizons_CA$hz_thickness <- horizons_CA$hzdepb_r - horizons_CA$hzdept_r
horizons_CA$kg_OCm2 <- (horizons_CA$om_r / OM_to_OC) * horizons_CA$dbthirdbar_r * (100 - ifelse(horizons_CA$fragvol_r_sum > 100, NA, horizons_CA$fragvol_r_sum)) * horizons_CA$hz_thickness / 1000
summary(horizons_CA$kg_OCm2) #39630

components_OC <- data.frame(cokey=unique(horizons_CA$cokey), kg_OCm2=as.numeric(tapply(horizons_CA$kg_OCm2, horizons_CA$cokey, sum, na.rm=TRUE)), pedon_depth=as.numeric(tapply(horizons_CA$hz_thickness, horizons_CA$cokey, sum, na.rm=TRUE)))
dim(components_OC)
summary(components_OC$kg_OCm2)
summary(components_OC$pedon_depth)
components_OC[components_OC$pedon_depth==999,]
sum(components_OC$pedon_depth < 30)
hist(components_OC$kg_OCm2[components_OC$kg_OCm2 < 100])
hist(components_OC$pedon_depth)

#divide horizon data into those with lithic or paralithic contacts or not
i <- horizons_CA$cokey %in% cokeys_rock
horizons_CA_Cr_R <- horizons_CA[i, ]
horizons_CA_no_Cr_R <- horizons_CA[!i, ]

#now sum OC by different depth intervals
sum_kg_OC_m2 <- function(df, depth, OC) {
  ifelse(is.na(df$hz_thickness) & is.na(df[[OC]]), NA,
         ifelse(df$hzdepb_r <= depth, df[[OC]],
                ifelse(df$hzdept_r < depth, 
                       ((depth - df$hzdept_r)/df$hz_thickness) * df[[OC]], NA)))
}
horizons_CA$kg_OCm2_30cm <- sum_kg_OC_m2(horizons_CA, 30, 'kg_OCm2')
summary(horizons_CA$kg_OCm2_30cm)
sum_OC_soilthickness <- function(df, depth, OC) { 
  ifelse(df$hzdepb_r <= depth & (!is.na(df[[OC]]) & df[[OC]] != 0), 
         df$hz_thickness,
         ifelse(df$hzdept_r < depth & (!is.na(df[[OC]]) & df[[OC]] != 0),
                depth - df$hzdept_r, 0))
}
horizons_CA$OC_30cm_thickness <- sum_OC_soilthickness(horizons_CA, 30, 'kg_OCm2')
sum_OC_bycokey <- function(df, var) {
  y <- aggregate(x=df[[var]], by=list(df$cokey), FUN=sum_modified)
  colnames(y) <- c('cokey', paste0(var, '_comp'))
  return(y)
}
sum_modified <- function(x) {
  if(all(is.na(x))) {
    return(NA)
  }
  else {sum(x, na.rm = TRUE)}
}
comps_OC_unmodified <- Reduce(function(...) merge(..., all=TRUE), list(sum_OC_bycokey(horizons_CA, 'kg_OCm2_30cm'), sum_OC_bycokey(horizons_CA, 'OC_30cm_thickness')))
summary(comps_OC_unmodified)
dim(comps_OC_unmodified)
