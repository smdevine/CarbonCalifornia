#revised 11/21/19 because coords are WGS84 not NAD83
workDir <- 'C:/Users/smdevine/Desktop/PostDoc/SoilC/SSURGO'
library(raster)
list.files(workDir)
kssl_sites <- read.csv(file.path(workDir, 'ca_kssl_sites.csv'), stringsAsFactors = FALSE)
dim(kssl_sites)
head(kssl_sites)
colnames(kssl_sites)
kssl_sites_shp <- SpatialPointsDataFrame(coords = kssl_sites[,c('lon', 'lat')], data = kssl_sites[,1:13], proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
+towgs84=0,0,0")))
plot(kssl_sites_shp)
shapefile(kssl_sites_shp, file.path(workDir, 'shapefiles', 'kssl_CA.shp'), overwrite=TRUE)
kssl_horizons <- read.csv(file.path(workDir, 'ca_kssl_horizons.csv'), stringsAsFactors = FALSE)
dim(kssl_horizons)
colnames(kssl_horizons)
summary(kssl_horizons$oc)
summary(kssl_horizons$db_13b)
length(unique(kssl_horizons$pedon_key))
length(unique(kssl_sites$pedon_id))
kssl_horizons$hz_thickness <- kssl_horizons$hzn_bot - kssl_horizons$hzn_top
summary(kssl_horizons$hz_thickness)
hist(kssl_horizons$hz_thickness)
kssl_horizons$kgOC_m2 <- kssl_horizons$oc * kssl_horizons$db_13b * (100 - kssl_horizons$frags) * kssl_horizons$hz_thickness / 1000
summary(kssl_horizons$kgOC_m2)
sum(!is.na(kssl_horizons$db_13b))
sum(grepl('UCD', kssl_horizons$labsampnum)) #4353 UCD horizons, all missing BD
ucdlabnumbers <- kssl_horizons$labsampnum[grepl('UCD', kssl_horizons$labsampnum)]
kssl_horizons_UCD <- kssl_horizons[grepl('UCD', kssl_horizons$labsampnum),]
summary(kssl_horizons_UCD$db_13b)
length(unique(kssl_horizons_UCD$pedon_key)) #782 pedons

#examples of how merge can be accomplished to link horizons to sites
kssl_horizons[kssl_horizons$labsampnum == ucdlabnumbers[1],]
kssl_horizons[kssl_horizons$pedon_key==52474,]
kssl_sites[kssl_sites$pedon_key==52474,]

kssl_sites[kssl_sites$pedon_id=='56-CA-11-021',]
kssl_horizons[kssl_horizons$pedon_key==52560,]

kssl_sites[kssl_sites$pedon_id=='72-CA-04-041x',]
kssl_horizons[kssl_horizons$pedon_key==52952,]

kssl_sites[kssl_sites$pedon_id=='62-CA-12-028x',]
kssl_horizons[kssl_horizons$pedon_key==52713,]

kssl_sites[kssl_sites$pedon_id=='61-CA-45-111',]
kssl_horizons[kssl_horizons$pedon_key==52690,]

kssl_sites[kssl_sites$pedon_id=='S-77-CA-47-136x',]
kssl_horizons[kssl_horizons$pedon_key==53243,]

#horizon data with DB
kssl_horizons_withDB <- kssl_horizons[!is.na(kssl_horizons$db_13b),]
length(unique(kssl_horizons_withDB$pedon_key)) #941 pedons have at least 1 hz with DB
kssl_sites_withDB <- kssl_sites[kssl_sites$pedon_key %in% unique(kssl_horizons_withDB$pedon_key),]
dim(kssl_sites_withDB) #941 indeed
kssl_sites_withDB_shp <- SpatialPointsDataFrame(coords = kssl_sites_withDB[,c('lon', 'lat')], data = kssl_sites_withDB[,1:13], proj4string = CRS(as.character("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")))
#shapefile(kssl_sites_withDB_shp, file.path(workDir, 'shapefiles', 'withDB.shp'))
plot(kssl_sites_withDB_shp)

#data from UCD Soil-Veg
kssl_sites_UCD <- kssl_sites[kssl_sites$pedon_key %in% unique(kssl_horizons_UCD$pedon_key),]
dim(kssl_sites_UCD)
kssl_sites_UCD_shp <- SpatialPointsDataFrame(coords = kssl_sites_UCD[,c('lon', 'lat')], data = kssl_sites_UCD[,1:13], proj4string = CRS(as.character("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")))
#shapefile(kssl_sites_UCD_shp, file.path(workDir, 'shapefiles', 'UCDsites.shp'))
plot(kssl_sites_UCD_shp, col='red', add=TRUE)

#check SOC data
pedons_with_data <- unique(kssl_horizons$pedon_key[!is.na(kssl_horizons$kgOC_m2)]) #622 have at least one horizon with SOC stock
kssl_horizons_SOC <- kssl_horizons[kssl_horizons$pedon_key %in% pedons_with_data,]
dim(kssl_horizons_SOC) #4205 horizons
length(unique(kssl_horizons_SOC$pedon_key)) #622 profiles
summary(kssl_horizons_SOC$kgOC_m2) #1,373 horizons are NA
summary(kssl_horizons_SOC)
head(kssl_horizons_SOC, 20)
#sort to verify
kssl_horizons_SOC <- kssl_horizons_SOC[order()]
tapply(kssl_horizons_SOC)
profile_SOC <- data.frame(pedon_key=unique(kssl_horizons_SOC$pedon_key), kgOC_m2 = as.numeric(tapply(kssl_horizons_SOC$kgOC_m2, kssl_horizons_SOC$pedon_key, sum, na.rm=TRUE)), SOC_thickenss=as.numeric(tapply(kssl_horizons_SOC$hz_thickness[!is.na(kssl_horizons_SOC$kgOC_m2)], kssl_horizons_SOC$pedon_key[!is.na(kssl_horizons_SOC$kgOC_m2)], sum)), profile_thickness=as.numeric(tapply(kssl_horizons_SOC$hz_thickness, kssl_horizons_SOC$pedon_key, sum)))
dim(profile_SOC) #622 profiles
summary(profile_SOC)
sum(profile_SOC$kgOC_m2 > 0) #9 are 0
sum(profile_SOC$SOC_thickenss==0)
sum(profile_SOC$SOC_thickenss == profile_SOC$profile_thickness) #163 of 622 profiles have complete SOC profiles
kssl_sites_SOC <- kssl_sites[kssl_sites$pedon_key %in% unique(kssl_horizons_SOC$pedon_key),]
kssl_sites_SOC$kgOC_m2 <- round(profile_SOC$kgOC_m2[match(kssl_sites_SOC$pedon_key, profile_SOC$pedon_key)], digits=2)
kssl_sites_SOC$SOC_thickeness <- profile_SOC$SOC_thickenss[match(kssl_sites_SOC$pedon_key, profile_SOC$pedon_key)]
kssl_sites_SOC_shp <- SpatialPointsDataFrame(coords = kssl_sites_SOC[,c('lon', 'lat')], data = kssl_sites_SOC[,c(1:13, 16:17)], proj4string = CRS(as.character("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")))

shapefile(kssl_sites_SOC_shp, file.path(workDir, 'shapefiles', 'SOCsites.shp'), overwrite=TRUE)

#couple of checks on
summary(kssl_sites_SOC$kgOC_m2)
kssl_horizons[kssl_horizons$pedon_key==13649,] #15.89 kg OC m^-2 is correct
sum(kssl_sites_SOC$kgOC_m2 < 40) #621 of 622 are less than 40 kg OC m-2
hist(kssl_sites_SOC$kgOC_m2[kssl_sites_SOC$kgOC_m2 < 40])