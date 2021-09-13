#TO-DO:
#re-make violin plots; remove point with 0-30 cm horizonation problems
#re-ran 8/24/21 with crops trimmmed to "speciality" only (no field crops or alfalfa)
laptop <- FALSE
library(vioplot)
library(raster)
library(aqp)
library(soilDB)
library(lattice)
# library(corrplot)
# library(cluster)
# library(factoextra)
# library(fpc)
# library(fmsb)
library(extrafont)
library(extrafontdb)
om_to_oc <- 1.72
crit_pH <- 7.8
clus_7_names <- c('6. Fine salt-affected', '3. Low OM with restrictive horizons', '4. High OM with restrictive horizons', '1. Coarse with no restrictions', '2. Loamy with no restrictions', '7. Shrink-swell', '5. Coarse-loamy salt-affected')
clus_7_colors <- c('deepskyblue', 'olivedrab3', 'firebrick3', 'lightgoldenrod', 'tan4', 'violetred', 'lightblue1') #olivedrab3
order_lgnd_7 <- c(4,5,2,3,7,1,6)
# font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
#have to install older version of package that extrafont uses to get fonts to load correctly
# library(remotes)
# install_version("Rttf2pt1", version = "1.3.8")
loadfonts(device = 'win')
if (laptop) {
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
} else { #on UCD desktop
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
}
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/Figures' #was valley_final
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/SoilC/Figures' #was valley_final
  LandIQDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/crops'
  ResultsDir <- 'C:/Users/smdevine/Desktop/PostDoc/SoilC/Crops'
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}

crops <- shapefile(file.path(LandIQDir, 'i15_Crop_Mapping_2014_Final_LandIQonAtlas.shp'))
crs(crops)
unique(crops$Crop2014)
sum(crops$Acres) #14202468
sum(area(crops) / 10000 * 2.47105) #14202446
crops_simplifed <- c('Alfalfa and pasture', 'Almonds', 'Deciduous orchard', 'Deciduous orchard', 'Legumes', 'Small fruit and vegetables', 'Small fruit and vegetables', 'Deciduous orchard', 'Citrus', 'Small fruit and vegetables', 'Cereals and grain hay', 'Cotton', 'Specialty crops', 'Specialty crops', 'Grapes', 'Specialty crops', 'Idle', 'Specialty crops', 'Small fruit and vegetables', 'Managed wetland', 'Small fruit and vegetables', 'Deciduous orchard', 'Cereals and grain hay', 'Cereals and grain hay', 'Alfalfa and pasture', 'Specialty crops', 'Small fruit and vegetables', 'Alfalfa and pasture', 'Olives', 'Small fruit and vegetables', 'Deciduous orchard', 'Deciduous orchard', 'Small fruit and vegetables', 'Pistachios', 'Deciduous orchard', 'Deciduous orchard', 'Small fruit and vegetables', 'Rice', 'Oil crops', 'Small fruit and vegetables', 'Oil crops', 'Tomatoes', 'Urban', 'Walnuts', 'Cereals and grain hay', 'Rice', 'Young perennials')
crops_simplifed2 <- c('Alfalfa and pasture', 'Orchard', 'Orchard', 'Orchard', 'Annual field crops', 'Small fruit and vegetables', 'Small fruit and vegetables', 'Orchard', 'Orchard', 'Small fruit and vegetables', 'Annual field crops', 'Annual field crops', 'Specialty crops', 'Specialty crops', 'Grapes', 'Specialty crops', 'Idle', 'Specialty crops', 'Small fruit and vegetables', 'Managed wetland', 'Small fruit and vegetables', 'Orchard', 'Annual field crops', 'Annual field crops', 'Alfalfa and pasture', 'Specialty crops', 'Small fruit and vegetables', 'Alfalfa and pasture', 'Orchard', 'Small fruit and vegetables', 'Orchard', 'Orchard', 'Small fruit and vegetables', 'Orchard', 'Orchard', 'Orchard', 'Small fruit and vegetables', 'Annual field crops', 'Annual field crops', 'Small fruit and vegetables', 'Annual field crops', 'Annual field crops', 'Urban', 'Orchard', 'Annual field crops', 'Annual field crops', 'Young perennials')
length(crops_simplifed)
length(crops_simplifed2)
unique(crops_simplifed) #20
unique(crops_simplifed2) #10
crop_cats <- data.frame(crop=unique(crops$Crop2014)[order(unique(crops$Crop2014))], crop_cat1=crops_simplifed, crop_cat2=crops_simplifed2, stringsAsFactors = FALSE)
head(crop_cats)
crops$crop_cat <- crop_cats$crop_cat2[match(crops$Crop2014, crop_cats$crop)]
unique(crops$crop_cat)
crops_specialty <- crops[crops$crop_cat %in% c('Orchard', 'Grapes', 'Small fruit and vegetables', 'Specialty crops', 'Young perennials'),]
sum(crops_specialty$Acres) #4198063
unique(crops_specialty$crop_cat)
shapefile(crops_specialty, file.path(LandIQDir, 'speciality crops only', 'crops_specialty.shp'))
sum(area(crops_specialty))
(16988959243/10000)*2.47105 #4,198,057 acres

#performed intersection between specialty crops and soils in ArcGIS Desktop 10.5 on 8/24/21 using Analysis:Overlay:Intersect
valley_crops_mu <- shapefile(file.path(ResultsDir, 'specialty_crop_SHRs.shp')) #539297 little polygons
names(valley_crops_mu)
valley_crops_mu$area_ac <- area(valley_crops_mu) / 10000 * 2.47105
sum(valley_crops_mu$area_ac) #3490917
sum(valley_crops_mu$area_ac) / sum(crops_specialty$Acres) #defined SHRs cover 83.2% of specialty crop acreage
length(unique(valley_crops_mu$Crop2014)) #30 crops
unique(valley_crops_mu$Crop2014)
valley_crops_mu$SHR7name <- clus_7_names[as.integer(valley_crops_mu$cluster_7)]
shapefile(valley_crops_mu, file.path(ResultsDir, 'specialty_crop_SHRs.shp'), overwrite=TRUE)

#get acreage summary by mukey
acres_by_mukey <- tapply(valley_crops_mu$area_ac, valley_crops_mu$mukey, sum)

#get data.frame listing all mukeys in specialty crop zone
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', 'valley30cm_by_mukey_cluster_FINAL.csv'), stringsAsFactors = FALSE)
valley30cm_by_mukey_specialty <- valley30cm_by_mukey[valley30cm_by_mukey$mukey %in% valley_crops_mu$mukey,]
dim(valley30cm_by_mukey_specialty)
sum(valley30cm_by_mukey_specialty$area_ac)
sum(valley_crops_mu$area_ac)
valley30cm_by_mukey_specialty$area_ac <- acres_by_mukey[match(valley30cm_by_mukey_specialty$mukey, names(acres_by_mukey))]

#violin plots
kssl_SHR_shp <- shapefile(file.path(ksslDir, 'shapefiles', 'kssl_SHR7.shp')) #369 points
kssl_SHR_shp <- spTransform(kssl_SHR_shp, crs(valley_crops_mu))
kssl_specialty_crop <- kssl_SHR_shp[!is.na(over(kssl_SHR_shp, valley_crops_mu)$cluster_7),]
kssl_specialty_crop #only 92 points left
kssl_specialty_crop$pedn_ky

#now read-in from creation below
kssl_specialty_crop <- shapefile(file.path(ResultsDir, 'points', 'kssl_specialty_crop.shp'))

kssl_points_30cm <- read.csv(file.path(ksslDir, 'kssl_cluster_30cm_FINAL.csv'), stringsAsFactors = FALSE) #oc and kgSOC_m2 updated 2/11/20
#replace OC with totC when pH sufficiently low
kssl_points_30cm$oc_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)] <- kssl_points_30cm$c_tot_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)]
sum(is.na(kssl_points_30cm$oc)) #NAs reduced from 88 to 56
#replace estimated om with oc times assumption (see 'aov_soil_properties_KSSL_CDFA.R' for details comparing KSSL OM est with this)
kssl_points_30cm$om_30cm <- kssl_points_30cm$oc_30cm * om_to_oc
kssl_points_30cm$xdim_vioplot7 <- match(kssl_points_30cm$cluster_7, order_lgnd_7)
kssl_points_30cm$lep_30cm <- kssl_points_30cm$lep_30cm*100 #to match SSURGO scale
# kssl_points_30cm$mukey <- kssl_ssurgo_extract$mukey[match(kssl_points_30cm$pedon_key, kssl_ssurgo_extract$pedon_key)]
# kssl_points_30cm$dom_order <- dom_order_by_mukey$dom_order[match(kssl_points_30cm$mukey, dom_order_by_mukey$mukey)]
kssl_specialty_crops_30cm <- kssl_points_30cm[kssl_points_30cm$pedon_key %in% kssl_specialty_crop$pedn_ky,]
dim(kssl_specialty_crops_30cm) #92
table(kssl_specialty_crop$SHR7name)
kssl_specialty_crop$om_30cm <- kssl_specialty_crops_30cm$om_30cm[match(kssl_specialty_crop$pedn_ky, kssl_specialty_crops_30cm$pedon_key)]
kssl_specialty_crop$oc_30cm <- kssl_specialty_crops_30cm$oc_30cm[match(kssl_specialty_crop$pedn_ky, kssl_specialty_crops_30cm$pedon_key)]
kssl_specialty_crops_30cm$SHR7name <- clus_7_names[kssl_specialty_crops_30cm$cluster_7]
shapefile(kssl_specialty_crop, file.path(ResultsDir, 'points', 'kssl_specialty_crop.shp'), overwrite=TRUE)
write.csv(kssl_specialty_crops_30cm, file.path(ResultsDir, 'points', 'kssl_specialty_crop_30cm.csv'), row.names = FALSE)

#write all KSSL pts in SHRs to shapefile with 30-cm data
all(kssl_SHR_shp$pedn_ky, kssl_points_30cm$pedon_key)
kssl_SHR_30cm_shp <- merge(kssl_SHR_shp, kssl_points_30cm[,1:25], by.x='pedn_ky', by.y='pedon_key')
kssl_SHR_30cm_shp$om_30cm
shapefile(kssl_SHR_30cm_shp, file.path(ResultsDir, 'points', 'kssl_SHR_30cm.shp'), overwrite=TRUE)
# kssl_SHR_shp[colnames(kssl_points_30cm)[2:25]] <- kssl_points_30cm[,2:25]


#now bring in CDFA points (aka "Kerri's" points)
soil_data <- read.csv(file.path(kerriDir, 'CDFA Soil Survey All Data_copy.csv'), stringsAsFactors = FALSE, colClasses = c(Concatenate='character'), na.strings = c('#N/A', 'pOOR DATA', 'POOR DATA', 'NO Sample', 'Missing Data', '?', ""))
dim(soil_data)
soil_data$GPS.N <- as.numeric(gsub("N", "", soil_data$GPS.N))
soil_data$GPS.W <- -as.numeric(gsub("W", "", soil_data$GPS.W))
soil_data_pts <- data.frame(ID=unique(soil_data$Concatenate),stringsAsFactors = FALSE)
soil_data_pts$Lat_WGS84 <- soil_data$GPS.N[match(soil_data_pts$ID, soil_data$Concatenate)]
soil_data_pts$Lon_WGS84 <- soil_data$GPS.W[match(soil_data_pts$ID, soil_data$Concatenate)]
soil_data_pts <- SpatialPointsDataFrame(coords = soil_data_pts[,c('Lon_WGS84', 'Lat_WGS84')], data = soil_data_pts['ID'], proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
soil_data_pts <- spTransform(soil_data_pts, crs(valley_crops_mu))
soil_data_pts$ID <- as.integer(soil_data_pts$ID)

kerri_specialty_crops <- soil_data_pts[!is.na(over(soil_data_pts, valley_crops_mu)$cluster_7), ]
kerri_specialty_crops #100 features
names(kerri_specialty_crops)

kerri_points_30cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
dim(kerri_points_30cm) #127 features
# kerri_points_30cm$xdim_vioplot9 <- match(kerri_points_30cm$cluster_9, order_lgnd_9)
kerri_points_30cm$xdim_vioplot7 <- match(kerri_points_30cm$cluster_7, order_lgnd_7)
kerri_points_30cm$SHR7name <- clus_7_names[kerri_points_30cm$cluster_7]
# kerri_points_30cm$xdim_vioplot6 <- match(kerri_points_30cm$cluster_6, order_lgnd_6)
# kerri_points_30cm$xdim_vioplot5 <- match(kerri_points_30cm$cluster_5, order_lgnd_5)
colnames(kerri_points_30cm)
sum(is.na(kerri_points_30cm$cluster_7)) #24 were already NA
kerri_points_30cm$om_30cm <- kerri_points_30cm$totC_30cm * om_to_oc
kerri_specialty_crops$om_30cm <- kerri_points_30cm$om_30cm[match(kerri_specialty_crops$ID, kerri_points_30cm$Concatenate)]
kerri_specialty_crops$oc_30cm <- kerri_points_30cm$totC_30cm[match(kerri_specialty_crops$ID, kerri_points_30cm$Concatenate)]
kerri_specialty_crops$SHR7name <- kerri_points_30cm$SHR7name[match(kerri_specialty_crops$ID, kerri_points_30cm$Concatenate)]
shapefile(kerri_specialty_crops, file.path(ResultsDir, 'points', 'kerri_specialty_crop.shp'), overwrite=TRUE)

kerri_specialty_crops_30cm <- kerri_points_30cm[kerri_points_30cm$Concatenate%in% kerri_specialty_crops$ID,]

#bind cdfa and kssl data together
shp1 <- kssl_specialty_crop[c('pedn_ky', 'om_30cm', 'oc_30cm', 'SHR7name')]
names(shp1)[1] <- 'ID'
shp1$source <- 'KSSL'
shp2 <- kerri_specialty_crops
shp2$source <- 'CDFA'
shp_all <- rbind(shp1, shp2)
hist(shp_all$om_30cm)
shapefile(shp_all, file.path(ResultsDir, 'points', 'all_pts_specialty_crop.shp'))

# kerri_ssurgo_extract <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_pts_ssurgo_30cm_extract_FINAL.csv'), stringsAsFactors = FALSE)
# kerri_points_30cm$mukey <- kerri_ssurgo_extract$mukey[match(kerri_points_30cm$Concatenate, kerri_ssurgo_extract$Concatenate)]

#sig_labels, fig_label
vioplot_mod_clus7_validation <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, kssl_df, kssl_varname, cdfa_pts, cdfa_varname, legend_plot, legendloc, legend_cex, legend_text, fig_height) {
  plot_order2 <- (1:7)[plot_order]
  tiff(file = file.path(FiguresDir, 'violin plots', fname), family = 'Times New Roman', width = 4.25, height = fig_height, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_7==plot_order2[1]], times=round(df$area_ac[df$cluster_7==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[2]], times=round(df$area_ac[df$cluster_7==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[3]], times=round(df$area_ac[df$cluster_7==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[4]], times=round(df$area_ac[df$cluster_7==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[5]], times=round(df$area_ac[df$cluster_7==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[6]], times=round(df$area_ac[df$cluster_7==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[7]], times=round(df$area_ac[df$cluster_7==plot_order2[7]]/area_fact, 0)), col=clus_7_colors[plot_order], wex=1.2, cex=0.8, rectCol = 'gray', ylim = ylim_vioplot, ylab = NULL)
  mtext('Soil health region', side = 1, line = 2)
  mtext(ylab, side = 2, line = 2)
  points(x=kssl_df$xdim_vioplot7[!is.na(kssl_df[[kssl_varname]])]-0.15, y=kssl_df[[kssl_varname]][!is.na(kssl_df[[kssl_varname]])]*if(kssl_varname=='totC_30cm'){1.72} else{1}, cex=0.5, pch=1, col='black')
  kssl_means <- data.frame(mean=tapply(kssl_df[[kssl_varname]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, kssl_df$cluster_7, mean, na.rm=TRUE))
  kssl_means$xdim_vioplot <- match(row.names(kssl_means), plot_order)
  points(x=kssl_means$xdim_vioplot-0.15, y=kssl_means$mean, pch=8, cex=0.5, col='orange')
  points(x=cdfa_pts$xdim_vioplot7[!is.na(cdfa_pts[[cdfa_varname]])]+0.15, y=cdfa_pts[[cdfa_varname]][!is.na(cdfa_pts[[cdfa_varname]])], cex=0.5, pch=4, col='black')
  cdfa_means <- data.frame(mean=tapply(cdfa_pts[[cdfa_varname]], cdfa_pts$cluster_7, mean, na.rm=TRUE))
  cdfa_means$xdim_vioplot <- match(row.names(cdfa_means), plot_order)
  points(x=cdfa_means$xdim_vioplot+0.15, y=cdfa_means$mean, pch=8, cex=0.5, col='darkblue')
  # text(x=1:7, y=ylim_vioplot[1], labels = sig_labels, adj=0.5, cex=0.9)
  # legend('topright', fig_label, bty='n', inset=0.005)
  if(legend_plot) {
    legend(x=legendloc, legend=c('SSURGO violin plots', 'KSSL data', 'KSSL mean', 'Napa-Lodi data', 'Napa-Lodi mean'), pch=c(NA,1,8,4,8), col = c(NA, 'black', 'orange', 'black', 'darkblue'), pt.cex=legend_cex, bty='n')
  }
  dev.off()
}

vioplot_mod_clus7_validation(valley30cm_by_mukey_specialty, 'om_30cm', ylim_vioplot = c(0,10), plot_order = order_lgnd_7, area_fact = 10, ylab='Organic matter (%)', fname='specialty_crop_OM_vioplots_KSSL_CDFA.tif', mar=c(0.02, 3.25, 0.25, 0.25), kssl_df = kssl_specialty_crops_30cm, cdfa_pts=kerri_specialty_crops_30cm, legend_plot=TRUE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "om_30cm", cdfa_varname = 'om_30cm', fig_height = 3.75) #sig_labels = c('A', 'C', 'A', 'C', 'A', 'AB', 'BC'), fig_label = 'b'

#now make SOM profile plot for KSSL points within specialty crop zone
kssl_horizons_SHR <- read.csv(file.path(ksslDir, 'kssl_horizons_SHRonly.csv'), stringsAsFactors = FALSE) #this had 0-30 cm bad horizonation fixed in kssl_validation_FINAL.R
colnames(kssl_horizons_SHR)
summary(kssl_horizons_SHR$oc) #497 NAs
summary(kssl_horizons_SHR$oc_est)
kssl_horizons_SHR$SOM <- kssl_horizons_SHR$oc_est*om_to_oc
summary(kssl_horizons_SHR$SOM)
kssl_horizons_SHR$SHR7name <- clus_7_names[kssl_horizons_SHR$SHR7code]
kssl_horizons_SHR$SHR7name <- as.factor(kssl_horizons_SHR$SHR7name)
kssl_specialty_crop_horizons_SHR <- kssl_horizons_SHR[kssl_horizons_SHR$pedon_key %in% kssl_specialty_crop$pedn_ky, ]
length(unique(kssl_specialty_crop_horizons_SHR$pedon_key)) #92
kssl_specialty_crop_horizons_SHR[kssl_specialty_crop_horizons_SHR$SHR7name=='6. Fine salt-affected',]
#get rid of horizonation problems identified below
horizons_to_exclude <- c('16N03041', 'UCD00571', 'UCD03484', '40A23760', '40A23759', '40A23007', '91P04452', '91P04453', '91P04454', '91P04451', '91P02049', '91P02050') 
pedons_to_exclude <- c(17781, 17785, 17778, 17772, 15471, 15468, 15467, 14403, 10997)
kssl_specialty_crop_horizons_SHR <- kssl_specialty_crop_horizons_SHR[!(kssl_specialty_crop_horizons_SHR$labsampnum %in% horizons_to_exclude), ]
kssl_specialty_crop_horizons_SHR <- kssl_specialty_crop_horizons_SHR[!(kssl_specialty_crop_horizons_SHR$pedon_key %in% pedons_to_exclude),]
kssl_specialty_crop_horizons_SHR$hzn_bot[kssl_specialty_crop_horizons_SHR$labsampnum=='40A23003'] <- 110
kssl_specialty_crop_horizons_SHR$hzn_bot[kssl_specialty_crop_horizons_SHR$labsampnum=='40A23356'] <- 56
kssl_specialty_crop_horizons_SHR$hzn_bot[kssl_specialty_crop_horizons_SHR$labsampnum=='83P01353'] <- 56
kssl_specialty_crop_horizons_SHR$hzn_bot[kssl_specialty_crop_horizons_SHR$labsampnum=='79P00464'] <- 25

depths(kssl_specialty_crop_horizons_SHR) <- pedon_key ~ hzn_top + hzn_bot
site(kssl_specialty_crop_horizons_SHR) <- ~ SHR7name
specialty.hz.ck <- checkHzDepthLogic(kssl_specialty_crop_horizons_SHR)
sum(specialty.hz.ck$valid) #ony 74 valid out of 92
lapply(specialty.hz.ck$pedon_key[specialty.hz.ck$valid==FALSE], function(x) print(kssl_horizons_SHR[kssl_horizons_SHR$pedon_key==x,]))

all.hz.ck <- checkHzDepthLogic() #check all horizonation
#change 40A23003      hzn_bot 91 to 110 
#change 40A23356      hzn_bot 53 to 56
#change 83P01353      hzn_bot 55 to 56     

kssl_specialty_crop_horizons_SHR.slab <- slab(kssl_specialty_crop_horizons_SHR, SHR7name ~ SOM, slab.structure = 1)
dim(kssl_specialty_crop_horizons_SHR.slab)
str(kssl_specialty_crop_horizons_SHR.slab)
levels(kssl_specialty_crop_horizons_SHR.slab$variable)
levels(kssl_specialty_crop_horizons_SHR.slab$variable) <- 'Soil organic matter (%)'
tps <- list(superpose.line=list(col=clus_7_colors[order(clus_7_names)], lwd=2))
class(kssl_specialty_crop_horizons_SHR.slab)
kssl_specialty_crop_horizons_SHR.slab <- kssl_specialty_crop_horizons_SHR.slab[kssl_specialty_crop_horizons_SHR.slab$bottom<=150,]

tiff(file = file.path(FiguresDir, 'SOM depth plots', 'SOM_KSSL_profiles_specialty_crop.tif'), family = 'Times New Roman', width = 9, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
xyplot(top ~ p.q50 | variable, groups=SHR7name, data=kssl_specialty_crop_horizons_SHR.slab, ylab='Depth (cm)', xlab='specialty crop zone KSSL median bounded by 25th and 75th percentiles', lower=kssl_specialty_crop_horizons_SHR.slab$p.q25, upper=kssl_specialty_crop_horizons_SHR.slab$p.q75, ylim=c(155,-5), xlim = c(0,4), panel=panel.depth_function, alpha=0.4, sync.colors=TRUE, prepanel=prepanel.depth_function, par.strip.text=list(cex=0.8), strip=strip.custom(bg=grey(0.85)), layout=c(1,1), scales=list(x=list(alternating=1, relation='free'), y=list(alternating=3)), par.settings=tps, auto.key=list(columns=3, lines=TRUE, points=FALSE, lwd=2))
dev.off()
kssl_specialty_crop_horizons_SHR.slab[kssl_specialty_crop_horizons_SHR.slab$SHR7name=='1. Coarse with no restrictions',]

tapply(kssl_specialty_crops_30cm$om_30cm, kssl_specialty_crops_30cm$SHR7name, summary)
tapply(kssl_specialty_crops_30cm$om_30cm, kssl_specialty_crops_30cm$SHR7name, function(x) sum(!is.na(x)))
clus_7_colors[order(clus_7_names)]
