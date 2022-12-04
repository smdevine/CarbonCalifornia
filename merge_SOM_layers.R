mainDir <- 'C:/Users/smdevine/Desktop/post doc/soil organic carbon'
SSURGOdir <- file.path(mainDir, 'SSURGO data')
list.files(file.path(SSURGOdir, 'summaries'))
summary_100cm <- read.csv(file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_100cm_2022-10-11.csv'), stringsAsFactors = FALSE)
summary_30cm <- read.csv(file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_30cm_2022-10-11.csv'), stringsAsFactors = FALSE)
summary_30_50cm <- read.csv(file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_30_50cm_2022-10-11.csv'), stringsAsFactors = FALSE)
summary_50_75cm <- read.csv(file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_50_75cm_2022-10-11.csv'), stringsAsFactors = FALSE)
summary_75_100cm <- read.csv(file.path(SSURGOdir, 'summaries', 'CA_mapunit_SOM_75_100cm_2022-10-11.csv'), stringsAsFactors = FALSE)
dim(summary_100cm) #19036     7
lapply(summary_100cm[,4:7], summary)
lapply(summary_30cm[,4:7], summary)
all(summary_100cm$mukey==summary_30cm$mukey)
all(summary_30cm$mukey==summary_30_50cm$mukey)
all(summary_30cm$mukey==summary_50_75cm$mukey)
all(summary_30cm$mukey==summary_75_100cm$mukey)
sum(summary_100cm$kgSOC.m2_100cm < summary_30cm$kgSOC.m2_30cm, na.rm = TRUE) #32

sum(summary_30cm$compct_SOM_30cm==0, na.rm = TRUE)
summary_30cm$compct_SOM_30cm[is.na(summary_30cm$compct_SOM_30cm)] <- 0
summary_30cm$compct_SOCkg_30cm[is.na(summary_30cm$compct_SOCkg_30cm)] <- 0
summary(summary_30cm$SOM_30cm)
summary_30cm$SOM_30cm <- round(summary_30cm$SOM_30cm, 3)
summary_30cm$kgSOC.m2_30cm <- round(summary_30cm$kgSOC.m2_30cm, 3)

summary_30_50cm$compct_SOM_30_50cm[is.na(summary_30_50cm$compct_SOM_30_50cm)] <- 0
summary_30_50cm$compct_SOCkg_30_50cm[is.na(summary_30_50cm$compct_SOCkg_30_50cm)] <- 0
summary_30_50cm$SOM_30_50cm <- round(summary_30_50cm$SOM_30_50cm, 3)
summary_30_50cm$kgSOC.m2_30_50cm <- round(summary_30_50cm$kgSOC.m2_30_50cm, 3)

summary_50_75cm$compct_SOM_50_75cm[is.na(summary_50_75cm$compct_SOM_50_75cm)] <- 0
summary_50_75cm$compct_SOCkg_50_75cm[is.na(summary_50_75cm$compct_SOCkg_50_75cm)] <- 0
summary_50_75cm$SOM_50_75cm <- round(summary_50_75cm$SOM_50_75cm, 3)
summary_50_75cm$kgSOC.m2_50_75cm <- round(summary_50_75cm$kgSOC.m2_50_75cm, 3)

summary_75_100cm$compct_SOM_75_100cm[is.na(summary_75_100cm$compct_SOM_75_100cm)] <- 0
summary_75_100cm$compct_SOCkg_75_100cm[is.na(summary_75_100cm$compct_SOCkg_75_100cm)] <- 0
summary_75_100cm$SOM_75_100cm <- round(summary_75_100cm$SOM_75_100cm, 3)
summary_75_100cm$kgSOC.m2_75_100cm <- round(summary_75_100cm$kgSOC.m2_75_100cm, 3)

summary_100cm$compct_SOM_100cm[is.na(summary_100cm$compct_SOM_100cm)] <- 0
summary_100cm$compct_SOCkg_100cm[is.na(summary_100cm$compct_SOCkg_100cm)] <- 0
summary_100cm$SOM_100cm <- round(summary_100cm$SOM_100cm, 3)
summary_100cm$kgSOC.m2_100cm <- round(summary_100cm$kgSOC.m2_100cm, 3)

summary_all <- cbind(summary_30cm, summary_30_50cm[4:7], summary_50_75cm[4:7], summary_75_100cm[4:7], summary_75_100cm[4:7], summary_100cm[4:7])
lapply(summary_all, class)

sum(summary_100cm$compct_SOM_100cm==summary_30cm$compct_SOM_30cm & summary_100cm$compct_SOM_100cm!=0) #16520
sum(summary_30cm$compct_SOM_30cm==summary_30_50cm$compct_SOM_30_50cm & summary_30cm$compct_SOM_30cm==summary_50_75cm$compct_SOM_50_75cm & summary_30cm$compct_SOM_30cm==summary_75_100cm$compct_SOM_75_100cm & summary_100cm$compct_SOM_100cm!=0) #10679

head(summary_all[which(summary_30cm$compct_SOM_30cm==summary_30_50cm$compct_SOM_30_50cm & summary_30cm$compct_SOM_30cm==summary_50_75cm$compct_SOM_50_75cm & summary_30cm$compct_SOM_30cm==summary_75_100cm$compct_SOM_75_100cm & summary_100cm$compct_SOM_100cm!=0),])
1.836628+0.6887355+0.5739462+0.5739462

summary_all$SOM_compct_eql <- ifelse(summary_30cm$compct_SOM_30cm==summary_30_50cm$compct_SOM_30_50cm & summary_30cm$compct_SOM_30cm==summary_50_75cm$compct_SOM_50_75cm & summary_30cm$compct_SOM_30cm==summary_75_100cm$compct_SOM_75_100cm & summary_100cm$compct_SOM_100cm!=0, TRUE, FALSE)
table(summary_all$SOM_compct_eql)
sum(summary_all$hectares[summary_all$SOM_compct_eql==TRUE])
#15,988,385 hectares with same component percentage for % SOM in all layers

summary_all$SOCkg_compct_eql <- ifelse(summary_30cm$compct_SOCkg_30cm==summary_30_50cm$compct_SOCkg_30_50cm & summary_30cm$compct_SOCkg_30cm==summary_50_75cm$compct_SOCkg_50_75cm & summary_30cm$compct_SOCkg_30cm==summary_75_100cm$compct_SOCkg_75_100cm & summary_100cm$compct_SOCkg_100cm!=0, TRUE, FALSE)
table(summary_all$SOCkg_compct_eql)
sum(summary_all$hectares[summary_all$SOCkg_compct_eql==TRUE])
#15,860,293 hectares with same component percentage for SOC kg in all layers

sum(summary_all$hectares[!is.na(summary_all$SOM_30cm)])
#31,166,033 hectares with 0-30 SOM estimate
sum(summary_all$hectares[!is.na(summary_all$kgSOC.m2_30cm)])
#31,079,945 hectares with 0-30 kg SOC estimate

sum(summary_all$hectares[!is.na(summary_all$SOM_100cm)])
#29,630,957 hectares with 0-100 SOM estimate
sum(summary_all$hectares[!is.na(summary_all$kgSOC.m2_100cm)])
#29,507,050 hectares with 0-100 kg SOC estimate

sum(summary_all$hectares) #41,477,744

write.csv(summary_all, file.path(SSURGOdir, 'summaries', paste0('CA_mapunit_SOM_0_100cm_all_layers_', Sys.Date(), '.csv')), row.names=FALSE)
