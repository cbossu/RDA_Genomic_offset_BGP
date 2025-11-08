library(raster)
library(sp)
library(sf)
library(terra)
library(xlsx)
library(sf)
library(tidyverse)

#Load bioclim rasters for wc2.1

for (i in 1:19) {
    fileTemp <- paste("rasters/bioclim/wc2.1_30s_bio_", i ,".tif", sep="")
    name <- paste("bio", i, sep="")
    assign(name, raster(fileTemp))
    rm(name, fileTemp)
  }
rm(i)
##Load srtm raster
srtm <- raster("rasters/srtm/srtm.tif")

##loadNDVI
ndvi <- raster("rasters/NDVI/ndvimax.tif")
ndvistd <- raster("rasters/NDVI/ndvistd.tif")

##Load treecover, wind
qscat <- raster("rasters/other/qscat.tif")
tree <- raster("rasters/other/tree.tif")

#hii <-raster("rasters/other/hii.tiff")
##Combine rasters into one stack
compareRaster(srtm, ndvi, ndvistd, qscat, tree)
srtmstack=stack(srtm, ndvi, ndvistd, qscat, tree)
names(srtmstack) <- c("srtm", "ndvi", "ndvistd", "qscat", "tree")
bioclim= stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10,
               bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)
names(bioclim) <- c("bio01", "bio02", "bio03", "bio04", "bio05",
                    "bio06", "bio07", "bio08", "bio09", "bio10",
                    "bio11", "bio12", "bio13", "bio14", "bio15", 
                    "bio16", "bio17", "bio18", "bio19")

##Load bioclim rasters for both pathways
##Load ssp126 raster
ssp126 <- raster::stack("future_bioclim_2041-2060_ssp126/wc2.1_2.5m_bioc_CanESM5_ssp126_2041-2060.tif")
names(ssp126) <- c("bio01", "bio02", "bio03", "bio04", "bio05",
                   "bio06", "bio07", "bio08", "bio09", "bio10",
                   "bio11", "bio12", "bio13", "bio14", "bio15", 
                   "bio16", "bio17", "bio18", "bio19")
##Load ssp585 raster
ssp585 <- raster::stack("future_bioclim_2041-2060_ssp585/wc2.1_2.5m_bioc_CanESM5_ssp585_2041-2060.tif")
names(ssp585) <- c("bio01", "bio02", "bio03", "bio04", "bio05",
                   "bio06", "bio07", "bio08", "bio09", "bio10",
                   "bio11", "bio12", "bio13", "bio14", "bio15", 
                   "bio16", "bio17", "bio18", "bio19")

##remove individual rasters
rm(list= ls()[!(ls() %in% c('srtmstack','bioclim', 'ssp126', 'ssp585'))])

cawaSampled <- read.xlsx("cawa_seq_metadata.xlsx", sheetIndex=1) %>% rename(long=X.long,lat=X.lat)
cawaSampled.xy <- as.matrix(cawaSampled[,c("long","lat")])

##Extract current variables for sampled locations
cawaBioClim <- raster::extract(bioclim,  cawaSampled.xy)
cawaSrtm <- raster::extract(srtmstack,  cawaSampled.xy)
##Extract future variables for sampled locations
cawaSSP126 <- raster::extract(ssp126,  cawaSampled.xy)
cawaSSP585 <- raster::extract(ssp585,  cawaSampled.xy)
##Link data
cawaEnv<- data.frame(cawaSampled[, c("long","lat")], cawaBioClim, cawaSrtm)
#write.table(cawaEnv,file="cawaEnv_wc2.1.txt",row.names=FALSE,quote=F,sep='\t')
cawa126<- data.frame(cawaSampled[, c("long","lat")], cawaSSP126, cawaSrtm)
#write.table(cawa126, file="cawaFuture126_CH.txt",row.names=FALSE,quote=F,sep='\t')
cawa585<- data.frame(cawaSampled[, c("long","lat")], cawaSSP585, cawaSrtm)
#write.table(cawa585, file="cawaFuture585_CH.txt",row.names=FALSE,quote=F,sep='\t')

##Get breeding range for CAWA
cawa <- st_read("CAWA_eBird/CAWA.range_smooth.sf.WGS84.Ebird.shp") %>% dplyr::filter(season=="breeding")

##Limit variables to JUST the breeding range of Canada Warbler
bioclim_crop=crop(bioclim,cawa)
bioclim_cawa=mask(bioclim_crop,cawa)
#plot(bioclim_cawa)
srtmstack_crop=crop(srtmstack,cawa)
srtmstack_cawa=mask(srtmstack_crop,cawa)
#plot(srtmstack_cawa)
envVars <- stack(bioclim_cawa, srtmstack_cawa)

##Future vars
##Uses current vars for srtm, tree, qscat, and ndvi b/c no great predictive rasters for those

##SSP126
ssp126_crop=crop(ssp126,cawa)
ssp126_cawa=mask(ssp126_crop,cawa)
#plot(ssp126_cawa)
##Need slightly different extent for future
#bigger_raster <- projectRaster(raster5k, smaller_raster)- envVars is larger
envVars_proj <- projectRaster(envVars, ssp126_cawa)
srtmstack_cawa_proj <- projectRaster(srtmstack_cawa, ssp126_cawa)

#srtm_2.5 <- raster(extent(ssp126_cawa))
#srtm_2.5 <-  resample(srtmstack_cawa, ssp126_cawa, method = 'bilinear')
##check
#compareRaster(srtm_2.5, ssp126_cawa)
compareRaster(srtmstack_cawa_proj, ssp126_cawa)
future126 <- stack(ssp126_cawa, srtmstack_cawa_proj)

##SSP585
ssp585_crop=crop(ssp585,cawa)
ssp585_cawa=mask(ssp585_crop,cawa)
#plot(ssp585_cawa)
future585 <- stack(ssp585_cawa, srtmstack_cawa_proj)

