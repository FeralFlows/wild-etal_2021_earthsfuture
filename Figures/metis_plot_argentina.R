##################################################################
# Plot Demeter Output (land allocation) to spatial mapping       #
# Author: Mengqi Zhao                                            #
# Email: mengqiz@umd.edu                                         #
# Last Update: 2020-10-21                                        #  
##################################################################



rm(list = ls())

if("devtools" %in% rownames(installed.packages()) == F){install.packages("devtools")}
library(devtools)
if("metis" %in% rownames(installed.packages()) == F){devtools::install_github(repo="JGCRI/metis")}
library(metis)
if('gcamdata' %in% rownames(installed.packages()) == F){devtools::install_github(repo='JGCRI/gcamdata')}
library(gcamdata)
if('rgdal' %in% rownames(installed.packages()) == F){install.packages('rgdal')}
library(rgdal)
if('raster' %in% rownames(installed.packages()) == F){install.packages('raster')}
library(raster)
if('dplyr' %in% rownames(installed.packages()) == F){install.packages('dplyr')}
library(dplyr)
if('tidyr' %in% rownames(installed.packages()) == F){install.packages('tidyr')}
library(tidyr)
if('data.table' %in% rownames(installed.packages()) == F){install.packages('data.table')}
library(data.table)
if('snow' %in% rownames(installed.packages()) == F){install.packages('snow')}
library(snow)
if('sp' %in% rownames(installed.packages()) == F){install.packages('sp')}
library(sp)
if('stringr' %in% rownames(installed.packages()) == F){install.packages('stringr')}
library(stringr)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work.dir <- getwd()
demeter_output <- 'E:/NEXO-UA/Demeter/example/outputs'
select_folders <- c('reference_2020-10-28_20h29m05s',
                    'impacts_2020-10-29_08h31m08s',
                    'policy_2020-10-29_09h21m42s')
scenarios <- c('reference', 'impacts', 'policy')
select_output <- paste(demeter_output, select_folders, sep = '/')
output_files <- list.files(select_output, pattern = '0p5deg', recursive = TRUE, full.names = TRUE)




# Create Basins and hydrologic watersheds within Argentina (better solution)
m_basin <- metis::mapGCAMBasins
m_hydroshed <- metis::mapHydroShed3
m_country <- metis::mapCountries

m_argentina <- m_country[m_country@data$subRegion %in% c("Argentina"),]; sp::plot(m_argentina)
m_argentina_ext <- m_country[m_country@data$subRegion %in% c("Argentina", 'Chile', 'Bolivia', 'Paraguay', 'Brazil', 'Uruguay'),]; sp::plot(m_argentina_ext)

# Basins within Argentina
m_argentina_basin <- sp::spTransform(m_argentina, raster::crs(m_basin))
m_argentina_basin <- raster::crop(m_basin, m_argentina_basin)
m_argentina_basin@data <-  droplevels(m_argentina_basin@data)
sp::plot(m_argentina_basin)

# Basins within extended Argentina Region
m_argentina_ext_basin <- sp::spTransform(m_argentina_ext, raster::crs(m_basin))
m_argentina_ext_basin <- raster::crop(m_basin, m_argentina_ext_basin)
m_argentina_ext_basin@data <- droplevels(m_argentina_ext_basin@data)
sp::plot(m_argentina_ext_basin)

# Countrys within extended Argentina Region
m_argentina_ext_country <- sp::spTransform(m_argentina_ext, raster::crs(m_country))
m_argentina_ext_country <- raster::crop(m_country, m_argentina_ext_country)
m_argentina_ext_country@data <- droplevels(m_argentina_ext_country@data)
sp::plot(m_argentina_ext_country)

# Hydrosheds within Argentina
m_argentina_hydroshed <- sp::spTransform(m_argentina, raster::crs(m_hydroshed))
m_argentina_hydroshed <- raster::crop(m_hydroshed, m_argentina_hydroshed)
m_argentina_hydroshed@data <- droplevels(m_argentina_hydroshed@data)
sp::plot(m_argentina_hydroshed)

# Hydrosheds within extended Argentina
m_argentina_ext_hydroshed <- sp::spTransform(m_argentina, raster::crs(m_hydroshed))
m_argentina_ext_hydroshed <- raster::crop(m_hydroshed, m_argentina_ext_hydroshed)
m_argentina_ext_hydroshed <- raster::bind(m_argentina_ext_country, m_argentina_ext_hydroshed)
m_argentina_ext_hydroshed@data <- droplevels(m_argentina_ext_hydroshed@data)
sp::plot(m_argentina_ext_hydroshed)

# Plot Argentina Maps at basin and hydroshed level
metis.map(m_argentina_basin, fileName = 'Argentina_basin', labels = T, labelsSize = 0.6)
metis.map(m_argentina_ext_basin, fileName = 'Argentina_ext_basin', labels = T, labelsSize = 0.6)
metis.map(m_argentina_ext_country, fileName = 'Argentina_ext_country', labels = T, labelsSize = 0.6)
metis.map(m_argentina_hydroshed, fileName = 'Argentina_hydroshed', labels = T, labelsSize = 0.2)
metis.map(m_argentina_ext_hydroshed, fileName = 'Argentina_ext_hydroshed', labels = F, labelsSize = 0.2)

# Shapefiles for selected region
shape <- m_argentina_hydroshed
shape_ext <- m_argentina_ext_hydroshed # extended region for plotting
boundary <- raster::extent(bbox(m_argentina_hydroshed))

# Create bounday shapefile in SpatialPolygonsDataFrame format
x1 <- boundary@xmin - 1.0
x2 <- boundary@xmax + 1.0
y1 <- boundary@ymin - 1.0
y2 <- boundary@ymax + 1.0
shape_boundary <- sp::Polygon(cbind(c(x1,x1,x2,x2,x1),c(y1,y2,y2,y1,y1)))
shape_boundary <- sp::Polygons(list(shape_boundary), ID = "A")
shape_boundary <- sp::SpatialPolygons(list(shape_boundary))
sp::proj4string(shape_boundary) <- m_argentina_ext_hydroshed@proj4string
df_boundary <- matrix(data = NA)
rownames(df_boundary) <- 'A'
shape_boundary <- sp::SpatialPolygonsDataFrame(shape_boundary, data = as.data.frame(df_boundary))

region <- 'Argentina'
subRegType_i <- 'hydroshed' # hydroshed
gcm <- 'MIROC-ESM-CHEM'
rcp <- 'rcp6p0'
main_folder <- paste(region, gcm, rcp, sep = '_')


#------------------------------- For Demeter -------------------------------#
#Read Demeter data in paralell####
cl <- makeSOCKcluster(4)

#function to read before run####
read_files <- function(file, header = TRUE, sep = 'auto', ...){
  data <- data.table::fread(file, header = header, sep = sep, ...)
  filename <- basename(file[1])
  scenario <- strsplit(basename(filename), '\\_|\\.')[[1]][2]
  year <- strsplit(basename(filename), '\\_|\\.')[[1]][4]
  data$year <- year
  data$scenario <- scenario
  # data$file_type <- filename
  return(data)
}

# basename is a function to choose file by name
all_data <- parLapply(cl, output_files, read_files)
stopCluster(cl)

col_gather <- c('water', 'forest', 'shrub', 'grass',
                'urban', 'snow', 'sparse', 'corn_irr',
                'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
                'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
                'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
                'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
                'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
                'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
                'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')

df_data <- rbindlist(all_data, idcol = FALSE) %>% 
  dplyr::rename(gridcode = pkey_0p5_deg, lon = longitude, lat = latitude) %>%
  gather(key = 'crop', value = 'value', col_gather)
rm(all_data)
# Filter data for Argentina region
df_data_argentina <- df_data %>% 
  filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)

#------------------------------- DEMETER MAP PLOTTING -------------------------------#
# Plot all the crops for 3 scenarios: reference, impacts, and policy
# The output figures include crop land allocation fraction under individual scenario,
# and the absolute difference and percent difference between (1) impacts and reference;
# and (2) policy and reference
# Note: this plotting may take hours. You may select less crop types in order to reduce run time
if(T){
  crop_i <-c('water', 'forest', 'shrub', 'grass',
             'urban', 'snow', 'sparse', 'corn_irr',
             'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
             'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
             'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
             'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
             'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
             'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
             'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  # crop_i <-c('corn_irr',
  #            'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
  #            'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
  #            'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
  #            'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
  #            'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
  #            'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd',
  #            'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  # crop_i <- c('biomass_grass_rfd', 'biomass_tree_rfd')
  
  # crop_i <- c('biomass_grass_irr', 'biomass_tree_irr')
  
  nameAppend_i <- c('_LandUse')
  df <- df_data_argentina
  scenario_i <- c('reference', 'impacts', 'policy')
  year_i <- c(2020, 2030, 2040, 2050)

  data_2 <- df %>% 
    filter(scenario %in% scenario_i & year %in% year_i & crop %in% crop_i) %>% 
    tidyr::unite(class, c(crop, scenario, year), sep = '-', remove = TRUE) %>% 
    dplyr::select(-gridcode); str(data_2)
  
  poly_table <- metis.grid2poly(gridFiles = data_2,
                                subRegShape = shape,
                                aggType = 'depth',
                                subRegCol = 'subRegion',
                                subRegType = subRegType_i,
                                nameAppend = '_Demeter',
                                folderName = main_folder,
                                saveFiles = T)
  
  # Find and add missing polygons
  poly_table_temp <- poly_table[0,]
  for (crop in crop_i){
    for (scenario in scenario_i){
      for (year in year_i){
        class_temp <- paste(crop, scenario, year, sep = '-')
        temp <- poly_table %>% 
          filter(class == class_temp) %>% 
          as.data.frame()
        missing <- dplyr::setdiff(shape$subRegion, temp$subRegion) %>% 
          as.data.frame() %>% 
          dplyr::rename(subRegion = '.')
        n_miss <- length(missing$subRegion)
        print(paste0('Missing number of subregions for ', crop, ' | ', scenario, ' | ', year, ' are: ', n_miss))
        if(n_miss > 0){
          df_missing <- data.frame(subRegion = missing$subRegion,
                                   value = rep(1e-10, times = n_miss),
                                   class = rep(class_temp, times = n_miss),
                                   scenario = rep('scenario', times = n_miss),
                                   x = rep('x', times = n_miss),
                                   param = rep('param', times = n_miss),
                                   aggType = rep('depth', times = n_miss),
                                   subRegType = rep(subRegType_i, times = n_miss))
          temp <- rbind(temp, df_missing)
        }
        poly_table_temp <- rbind(poly_table_temp, temp)
      }
    }
  }
  
  # regular plot
  poly_table_2 <- poly_table_temp %>%
    tidyr::separate(class, sep = '-', into = c('class', 'scenario', 'x')) %>% 
    dplyr::mutate(param = class,
                  classPalette = 'pal_green',
                  value = if_else(value < 1e-6, 1e-6, value))
  
  metis.mapsProcess(polygonTable = poly_table_2,
                    subRegShape = shape,
                    subRegCol = 'subRegion',
                    subRegType = subRegType_i,
                    scenRef = 'reference',
                    nameAppend = nameAppend_i,
                    folderName = main_folder,
                    xRange = year_i,
                    scaleRange = c(0,1),
                    mapTitleOn = F,
                    legendFixedBreaks = 8,
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    extdendedLabelSize = 0.7,
                    cropToBoundary = F,
                    extension = T,
                    expandPercent = 25,
                    # classPalette = 'pal_green',
                    classPaletteDiff = 'pal_div_BrGn',
                    facetCols = 4,
                    animateOn = F,
                    pdfpng = 'pdf')
  
  if(T){
    # aggregate biomass_grass_irr and biomass_tree_irr together
    poly_table_biomass_irr <- poly_table_2 %>% 
      filter(class %in% c('biomass_grass_irr', 'biomass_tree_irr')) %>% 
      mutate(class = 'biomass_irr',
             param = 'biomass_irr',) %>% 
      group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      summarise(value = sum(value)) %>% 
      ungroup()
    
    poly_table_biomass_rfd <- poly_table_2 %>% 
      filter(class %in% c('biomass_grass_rfd', 'biomass_tree_rfd')) %>% 
      mutate(class = 'biomass_rfd',
             param = 'biomass_rfd') %>% 
      group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      summarise(value = sum(value)) %>% 
      ungroup()
    
    
    # plot by crop_irr, crop_rfd
    poly_table_crop <- poly_table_2
    poly_table_crop$param[grepl('irr', poly_table_crop$param)] <- 'crop_irr'
    poly_table_crop$param[grepl('rfd', poly_table_crop$param)] <- 'crop_rfd'
    poly_table_crop$class <- poly_table_crop$param
    poly_table_crop <- poly_table_crop %>% 
      dplyr::filter(class %in% c('crop_irr', 'crop_rfd')) %>% 
      dplyr::group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      dplyr::summarise(value = sum(value)) %>% 
      dplyr::ungroup()
    
    # plot by crop total
    poly_table_crop_all <- poly_table_crop
    poly_table_crop_all[grepl('crop_irr|crop_rfd', poly_table_crop_all)] <- 'crop_all'
    poly_table_crop_all <- poly_table_crop_all %>% 
      dplyr::group_by(aggType, subRegType, scenario, x, subRegion, param, class, classPalette) %>% 
      dplyr::summarise(value = sum(value)) %>% 
      dplyr::ungroup()
    
    # plot by land allocation type
    # c('water', 'forest', 'shrub', 'grass', 'urban', 'snow', 'sparse')
    poly_table_landAlloc <- poly_table_2 %>% 
      dplyr::filter(class %in% c('water', 'forest', 'shrub', 'grass', 'urban', 'snow', 'sparse')) %>% 
      dplyr::bind_rows(poly_table_crop_all %>% dplyr::mutate(class = 'crop')) %>% 
      dplyr::mutate(param = 'landAlloc_All')
  }
  
  poly_table_list <- list(poly_table_biomass_irr, poly_table_biomass_rfd, poly_table_crop, poly_table_crop_all, poly_table_landAlloc)
  
  # Plot all the crop land allocations by crop types
  for(poly_table_i in poly_table_list){
    metis.mapsProcess(polygonTable = poly_table_i,
                      subRegShape = shape,
                      subRegCol = 'subRegion',
                      subRegType = subRegType_i, # 'hydroshed',
                      scenRef = 'reference',
                      nameAppend = nameAppend_i,
                      folderName = main_folder,
                      xRange = year_i,
                      mapTitleOn = F,
                      legendFixedBreaks = 8,
                      boundaryRegShape = shape_boundary,
                      extendedLabels = T,
                      cropToBoundary = F,
                      extension = T,
                      expandPercent = 25,
                      # classPalette = 'pal_green',
                      classPaletteDiff = 'pal_div_BrGn',
                      facetCols = 4,
                      animateOn = F)
  }
  
  
}


# Spatial mapping function for Tethys and Xanthos output in km3
sp_mapping <- function(gridFiles, subRegShape, subRegShapeExt, boundaryRegShape,
                       aggType, subRegCol, subRegType, nameAppend, folderName, years, ...){
  # Aggregate Tethys 0.5 degree cells to selected polygons (e.g. hydroshed)
  poly_table <- metis.grid2poly(gridFiles = gridFiles,
                                subRegShape = subRegShape,
                                aggType = aggType, # Tethys output in km3, so use 'vol' to aggregate
                                subRegCol = subRegCol,
                                subRegType = subRegType,
                                nameAppend = nameAppend,
                                folderName = folderName,
                                saveFiles = T)
  
  # Find missing polygons
  shape_area <- data.frame(subRegion = subRegShape@data$subRegion, area = raster::area(subRegShape)) # area is in m2
  
  poly_table_temp <- poly_table[0,]
  for(scenario_t in unique(poly_table$scenario)){
    for(year_t in unique(years)){
      for(param_t in unique(poly_table$param)){
        for(class_t in unique(poly_table$class)){
          temp <- poly_table %>% 
            filter(class == class_t & param == param_t & scenario == scenario_t & x == year_t) %>% 
            as.data.frame()
          missing <- dplyr::setdiff(subRegShape$subRegion, temp$subRegion) %>% 
            as.data.frame() %>% 
            dplyr::rename(subRegion = '.')
          n_miss <- length(missing$subRegion)
          if(n_miss > 0){
            sprintf('-------------------- Adding %s missing polygons --------------------', n_miss)
            missing <- data.frame(subRegion = missing$subRegion,
                                  value = rep(1e-10, times = n_miss),
                                  year = paste('X', rep(year_t, times = n_miss), sep = ''),
                                  scenario = rep(scenario_t, times = n_miss),
                                  x = rep(year_t, times = n_miss),
                                  class = rep(class_t, times = n_miss),
                                  param = rep(param_t, times = n_miss),
                                  aggType = rep(aggType, times = n_miss),
                                  subRegType = rep(subRegType, times = n_miss)) %>% 
              dplyr::select(names(temp))
            
            # Add missing polygons and convert unit from volume to depth
            temp <- rbind(temp, missing)
          }
          poly_table_temp <- rbind(poly_table_temp, temp)
        }
      }
    }
  }
  
  
  if(aggType == 'vol'){
    poly_table_2 <- poly_table_temp %>% 
      left_join(shape_area, by = subRegCol) %>% 
      mutate(value = (value*1000^3/area)*1000) # convert from vol (km3) to depth (mm)
  }else{
    poly_table_2 <- poly_table_temp
    # facetCols <- ceiling(length(unique(poly_table_2$class)))
  }
  
  
  facetCols <- if (length(years) <= 4) {
    length(years)
  }
  else if (length(years) <= 10) {
    ceiling(length(years)/2)
  }
  else {
    ceiling(sqrt(length(years)))
  }
  
  # facetCols <- 3

  # Plot maps
  metis.mapsProcess(polygonTable = poly_table_2,
                    subRegShape = subRegShapeExt,
                    subRegCol = subRegCol,
                    subRegType = subRegType,
                    nameAppend = nameAppend,
                    folderName = folderName,
                    mapTitleOn = F,
                    legendFixedBreaks = 7,
                    # scaleRange = c(0,1),
                    # xRange = years,
                    boundaryRegShape = boundaryRegShape,
                    extendedLabels = T,
                    extension = T,
                    # expandPercent = 25,
                    cropToBoundary = F,
                    classPalette = 'pal_wet', # pal_wet for water, pal_green for cropland
                    animateOn = F,
                    facetCols = facetCols)
  
  return(poly_table_2)
}

# Calculate grid cell area and convert from water volumn to depth
grid_vol_to_depth <- function(grid_data, shape){
  grid <- sp::SpatialPointsDataFrame(sp::SpatialPoints(
    coords = (cbind(grid_data$lon, grid_data$lat))),
    data = grid_data %>% dplyr::select(lat, lon))
  
  sp::gridded(grid) <- TRUE
  grid_raster <- raster::stack(grid); grid_raster
  grid_poly <- raster::rasterToPolygons(grid_raster)
  sp::proj4string(grid_poly) <- sp::proj4string(shape)
  
  # Calculate grid cell area and left join based on lat and lon
  poly_area <- grid_poly@data %>%
    mutate(area=raster::area(grid_poly)) # area in m2
  grid_data <- grid_data %>%
    left_join(poly_area, by=c('lat', 'lon')) %>% 
    mutate(value = value*1000^3/area*1000) %>% 
    dplyr::select(-area)# depth in mm
  
  return(grid_data)
}

#------------------------------- Tethys Output -------------------------------#

years <- sprintf('%s', seq(from = 2005, to = 2050, by = 5))
# Total Water Demand (irrigation, domestic, electricity, livestock, manufacturing, mining, non-ag)
tethys_mapping <- function(wd_sector, years, boundary, scenario, gcm, rcp){
  tethys_output <- 'E:/NEXO-UA/Tethys/example/Output'
  gcam_v <- 'gcam_5p1'
  time_scale <- '2005-2050'
  folder_name <- paste(gcam_v, gcm, rcp, time_scale, scenario, sep = '_')
  file_name <- paste0(wd_sector, '_km3peryr', '.csv')
  file_path <- paste(tethys_output, folder_name, file_name, sep = '/')
  wd_data <- data.table::fread(file_path, header = TRUE)
  
  wd <- wd_data %>% 
    as.data.frame() %>% 
    gather(key = 'year', value = 'value', years) %>% 
    dplyr::select(-ID, -ilon, -ilat)
  
  # Filter data to Argentina area
  wd_argentina <- wd %>% 
    filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)
  
  # Convert volume to depth based on each grid cell
  wd_argentina_depth <- grid_vol_to_depth(wd_argentina, shape)
  
  nameAppend <- paste('_Tethys', scenario, sep = '_')
  
  demand <- sp_mapping(wd_argentina_depth, shape, shape, shape_boundary, 'depth',
                      'subRegion', 'hydroshed', nameAppend, main_folder, years)
  return(demand)
}
demand_reference <- tethys_mapping('wdtotal', years, boundary, 'Reference', gcm, rcp)
demand_impacts <- tethys_mapping('wdtotal', years, boundary, 'Impacts', gcm, rcp)
demand_policy <- tethys_mapping('wdtotal', years, boundary, 'Policy', gcm, rcp)

# Plot each water withdrawal sector in one frame
tethys_indv_mapping <- function(gcm, rcp, scenario, boundary, shape){
  tethys_output <- 'E:/NEXO-UA/Tethys/example/Output'
  gcam_v <- 'gcam_5p1'
  time_scale <- '2005-2050'
  folder_name <- paste(gcam_v, gcm, rcp, time_scale, scenario, sep = '_')
  folder_path <- paste(tethys_output, folder_name, sep = '/')
  file_list <- list.files(folder_path, pattern = 'km3peryr', recursive = TRUE, full.names = TRUE)
  file_list <- file_list[-c(7,8)]
  
  # Read Tethys files
  #Read Demeter data in paralell####
  cl <- makeSOCKcluster(4)
  clusterExport(cl, list('%>%', 'gather'), envir = .GlobalEnv)
  
  #function to read before run####
  read_files_tethys <- function(file, header = TRUE, sep = 'auto', ...){
    data <- data.table::fread(file, header = header, sep = sep, ...)
    filename <- basename(file[1])
    class <- tolower(strsplit(basename(filename), '\\_|\\.')[[1]][1])
    patterns <- c('wddom', 'wdelec', 'wdirr', 'wdliv', 'wdmfg', 'wdmin')
    replacement <- c('Domestic', 'Electric', 'Irrigation', 'Livestock', 'Manufacturing', 'Mining')
    for (i in seq_along(patterns)){
      class <- gsub(patterns[i], replacement[i], class, perl = TRUE)
    }
    data$class <- class
    
    years <- sprintf('%s', seq(from = 2005, to = 2050, by = 5))
    data <- data %>% 
      as.data.frame() %>% 
      tidyr::gather(key = 'year', value = 'value', years) %>% 
      dplyr::select(-ID, -ilon, -ilat)
    
    return(data)
  }
  
  # basename is a function to choose file by name
  wd_indv <- parLapply(cl, file_list, read_files_tethys)
  stopCluster(cl)
  
  wd_indv_df <- rbindlist(wd_indv, idcol = FALSE)
  
  # Filter data to boundary area
  wd_indv_boundary <- wd_indv_df %>% 
    dplyr::filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)
  
  # Convert volume to depth based on each grid cell
  wd_indv_boundary_depth <- grid_vol_to_depth(wd_indv_boundary, shape)
  
  nameAppend <- paste('_Tethys', scenario, 'Indv', sep = '_')
  
  demand <- sp_mapping(wd_indv_boundary_depth, shape, shape, shape_boundary, 'depth',
                       'subRegion', 'hydroshed', nameAppend, main_folder, years)
}

tethys_indv_mapping(gcm, rcp, 'Reference', boundary, shape)
tethys_indv_mapping(gcm, rcp, 'Impacts', boundary, shape)
tethys_indv_mapping(gcm, rcp, 'Policy', boundary, shape)



#------------------------------- Xanthos Output -------------------------------#
# for original Xanthos output
run <- 'clim_impacts'
runoff_var_name <- 'q_km3peryear'
time_scale <- '1950_2099'

xanthos_folder <- paste(run, gcm, rcp, sep = '_')
xanthos_output <- paste('E:/NEXO-UA/Xanthos/example/output', xanthos_folder, sep = '/')
runoff_file <- paste0(paste(runoff_var_name, gcm, rcp, time_scale, sep = '_'), '.csv')


# Read Xanthos runoff output in km3
runoff <- data.table::fread(file = paste(xanthos_output, runoff_file, sep = '/'), header = TRUE)
coord_0p5deg <- data.table::fread(file = 'coordinates.csv', header = FALSE) %>% 
  dplyr::rename(gridcode = V1,
         longitude = V2,
         latitude = V3,
         x = V4,
         y = V5)

runoff_years <- sprintf('%s', seq(from = 1950, to = 2099, by = 1))
years <- sprintf('%s', seq(from = 2005, to = 2050, by = 5))

runoff_argentina <- runoff %>% 
  mutate(lat = coord_0p5deg$latitude,
         lon = coord_0p5deg$longitude) %>% 
  gather(key = 'year', value = 'value', all_of(runoff_years)) %>% 
  filter(year %in% years, lat <= boundary@ymax &
           lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin) %>% 
  dplyr::select(-id)

runoff_argentina_depth <- grid_vol_to_depth(runoff_argentina, shape)

supply <- sp_mapping(runoff_argentina_depth, shape, shape, shape_boundary, 'depth',
                     'subRegion', 'hydroshed', '_Xanthos', main_folder, years)

# Plot Historical runoff from projections
# Run this step with years spanning from 1950 - 2010
years <- sprintf('%s', seq(from = 1950, to = 2010, by = 1))
supply <- sp_mapping(runoff_argentina_depth, shape, shape, shape_boundary, 'depth',
                     'subRegion', 'hydroshed', '_XanthosHist', main_folder, years)


if(F){
  # for historical water supply
  gcm <- 'watch+wfdei'
  rcp <- '1970_2010'
  time_scale <- '1970_2010'
  xanthos_folder <- paste(run, gcm, time_scale, sep = '_')
  xanthos_output <- paste('E:/NEXO-UA/Xanthos/example/output', xanthos_folder, sep = '/')
  runoff_file <- paste0(paste(runoff_var_name, gcm, time_scale, sep = '_'), '.csv')

  # Read Xanthos runoff output in km3
  runoff <- data.table::fread(file = paste(xanthos_output, runoff_file, sep = '/'), header = TRUE)
  coord_0p5deg <- data.table::fread(file = 'coordinates.csv', header = FALSE) %>% 
    dplyr::rename(gridcode = V1,
           longitude = V2,
           latitude = V3,
           x = V4,
           y = V5)
  
  years_hist <- c('2005', '2010')
  runoff_years <- sprintf('%s', seq(from = 1970, to = 2010, by = 1))
  runoff_argentina_hist <- runoff %>% 
    mutate(lat = coord_0p5deg$latitude,
           lon = coord_0p5deg$longitude) %>% 
    gather(key = 'year', value = 'value', all_of(runoff_years)) %>% 
    filter(year %in% runoff_years, lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin) %>% 
    dplyr::select(-id)
  
  runoff_argentina_hist_depth <- grid_vol_to_depth(runoff_argentina_hist, shape)
  
  supply_hist <- sp_mapping(runoff_argentina_hist_depth, shape, shape, shape_boundary, 'depth',
                       'subRegion', 'hydroshed', '_XanthosHist', 'Argentina_Hist', runoff_years)
}


#------------------------------- Water Scarcity -------------------------------#

scarcity_mapping <- function(demand, supply, scenario){
  poly_table_scarcity <- demand %>% 
    left_join(supply, by = c('subRegion',
                             'year',
                             'x',
                             'class',
                             'scenario',
                             'param',
                             'aggType',
                             'subRegType')) %>% 
    mutate(value = value.x/value.y) %>% 
    dplyr::select(-value.x, -value.y) 
  
  numeric2Cat_param <- list("param")
  numeric2Cat_breaks <- list(c(0,0.1,0.2,0.4,Inf))
  numeric2Cat_labels <- list(c("None","Low","Moderate","Severe"))
  numeric2Cat_palette <- list(c("None"="#3288BD","Low"="#ABDDA4","Moderate"="#FDAE61","Severe"="#9E0142"))
  numeric2Cat_legendTextSize <- list(c(1))
  numeric2Cat_list <-list(numeric2Cat_param=numeric2Cat_param,
                          numeric2Cat_breaks=numeric2Cat_breaks,
                          numeric2Cat_labels=numeric2Cat_labels,
                          numeric2Cat_palette=numeric2Cat_palette,
                          numeric2Cat_legendTextSize=numeric2Cat_legendTextSize)
  nameAppend <- paste('_Scarcity', scenario, sep = '_')
  
  # Plot maps
  metis.mapsProcess(polygonTable = poly_table_scarcity,
                    subRegShape = shape,
                    subRegCol = 'subRegion',
                    subRegType = 'hydroshed',
                    nameAppend = nameAppend,
                    folderName = main_folder,
                    mapTitleOn = F,
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    extension = T,
                    cropToBoundary = F,
                    facetCols = 5,
                    animateOn = F,
                    numeric2Cat_list = numeric2Cat_list)
}

scarcity_mapping(demand_reference, supply, 'Reference')
scarcity_mapping(demand_impacts, supply, 'Impacts')
scarcity_mapping(demand_policy, supply, 'Policy')



#------------------------------- Demeter Historical Land Allocation -------------------------------#
modis <- data.table::fread(file = 'E:/NEXO-UA/Demeter/example/inputs/observed/gcam_reg32_basin235_modis_mirca_2000_0p5deg_MZ.csv',
                           header = TRUE)
df_modis <- modis %>% 
  rename(lat = latitude,
         lon = longitude) %>% 
  dplyr::select(-pkey_0p5_deg) %>% 
  gather(key = 'class', value = 'value', append(col_gather, 'crops')) %>% 
  filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)

df_modis_select <- df_modis %>% 
  filter(class %in% c('grass', 'forest', 'shrub', 'urban', 'crops'))


poly_table_modis <- sp_mapping(df_modis_select, shape, shape_ext, shape_boundary, 'vol',
           'subRegion', 'hydroshed', '_DemeterModis', 'Argentina', '2000')
