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
# m_region <- metis::mapGCAMReg32

m_argentina <- m_country[m_country@data$subRegion %in% c("Argentina"),]; sp::plot(m_argentina)
m_argentina_ext <- m_country[m_country@data$subRegion %in% c("Argentina", 'Chile', 'Bolivia', 'Paraguay', 'Brazil', 'Uruguay'),]; sp::plot(m_argentina_ext)
# m_chile <- m_country[m_country@data$subRegion %in% c("Chile"),]; sp::plot(m_chile)
# m_argentina <- m_region[m_region@data$subRegion %in% c("Argentina"),]; sp::plot(m_argentina)

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
  rename(gridcode = pkey_0p5_deg, lon = longitude, lat = latitude) %>%
  gather(key = 'crop', value = 'value', col_gather)

# Filter data for Argentina region
df_data_argentina <- df_data %>% 
  filter(lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin)

#------------------------------- DEMETER MAP PLOTTING -------------------------------#
# Plot all the crops for 3 scenarios: reference, impacts, and policy
# The output figures include crop land allocation fraction under individual scenario,
# and the absolute difference and percent difference between (1) impacts and reference;
# and (2) policy and reference
# Note: this plotting may take hours. You may select less crop types in order to reduce run time
if(F){
  # crop_i <-c('water', 'forest', 'shrub', 'grass',
  #            'urban', 'snow', 'sparse', 'corn_irr',
  #            'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
  #            'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
  #            'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
  #            'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
  #            'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
  #            'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
  #            'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  # crop_i <-c('corn_irr',
  #            'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
  #            'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
  #            'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
  #            'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
  #            'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
  #            'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd',
  #            'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  crop_i <- c('biomass_grass_irr', 'biomass_tree_irr')
  # # 
  # crop_i <- c('fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
  #              'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
  #              'wheat_rfd',
  #              'biomass_grass_irr', 'biomass_grass_rfd')

  # crop_i <- c('shrub', 'sparse')
  # crop_i <- c('misccrop_irr', 'fodderherb_irr')
  region <- 'Argentina'
  gcm <- 'MIROC-ESM-CHEM'
  rcp <- 'rcp6p0'
  main_folder <- paste(region, gcm, rcp, sep = '_')
  nameAppend_i <- c('')
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
                                subRegType = 'hydroshed',
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
          rename(subRegion = '.')
        n_miss <- length(missing$subRegion)
        if(n_miss > 0){
          df_missing <- data.frame(subRegion = missing$subRegion,
                                   value = rep(1e-10, times = n_miss),
                                   class = rep(class_temp, times = n_miss),
                                   scenario = rep('scenario', times = n_miss),
                                   x = rep('x', times = n_miss),
                                   param = rep('param', times = n_miss),
                                   aggType = rep('depth', times = n_miss),
                                   subRegType = rep('hydroshed', times = n_miss))
          temp <- rbind(temp, df_missing)
        }
        poly_table_temp <- rbind(poly_table_temp, temp)
      }
    }
  }
  
  # regular plot
  poly_table_2 <- poly_table_temp %>%
    tidyr::separate(class, sep = '-', into = c('class', 'scenario', 'x')) %>% 
    mutate(param = class)
  
  if(F){
    # aggregate biomass_grass_irr and biomass_tree_irr together
    poly_table_biomass <- poly_table_2 %>% 
      mutate(class = 'biomass_irr',
             param = 'biomass_irr') %>% 
      group_by(aggType, subRegType, scenario, x, subRegion, param, class) %>% 
      summarise(value = sum(value)) %>% 
      ungroup()
    
    
    # plot by crop_irr, crop_rfd
    poly_table_crop <- poly_table_2
    poly_table_crop$param[grepl('irr', poly_table_crop$param)] <- 'crop_irr'
    poly_table_crop$param[grepl('rfd', poly_table_crop$param)] <- 'crop_rfd'
    poly_table_crop$class <- poly_table_crop$param
    poly_table_crop <- poly_table_crop %>% 
      group_by(aggType, subRegType, scenario, x, subRegion, param, class) %>% 
      summarise(value = sum(value)) %>% 
      ungroup()
    
    # plot by crop total
    poly_table_crop_all <- poly_table_crop
    poly_table_crop_all[grepl('crop_irr|crop_rfd', poly_table_crop_all)] <- 'crop_all'
    poly_table_crop_all <- poly_table_crop_all %>% 
      group_by(aggType, subRegType, scenario, x, subRegion, param, class) %>% 
      summarise(value = sum(value)) %>% 
      ungroup
  }
  
  
  metis.mapsProcess(polygonTable = poly_table_biomass,
                    subRegShape = shape_ext,
                    subRegCol = 'subRegion',
                    subRegType = 'hydroshed',
                    scenRef = 'reference',
                    nameAppend = nameAppend_i,
                    folderName = main_folder,
                    xRange = year_i,
                    # scaleRange = data.frame(param = unique(poly_table_2$param), min = c(0,0,0,0), max = c(1,1,1,1)),
                    mapTitleOn = F,
                    legendFixedBreaks = 6,
                    # scaleRange = c(0,1),
                    # scaleRangeDiffAbs = c(),
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    cropToBoundary = F,
                    extension = T,
                    expandPercent = 25,
                    classPalette = 'pal_green',
                    # classPaletteDiff = 'pal_RdBlu',
                    facetCols = 4,
                    animateOn = F)
  
}


if(F){
  # Plot crop land allocation fraction in group of four crops for years and scenarios
  # of your choice
  sp_alloc <- function(args_ls, df, scenario_i, year_i, shape){
    # browser()
    class_i <- args_ls[1:4]
    nameAppend_i <- args_ls[5]
    
    data_2 <- df %>% 
      filter(scenario %in% scenario_i & year %in% year_i & class %in% class_i) %>% 
      tidyr::unite(class, c(class, scenario, year), sep = '-', remove = TRUE) %>% 
      dplyr::select(-gridcode); str(data_2)
    
    poly_table <- metis.grid2poly(gridFiles = data_2,
                                  subRegShape = shape,
                                  aggType = 'depth',
                                  subRegCol = 'subRegion',
                                  subRegType = 'hydroshed',
                                  nameAppend = '_ArgentinaFrac',
                                  folderName = 'Argentina',
                                  saveFiles = T)
    
    # poly_table_2 <- poly_table %>%
    #   tidyr::separate(class, sep = '-', into = c('class', 'scenario', 'x')) %>% 
    #   mutate(param = class)
    
    metis.mapsProcess(polygonTable = poly_table,
                      subRegShape = shape,
                      subRegCol = 'subRegion',
                      subRegType = 'hydroshed',
                      nameAppend = nameAppend_i,
                      folderName = 'Argentina',
                      mapTitleOn = F,
                      legendFixedBreaks = 9,
                      cropToBoundary = F,
                      extension = T,
                      expandPercent = 25,
                      classPalette = 'pal_green',
                      facetCols = 8,
                      animateOn = F)
  }
  
  
  scenario_i <- 'reference'
  year_i <- c(2020, 2030, 2040, 2050)
  # class_i <- c('biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  
  args_ls <- list(C1 = c('water', 'forest', 'shrub', 'grass', '_Group1'),
                  C2 = c('urban', 'snow', 'sparse', 'corn_irr', '_Group2'),
                  C3 = c('fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr', '_Group3'),
                  C4 = c('oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr', '_Group4'),
                  C5 = c('root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd', '_Group5'),
                  C6 = c('fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd', '_Group6'),
                  C7 = c('oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd', '_Group7'),
                  C8 = c('root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland', '_Group8'),
                  C9 = c('biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd', '_Group9'))
  
  lapply(args_ls, sp_alloc, df_data_argentina, scenario_i, year_i, m_argentina_hydroshed)
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
  shape_area <- data.frame(subRegion = subRegShape@data$subRegion, area =  subRegShape@data$SUB_AREA) # area is in km2
  
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
            rename(subRegion = '.')
          n_miss <- length(missing$subRegion)
          if(n_miss > 0){
            sprintf('-------------------- Adding %s missing polygons --------------------', n_miss)
            missing <- data.frame(subRegion = missing$subRegion,
                                  value = rep(1e-10, times = n_miss),
                                  year = paste('X', rep(year, times = n_miss), sep = ''),
                                  scenario = rep(scenario, times = n_miss),
                                  x = rep(year, times = n_miss),
                                  class = rep(class, times = n_miss),
                                  param = rep(param, times = n_miss),
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
      mutate(value = (value/area)*1000*1000) # convert from km to mm
  }else{
    poly_table_2 <- poly_table_temp
    # facetCols <- ceiling(length(unique(poly_table_2$class)))
  }
  
  facetCols <- ifelse(length(years) <= 4, length(years), ceiling(length(years)/2))
  

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


#------------------------------- Tethys Output -------------------------------#
# wdirr_data_reference <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Reference/wdirr_km3peryr.csv', header = TRUE)
# wdtotal_data_reference <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Reference/wdtotal_km3peryr.csv', header = TRUE)
# wdtotal_data_impacts <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Impacts/wdtotal_km3peryr.csv', header = TRUE)
# wdtotal_data_policy <- data.table::fread('E:/NEXO-UA/Tethys/example/Output/gcam_5p1_ref_MIROC-ESM-CHEM_rcp6p0_2005-2050_Policy/wdtotal_km3peryr.csv', header = TRUE)
# years <- sprintf('%s', seq(from = 2005, to = 2050, by = 5))


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
  
  nameAppend <- paste('_Tethys', scenario, sep = '_')
  
  demand <- sp_mapping(wdtotal_reference, shape, shape_ext, shape_boundary, 'vol',
                      'subRegion', 'hydroshed', nameAppend, main_folder, years)
  return(demand)
}
demand_reference <- tethys_mapping('wdtotal', years, boundary, 'Reference', gcm, rcp)
demand_impacts <- tethys_mapping('wdtotal', years, boundary, 'Impacts', gcm, rcp)
demand_policy <- tethys_mapping('wdtotal', years, boundary, 'Policy', gcm, rcp)




#------------------------------- Xanthos Output -------------------------------#
# for original Xanthos output
run <- 'clim_impacts'
runoff_var_name <- 'q_km3peryear'
# gcm <- 'MIROC-ESM-CHEM'
# rcp <- 'rcp6p0'
time_scale <- '1950_2099'

xanthos_folder <- paste(run, gcm, rcp, sep = '_')
xanthos_output <- paste('E:/NEXO-UA/Xanthos/example/output', xanthos_folder, sep = '/')
runoff_file <- paste0(paste(runoff_var_name, gcm, rcp, time_scale, sep = '_'), '.csv')


# Read Xanthos runoff output in km3
runoff <- data.table::fread(file = paste(xanthos_output, runoff_file, sep = '/'), header = TRUE)
coord_0p5deg <- data.table::fread(file = 'coordinates.csv', header = FALSE) %>% 
  rename(gridcode = V1,
         longitude = V2,
         latitude = V3,
         x = V4,
         y = V5)

runoff_years <- sprintf('%s', seq(from = 1950, to = 2099, by = 1))
runoff_argentina <- runoff %>% 
  mutate(lat = coord_0p5deg$latitude,
         lon = coord_0p5deg$longitude) %>% 
  gather(key = 'year', value = 'value', all_of(runoff_years)) %>% 
  filter(year %in% years, lat <= boundary@ymax &
           lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin) %>% 
  dplyr::select(-id)



supply <- sp_mapping(runoff_argentina, shape, shape_ext, shape_boundary, 'vol',
                     'subRegion', 'hydroshed', '_Xanthos', main_folder, years)


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
    rename(gridcode = V1,
           longitude = V2,
           latitude = V3,
           x = V4,
           y = V5)
  
  years_hist <- c('2005', '2010')
  runoff_years <- sprintf('%s', seq(from = 1970, to = 2010, by = 1))
  runoff_argentina <- runoff %>% 
    mutate(lat = coord_0p5deg$latitude,
           lon = coord_0p5deg$longitude) %>% 
    gather(key = 'year', value = 'value', all_of(runoff_years)) %>% 
    filter(year %in% c('2005', '2010'), lat <= boundary@ymax & lat >= boundary@ymin & lon <= boundary@xmax & lon >= boundary@xmin) %>% 
    dplyr::select(-id)
  
  
  supply_hist <- sp_mapping(runoff_argentina, shape, shape_ext, shape_boundary, 'vol',
                       'subRegion', 'hydroshed', '_XanthosHist', 'Argentina', years_hist)
}


if(F){
  # Aggreagte gridcells to selected polygons
  poly_table_runoff <- metis.grid2poly(gridFiles = runoff_argentina,
                                       subRegShape = shape,
                                       aggType = 'vol', # Xanthos output in km3, so use 'vol' to aggregate
                                       subRegCol = 'subRegion',
                                       subRegType = 'hydroshed',
                                       nameAppend = '_ArgentinaXanthos',
                                       folderName = 'Argentina',
                                       saveFiles = T)
  
  # Find missing polygons
  missing <- dplyr::setdiff(shape$subRegion, poly_table$subRegion) %>% 
    as.data.frame() %>% 
    rename(subRegion = '.') %>% 
    left_join(shape_area, by = 'subRegion')
  if(length(missing$subRegion) > 0){
    rep <- length(years) * length(missing$subRegion)
    n_miss <- length(missing$subRegion)
    missing <- data.frame(subRegion = rep(missing$subRegion, times = length(years)),
                          value = rep(0, times = rep),
                          year = rep(years, each = n_miss),
                          class = rep('class', times = rep),
                          scenario = rep('scenario', times = rep),
                          x = paste('X', rep(years, each = n_miss), sep = ''),
                          param = rep('param', times = rep),
                          aggType = rep('vol', times = rep),
                          subRegType = rep('hyroshed', times = rep))
  }
  
  
  # Add missing polygons and convert unit from volume to depth
  poly_table_supply <- poly_table_runoff %>% 
    # dplyr::union(missing) %>% 
    left_join(shape_area, by = 'subRegion') %>% 
    # filter(subRegion %in% m_argentina_hydroshed@data$subRegion & !subRegion %in% m_chile_hydroshed$subRegion) %>% 
    # dplyr::select(subRegion %in% m_chile_hydroshed$subRegion) %>% 
    mutate(value = (value/area)*1000*1000) # convert from km to mm
  
  # Plot maps
  metis.mapsProcess(polygonTable = poly_table_supply,
                    subRegShape = shape_ext,
                    subRegCol = 'subRegion',
                    subRegType = 'hydroshed',
                    nameAppend = '_Xanthos',
                    folderName = 'Argentina',
                    mapTitleOn = F,
                    legendFixedBreaks = 7,
                    # scaleRange = c(0,254),
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    extension = T,
                    expandPercent = 25,
                    cropToBoundary = F,
                    classPalette = 'pal_wet',
                    facetCols = 5
  )
}



#------------------------------- Water Scarcity -------------------------------#
# demand_clean <- demand$shapeTbl %>% 
#   dplyr::select(-units, -region, -classPalette, -multiFacetCol, -multiFacetRow)
# supply_clean <- supply$shapeTbl %>% 
#   dplyr::select(-units, -region, -classPalette, -multiFacetCol, -multiFacetRow)

scarcity_mapping <- function(demand, supply, scenario){
  poly_table_scarcity <- demand %>% 
    left_join(supply, by = c('subRegion',
                             'year',
                             'x',
                             'class',
                             'scenario',
                             'param',
                             'aggType',
                             'subRegType',
                             'area')) %>% 
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
                    subRegShape = shape_ext,
                    subRegCol = 'subRegion',
                    subRegType = 'hydroshed',
                    nameAppend = nameAppend,
                    folderName = main_folder,
                    mapTitleOn = F,
                    # legendFixedBreaks = 4,
                    # scaleRange = c(0.4, 0),
                    boundaryRegShape = shape_boundary,
                    extendedLabels = T,
                    extension = T,
                    # expandPercent = 25,
                    cropToBoundary = F,
                    # classPalette = 'pal_ScarcityCat',
                    facetCols = 5,
                    animateOn = F,
                    numeric2Cat_list = numeric2Cat_list)
}

scarcity_mapping(demand_reference, supply, 'Reference')
scarcity_mapping(demand_impacts, supply, 'Impacts')
scarcity_mapping(demand_policy, supply, 'Policy')



# sp_scarcity <- sp::SpatialPointsDataFrame(sp::SpatialPoints(coords = (cbind(datax$lon, 
#                                                              datax$lat))), data = datax)
# metis.map(dataPolygon = poly_table_scarcity,
#           fillColumn = 'subRegion',
#           shpFile = shape_ext,
#           fillPalette = 'pal_ScarcityCat',
#           legendBreaks = c(0,0.1,0.2,0.4),
#           fileName = 'Argentina_Scarcity')



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





#------------------------------- NOT IN USE -------------------------------#
if(F){
  
  basins <- unique(m_argentina_basin@data$subRegion)
  value_col <- df_data_by_basin %>% 
    filter(basin_name %in% basins & year %in% 2005 & scenario %in% 'reference') %>% 
    select(value=corn_irr)
  
  
  data <- data.frame(subRegion = value_col$basin_name,
                     year = rep(2005, 12),
                     value = value_col$value)
  
  
  
  # metis.prepGrid (Fail)
  # This function is designed to be used with specific open-source downscaling models (Xanthos [18],
  # Demeter [19], and Tethys [20]) that downscale GCAM data to the grid level. The function takes
  # outputs from these various models and processes them into the format required for providing
  # input to the metis.mapsProcess.R function
  metis.prepGrid(demeterFolders = paste(demeter_output, '/reference_2020-10-05_15h02m50s',sep = ''),
                 demeterScenarios = 'reference',
                 demeterTimesteps = seq(from = 2005, to = 2050, by = 5),
                 demeterUnits = 'fraction',
                 dirOutputs = paste(getwd(), '/demeter', sep = ''))
  
  # metis.grid2poly (Works)
  # Function used to crop and aggregate gridded data by a given polygon shape file. If no grid is
  # provided, the function can still be used to produce regional and subregional maps.
  
  metis.grid2poly(gridFiles = 'E:/NEXO-UA/Results/metis/outputs/Argentina/Grid2Poly/demeter_test.csv',
                  subRegShape = m_argentina_hydroshed,
                  subRegCol = 'subRegion',
                  subRegType = 'hydrosheds',
                  aggType = 'depth',
                  folderName = 'Argentina',
                  saveFiles = T)
  
  data_hydroshed <- data.table::fread('E:/NEXO-UA/Results/metis/outputs/Argentina/Grid2Poly/poly_scenario_hydrosheds_param.csv',
                                      header = TRUE)
  
  str(data_hydroshed)
  col.names <- c('value', 'forest',
                 'shrub', 'grass', 'urban', 'snow', 'sparse', 'corn_irr',
                 'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
                 'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
                 'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
                 'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
                 'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
                 'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
                 'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd')
  
  data_hydroshed[, (col.names) := lapply(.SD, as.numeric), .SDcols = col.names]
  
  
  data <- data_hydroshed %>% 
    group_by(subRegion) %>% 
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
    ungroup()
  
  data <- data.frame(subRegion = as.character(m_argentina_hydroshed$subRegion),
                     year = rep(2005, length(m_argentina_hydroshed$subRegion)),
                     value = m_argentina_hydroshed$SUB_AREA)
  
  metis.mapsProcess(polygonTable = data_hydroshed,
                    subRegShape = m_argentina_hydroshed,
                    folderName = 'Argentina',
                    subRegCol = 'subRegion',
                    # subRegType = 'hydrosheds',
                    mapTitleOn = F,
                    legendFixedBreaks = 6,
                    cropToBoundary = F,
                    extension = T,
                    expandPercent = 15,
                    classPalette = 'pal_green')
  
  
  # read basin codes and coordinates
  basin <- data.table::fread(file = 'basin.csv', header = TRUE); str(basin)
  basin_names <- data.table::fread(file = 'gcam_basin_lookup.csv', header = TRUE) %>% 
    rename(basin = basin_id); str(basin_names)
  coord_5deg <- data.table::fread(file = 'coordinates.csv', header = FALSE) %>% 
    rename(gridcode = V1,
           longitude = V2,
           latitude = V3,
           x = V4,
           y = V5) %>% 
    append(basin) %>% 
    as.data.frame() %>% 
    left_join(basin_names, by = 'basin') %>% 
    select(-x, -y); str(coord_5deg)
  
  # assign the basin codes according to the coordinates
  df_data_2 <- df_data %>% 
    as_tibble() %>% 
    gcamdata::left_join_error_no_match(coord_5deg, by = c('gridcode', 'longitude', 'latitude'))
  
  df_data_by_basin <- df_data_2 %>% 
    group_by(scenario, year, basin, basin_name, glu_name) %>% 
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>% 
    select(-longitude, -latitude, -gridcode)
  str(df_data_by_basin)
  
  # Read in interception of World Basins and 32 region shape file
  map_inter_basin_32reg <- metis::mapIntersectGCAMBasin32Reg; head(map_inter_basin_32reg@data) 
  map_inter_basin_32reg <- map_inter_basin_32reg[map_inter_basin_32reg@data$subRegion_GCAMReg32 %in% c("Argentina"),]
  map_inter_basin_32reg@data <- droplevels(map_inter_basin_32reg@data)
  metis.map(map_inter_basin_32reg,labels=F) # View custom shape
  map_inter_basin_32reg@data <- map_inter_basin_32reg@data %>% 
    dplyr::rename(basins = subRegion_GCAMBasin)
  str(map_inter_basin_32reg@data)
  map_inter_basin_32reg@data <- map_inter_basin_32reg@data %>% 
    dplyr::mutate(subRegion=basins);  map_inter_basin_32reg@data
}