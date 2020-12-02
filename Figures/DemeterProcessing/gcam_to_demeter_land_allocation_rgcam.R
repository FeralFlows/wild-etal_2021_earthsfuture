##################################################################
# Take GCAM query and covert it to Demeter required input format #
# Author: Mengqi Zhao                                            #
# Email: mengqiz@umd.edu                                         #
# Last Update: 2020-10-02                                        #  
##################################################################

rm(list = ls())

if('devtools' %in% rownames(installed.packages()) == F){install.packages('devtools')}
library(devtools)
if('rgcam' %in% rownames(installed.packages()) == F){devtools::install_github(repo='JGCRI/rgcam')}
library(rgcam)
if('dplyr' %in% rownames(installed.packages()) == F){install.packages('dplyr')}
library(dplyr)
if('tidyr' %in% rownames(installed.packages()) == F){install.packages('tidyr')}
library(tidyr)
if('data.table' %in% rownames(installed.packages()) == F){install.packages('data.table')}
library(data.table)
if('DataCombine' %in% rownames(installed.packages()) == F){install.packages('DataCombine')}
library(DataCombine)
if('stringr' %in% rownames(installed.packages()) == F){install.packages('stringr')}
library(stringr)

gcam_path <- 'E:/NEXO-UA/GCAM-Workspace/gcam-core_LAC_v02_5Nov2019/output'
process_path <- 'E:/NEXO-UA/Results/downscaling/land/demeter'

conn <- localDBConn(paste(gcam_path, 'FinalRuns', sep='/'), 'IDBNexus')
prj <- addScenario(conn, proj=paste(process_path, 'IDBNexus.dat', sep='/'), scenario=c('Reference', 'Impacts', 'Policy'), queryFile=paste(process_path, 'query_demeter_33regions_3scenarios.xml', sep='/'))

scenarios <- listScenarios(prj)

# Read global land use mapping file
basin_to_country_mapping <- data.table::fread('E:/NEXO-UA/GCAM-Workspace/gcam-core_LAC_v02_5Nov2019/input/gcamdata/inst/extdata/water/basin_to_country_mapping.csv', header = TRUE)
GCAM_region_names <- data.table::fread('E:/NEXO-UA/GCAM-Workspace/gcam-core_LAC_v02_5Nov2019/input/gcamdata/inst/extdata/common/GCAM_region_names.csv', header = TRUE)

# Convert from list to data frame
gcam_df <- data.frame(t(sapply(prj, c)))
# Create data frame for each scenario 
gcam_impact <- gcam_df[[1]][[1]] %>% rename(glu_name = 'land-region', mgmt_tech = 'mgmt-tech') %>% filter(year<=2050) %>% tidyr::spread(year, value)
gcam_policy <- gcam_df[[2]][[1]] %>% rename(glu_name = 'land-region', mgmt_tech = 'mgmt-tech') %>% filter(year<=2050) %>%tidyr::spread(year, value)
gcam_ref <- gcam_df[[3]][[1]] %>% rename(glu_name = 'land-region', mgmt_tech = 'mgmt-tech') %>% filter(year<=2050) %>%tidyr::spread(year, value)

# Processing query
query_processing <- function(gcam, gcm, rcp){
  scenario_name <- unique(gcam$scenario)
  
  gcam <- gcam %>% 
    group_by(region, glu_name, crop, water) %>% 
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
    unite(landclass, crop, water) %>% 
    as.data.frame()
  
  rplc <- data.frame(from = c('_NA'), to = c(''))
  gcam <- FindReplace(data = gcam, Var = 'landclass', replaceData = rplc, from = 'from', to = 'to', exact = FALSE)
  
  # Add GLU code (same with basin id) based on basin to country mapping
  basin_to_country_mapping %>%
    mutate(GLU_code = str_remove(sub('\\GLU', '', basin_to_country_mapping$GLU_code), '^0+')) %>%
    select(-GCAM_basin_ID, -Basin_name, -ISO, -ISO_NUM, -Country_name) %>% 
    rename(glu_name = GLU_name,
           metric_id = GLU_code) %>%
    right_join(gcam, by = 'glu_name') ->
    gcam_glu
  
  # Add country code (gcam_region) based on GCAM_region_names
  GCAM_region_names %>% 
    rename(region_id = GCAM_region_ID) %>% 
    right_join(gcam_glu, by = 'region') %>% 
    as.data.frame() ->
    gcam_glu_region 
  
  
  # Reorder the columns in the output
  gcam_years <- sprintf('%s',seq(from = 2005, to = 2050, by = 5))
  col_order <- c('region', 'glu_name', 'region_id', 'metric_id', 'landclass', '1975', '1990', gcam_years)
  gcam_output <- gcam_glu_region[, col_order] %>% 
    arrange(region, glu_name)
  
  output_name <- paste('DemeterDownscaled_33Regions', gcm, rcp, sep = '_')
  save_path <- 'E:/NEXO-UA/Demeter/example/inputs/projected/'
  write.csv(gcam_output, file = paste(save_path, output_name, scenario_name, '.csv', sep = ''), row.names = FALSE)
}

query_processing(gcam_ref, 'MIROC-ESM-CHEM', 'rcp6p0')
query_processing(gcam_impact, 'MIROC-ESM-CHEM', 'rcp6p0')
query_processing(gcam_policy, 'MIROC-ESM-CHEM', 'rcp6p0')
