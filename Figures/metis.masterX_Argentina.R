
# metis.master.R
# Script to run different parts of the metis package.

rm(list=ls())

#----------------------------
# Install necessary packages
#----------------------------
if("devtools" %in% rownames(installed.packages()) == F){install.packages("devtools")}
library(devtools)
if("metis" %in% rownames(installed.packages()) == F){install_github(repo="JGCRI/metis")}
library(metis)
if("rgcam" %in% rownames(installed.packages()) == F){install_github(repo="JGCRI/rgcam")}
library(rgcam)
if("jgcricolors" %in% rownames(installed.packages()) == F){install_github(repo="JGCRI/jgcricolors")}
library(jgcricolors)
if("tibble" %in% rownames(installed.packages()) == F){install.packages("tibble")}
library(tibble)
if("ggplot2" %in% rownames(installed.packages()) == F){install.packages("ggplot2")}
library(ggplot2)
if("zoo" %in% rownames(installed.packages()) == F){install.packages("zoo")}
library(zoo)
if("dplyr" %in% rownames(installed.packages()) == F){install.packages("dplyr")}
library(dplyr)
if("dbplyr" %in% rownames(installed.packages()) == F){install.packages("dbplyr")}
library(dbplyr)
if("tidyr" %in% rownames(installed.packages()) == F){install.packages("tidyr")}
library(tidyr)
if("stringr" %in% rownames(installed.packages()) == F){install.packages("stringr")}
library(stringr)

library(plutus)


#----------------------------
# Read GCAM Data (metis.readgcam.R)
#---------------------------

setwd('E:/NEXO-UA/Results/metis/gcam_database')

# Connect to gcam database or project
gcamdatabasePath_i <- 'E:/NEXO-UA/GCAM-Workspace/gcam-core_LAC_v02_5Nov2019/output/FinalRuns' # Use if gcamdatabase is needed
gcamdatabaseName_i <- 'IDBNexus_MIROC-ESM-CHEM_rcp6p0' # "Reference_originalSW" Use if gcamdatabse is needed
dataProjPath_i <- getwd() # Path to dataProj file.
dataProj_i <-"dataProj_gcam5p1_MIROC-ESM-CHEM_rcp6p0.proj"  # Use if gcamdata has been saved as .proj file

# Get list of scenarios and rename if desired.
# conn <- rgcam::localDBConn(gcamdatabasePath_i,gcamdatabaseName_i) # if connecting directly to gcam database
# prj <- rgcam::addScenario(conn, proj=paste(dataProjPath_i, "IDBNexus.proj", sep='/'), scenario=c('Reference', 'Impacts', 'Policy'))
# dataProjLoaded <- loadProject(paste(dataProjPath_i, "/",dataProj_i , sep = ""))
#  listScenarios(dataProjLoaded)  # List of Scenarios in GCAM database

scenOrigNames_i = c("Reference", "Impacts", "Policy") #  c(, "DeepDecarb1Mkt_2DS")  #'GCAMOriginal',
scenNewNames_i = c("Reference", "Climate Impacts", "Climate Policy") # "Original",   c(, "DeepDecarb1Mkt_2DS"")  # Names to replace the original names for final figures.

# Choose Parameters or set to "All" for all params. For complete list see ?metis.readgcam
paramsList <-  metis.mappings()$mapParamQuery
# paramsSelect_i <- c('aggLandAlloc')
# paramsSelect_i <- c('watWithdrawBySec')
# paramsSelect_i <- c('watWithdrawBySec',
#                     'watWithdrawByCrop',
#                     'watIrrWithdrawBasin',
#                     'landAlloc',
#                     'landAllocByCrop',
#                     'elecByTechTWh',
#                     'elecNewCapCost',
#                     "elecCapByFuel",
#                     "elecFinalBySecTWh",
#                     "elecFinalByFuelTWh",
#                     "elecNewCapGW",
#                     "elecAnnualRetPrematureCost",
#                     "elecAnnualRetPrematureGW",
#                     "elecCumCapCost",
#                     "elecCumCapGW",
#                     "elecCumRetPrematureCost",
#                     "elecCumRetPrematureGW",
#                     'energyPrimaryByFuelEJ',
#                     "energyPrimaryRefLiqProdEJ",
#                     'energyFinalConsumBySecEJ',
#                     "energyFinalByFuelBySectorEJ",
#                     "energyFinalSubsecByFuelTranspEJ",
#                     "energyFinalSubsecByFuelBuildEJ",
#                     "energyFinalSubsecByFuelIndusEJ",
#                     "energyFinalSubsecBySectorBuildEJ",
#                     'agProdByCrop'
#                     )


if(T){
  ## If selecting all parameters, switch to if(T)
  paramesSelect_energy <- c("energyPrimaryByFuelEJ","energyPrimaryRefLiqProdEJ", "energyFinalConsumBySecEJ","energyFinalByFuelBySectorEJ","energyFinalSubsecByFuelTranspEJ", "energyFinalSubsecByFuelBuildEJ", "energyFinalSubsecByFuelIndusEJ","energyFinalSubsecBySectorBuildEJ", "energyPrimaryByFuelMTOE","energyPrimaryRefLiqProdMTOE", "energyFinalConsumBySecMTOE","energyFinalbyFuelMTOE","energyFinalSubsecByFuelTranspMTOE", "energyFinalSubsecByFuelBuildMTOE", "energyFinalSubsecByFuelIndusMTOE","energyFinalSubsecBySectorBuildMTOE", "energyPrimaryByFuelTWh","energyPrimaryRefLiqProdTWh", "energyFinalConsumBySecTWh","energyFinalbyFuelTWh","energyFinalSubsecByFuelTranspTWh", "energyFinalSubsecByFuelBuildTWh", "energyFinalSubsecByFuelIndusTWh","energyFinalSubsecBySectorBuildTWh")

  paramsSelect_electricity <- c("elecByTechTWh","elecCapByFuel","elecFinalBySecTWh","elecFinalByFuelTWh", "elecNewCapCost","elecNewCapGW","elecAnnualRetPrematureCost","elecAnnualRetPrematureGW","elecCumCapCost","elecCumCapGW","elecCumRetPrematureCost","elecCumRetPrematureGW")

  paramsSelect_transport <- c( "transportPassengerVMTByMode", "transportFreightVMTByMode", "transportPassengerVMTByFuel", "transportFreightVMTByFuel")

  paramsSelect_water <- c("watConsumBySec", "watWithdrawBySec", "watWithdrawByCrop", "watBioPhysCons", "watIrrWithdrawBasin","watIrrConsBasin")

  paramsSelect_socioecon <- c("gdpPerCapita", "gdp", "gdpGrowthRate", "pop")

  paramsSelect_ag <- c("agProdbyIrrRfd", "agProdBiomass", "agProdForest","agProdByCrop")

  paramsSelect_livestock <- c("livestock_MeatDairybyTechMixed","livestock_MeatDairybyTechPastoral","livestock_MeatDairybyTechImports", "livestock_MeatDairybySubsector")

  paramsSelect_land <- c("landIrrRfd", "landIrrCrop","landRfdCrop", "landAlloc","landAllocByCrop")

  paramsSelect_emissions <- c("emissLUC", "emissNonCO2BySectorGWPAR5","emissNonCO2BySectorGTPAR5","emissNonCO2BySectorOrigUnits", "emissNonCO2ByResProdGWPAR5", "emissBySectorGWPAR5FFI","emissMethaneBySourceGWPAR5", "emissByGasGWPAR5FFI", "emissByGasGWPAR5LUC", "emissBySectorGWPAR5LUC", "emissNonCO2ByResProdGTPAR5", "emissBySectorGTPAR5FFI","emissMethaneBySourceGTPAR5", "emissByGasGTPAR5FFI", "emissByGasGTPAR5LUC","emissBySectorGTPAR5LUC", "emissCO2BySectorNoBio")

  paramsSelect_all <- c(paramesSelect_energy, paramsSelect_electricity, paramsSelect_transport, paramsSelect_water, paramsSelect_socioecon, paramsSelect_ag, paramsSelect_livestock, paramsSelect_land, paramsSelect_emissions)

  paramsSelect_i <- paramsSelect_all
}


queriesSelect_i <- c("All")

# Select regions from the 32 GCAM regions.
regionsSelect_i <- c("Argentina")

# Read GCAM database
if(file.exists(paste(dataProjPath_i, 'outputs', dataProj_i, sep = "/"))){
  dataGCAM <- metis.readgcam(dataProjFile = paste(dataProjPath_i, 'outputs', dataProj_i, sep = "/"),
                                 scenOrigNames = scenOrigNames_i,
                                 scenNewNames = scenNewNames_i,
                                 regionsSelect = regionsSelect_i,
                                 paramsSelect = paramsSelect_i)
}else{
  dataGCAM <- metis.readgcam(gcamdatabase = paste(gcamdatabasePath_i, gcamdatabaseName_i, sep='/'),
                                 reReadData = TRUE,
                                 scenOrigNames = scenOrigNames_i,
                                 scenNewNames = scenNewNames_i,
                                 regionsSelect = regionsSelect_i,
                                 paramsSelect = paramsSelect_i,
                                 queryFile = NULL,
                                 nameAppend = gcamdatabaseName_i)
  file.rename(paste(getwd(),'outputs', 'dataProj.proj', sep = '/'),
              paste(getwd(), 'outputs', dataProj_i, sep = '/'))
}

# Calculate electricity investments
if(file.exists(paste(dataProjPath_i, 'outputs', paste0(gsub('.proj', '', dataProj_i), '_invest.proj'), sep = "/"))){
  invest <- plutus::gcamInvest(dataProjFile = paste(dataProjPath_i, 'outputs', paste0(gsub('.proj', '', dataProj_i), '_invest.proj'), sep = "/"),
                               reReadData = T,
                               scenOrigNames = scenOrigNames_i,
                               scenNewNames = scenNewNames_i,
                               regionsSelect = regionsSelect_i,
                               saveData = F)
}else{
  invest <- plutus::gcamInvest(gcamdatabase = paste(gcamdatabasePath_i, gcamdatabaseName_i, sep = '/'),
                               reReadData = T,
                               scenOrigNames = scenOrigNames_i,
                               scenNewNames = scenNewNames_i,
                               regionsSelect = regionsSelect_i,
                               saveData = F)
  file.rename(paste(getwd(),'outputs', 'dataProj.proj', sep = '/'),
              paste(getwd(), 'outputs', paste0(gsub('.proj', '', dataProj_i), '_invest.proj'), sep = '/'))
}

dataGCAM$data # To view the data read that was read.
invest$data
#------------------------------------------------------------------------------------------
# Charts Process (metis.chartsProcess.R)
#------------------------------------------------------------------------------------------

# Can also add data .csv outputs from metis.readgcam.R which are autmatically saved in
# ./metis/outputs/readGCAMTables/Tables_gcam
# for each of the regions selected.
# gcamDataTable_Argentina.csv, gcamDataTable_China.csv, gcamDataTable_Pakistan.csv
# This would be added to dataTables_i as:
# dataTables_i = c(paste(getwd(), "/outputs/readGCAMTables/Tables_local/local_Regional_Colombia.csv", sep = "")
#                  #paste(getwd(), "/outputs/readGCAMTables/Tables_gcam/gcamDataTable_Colombia.csv", sep = "")
# )

# Read in the data from the function metis.readgcam.
rTable_i <- rbind(dataGCAM$data, invest$data)
rTable_i <- rTable_i[complete.cases(rTable_i$class1),] # remove rows with NA values

# Get classes that have color codes in metis color palette pal_metis
class_metis <- unlist(attributes(metis.colors()$'pal_metis'), use.names=FALSE)
# Find parameters no have all classes in calss_metis, but the assigned palette is pal_metis
param_pal_16 <- unique(rTable_i$param[!(rTable_i$class1 %in% class_metis) & rTable_i$classPalette1 == 'pal_metis'])
param_pal_16_select <- c('livestock_MeatDairybyTechMixed', 'livestock_MeatDairybySubsector', 'livestock_MeatDairybyTechPastoral', 'watConsumBySec', 'transportPassengerVMTByMode', 'transportFreightVMTByMode', 'watConsumBySec')
rTable_i$classPalette1[rTable_i$param %in% param_pal_16_select] <- 'pal_16'

rTable_i$classPalette1[rTable_i$param %in% c('watWithdrawBySec', 'watBioPhysCons')] <- 'pal_metis'

# For GCAM 5.1
rTable_i$class1[grepl('biomass_tree|biomass_grass', rTable_i$class1)] <- 'biomass'

# For GCAM 5.3
rTable_i$class1[grepl('biomassTree|biomassGrass', rTable_i$class1)] <- 'biomass'
rTable_i$class1[grepl('naturalOtherTree|naturalOtherGrass', rTable_i$class1)] <- 'biomass'
rTable_i$class1[grepl('^landAlloc$', rTable_i$param) & grepl('RootTuber|biomass', rTable_i$class1)] <- 'crops'


charts <- metis.chartsProcess(rTable = rTable_i, # Default is NULL
                              #dataTables=dataTables_i, # Default is NULL
                              paramsSelect = paramsSelect_i, # Default is "All"
                              regionsSelect = regionsSelect_i, # Default is "All"
                              xCompare = c("2010", "2030", "2050"), # Default is c("2015","2030","2050","2100")
                              scenRef = "Reference", # Default is NULL
                              dirOutputs = paste(getwd(), "/outputs", sep = ""), # Default is paste(getwd(),"/outputs",sep="")
                              regionCompareOnly = 0, # Default 0. If set to 1, will only run comparison plots and not individual
                              scenarioCompareOnly = 0,
                              regionCompare = 0,
                              useNewLabels = 0,
                              folderName = "IDBNexus",
                              xRange = c(2020, 2030, 2040, 2050),
                              colOrder1 = c("Reference", "Climate Impacts", "Climate Policy"), #"Original",
                              colOrderName1 = "scenario",
                              pdfpng = 'pdf',
                              # plotBGcolor = 'white',
                              multiPlotFigsOnly = T) # Default 0. If set to 1, will only run comparison plots and not individual


if(T){
  # Only for when rTable = invest$data, and custom orders for Fuel in the plot
  df_plot <- invest$data %>%
    dplyr::select(scenario, region, subRegion, param, class1, x, units, value) %>%
    dplyr::filter(x %in% seq(2015, 2050, 5)) %>%
    dplyr::rename(sector = class1,
                  year = x)
  df_plot$sector[grepl('^Bioenergy$', df_plot$sector)] <- 'Biomass'
  df_plot$sector[grepl('^Bioenergy CCS$', df_plot$sector)] <- 'Biomass CCS'

  if(T){
    # Calculate 5 year moving average. For example, value in 2020 is the average value
    # between 2015, 2020, and 2025
    df_plot <- df_plot %>%
      dplyr::group_by(scenario, param, sector) %>%
      dplyr::mutate(value_15 = zoo::rollmean(value, k = 3, align = 'center', fill = NA),
                    value_10 = zoo::rollmean(value, k = 2, align = 'right', fill = NA),
                    value = if_else(year == 2050, value_10, value_15)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-value_15, -value_10)
  }

  # Filter GCAM output for elecByTechTWh
  elec <- dataGCAM$data %>%
    dplyr::select(scenario, region, subRegion, param, class1, class2, x, units, value) %>%
    dplyr::filter(x %in% seq(2015, 2050, 5), param %in% 'elecByTechTWh') %>%
    dplyr::rename(sector = class1,
                  technology = class2,
                  year = x) %>%
    dplyr::mutate(sector = stringr::str_to_title(sector),
                  sector = if_else(grepl('Refined Liquids', sector), 'Oil', sector),
                  agg_tech = if_else(grepl('CCS', technology), paste0(sector, ' CCS'), sector)) %>%
    dplyr::select(-sector, -technology) %>%
    dplyr::rename(sector = agg_tech) %>%
    dplyr::group_by(scenario, region, subRegion, param, units, sector, year) %>%
    dplyr::summarise(value = sum(value))

  df_plot <- rbind(df_plot, elec)
  df_plot <- df_plot[complete.cases(df_plot$sector),]

  ## Reference, Difference of Impacts-Reference, Policy-Reference
  value_ref <- df_plot %>%
    dplyr::filter(scenario %in% 'Reference') %>%
    dplyr::rename(value_ref = value) %>%
    dplyr::select(-scenario)
  df_plot_diff <- df_plot %>%
    dplyr::left_join(value_ref, by = c('region', 'subRegion', 'param', 'sector', 'units', 'year')) %>%
    dplyr::mutate(value_ref = if_else(is.na(value_ref), 0, value_ref),
                  value = if_else(scenario == 'Reference', value, value - value_ref),
                  value = if_else(is.na(value), 0, value)) %>%
    dplyr::select(-value_ref)

  ## Reorder scenarios and sectors
  df_plot_diff$scenario[grepl('Impacts', df_plot_diff$scenario)] <- 'Climate Impacts_Diff'
  df_plot_diff$scenario[grepl('Policy', df_plot_diff$scenario)] <- 'Climate Policy_Diff'

  df_plot$scenario <- factor(df_plot$scenario, levels = c('Reference', 'Climate Impacts', 'Climate Policy'))
  df_plot_diff$scenario <- factor(df_plot_diff$scenario, levels = c('Reference', 'Climate Impacts_Diff', 'Climate Policy_Diff'))

  sector_order <- c("Geothermal", "Solar", "Wind", "Hydro", "Nuclear", "Biomass CCS", "Biomass",
                   "Gas CCS", "Gas", "Oil", "Coal CCS", "Coal")
  df_plot$sector <- factor(df_plot$sector, levels = sector_order)
  df_plot_diff$sector <- factor(df_plot_diff$sector, levels = sector_order)

  ## Plot
  pal_gcam <- jgcricol()$pal_all
  # pal_gcam['Biomass'] <- '#BED96A'

  plot_ag <- function(df_plot, save_name, ts){
    units <- unique(df_plot$units)
    ggplot(df_plot, aes(x = year, y = value, fill = sector)) +
      geom_bar(position = "stack", stat = "identity", col = "black", lwd = 0.4) +
      scale_x_discrete(expand = c(0.13,0), limits = ts) +
      facet_grid(.~scenario) +
      scale_fill_manual(values=pal_gcam) +
      theme_bw() +
      xlab(NULL) +
      ylab(units) +
      theme(axis.text = element_text(size = 12, color = 'black'),
            axis.title = element_text(size = 14),
            panel.grid.minor.x = element_blank(),
            strip.text.x = element_text(size = 18, color = 'white'),
            strip.background = element_rect(fill = '#394546'),
            legend.title = element_text(size = 14),
            # legend.position = 'bottom',
            legend.text = element_text(size = 12),
            legend.key.size = unit(0.4, 'cm'))

    ggsave(save_name, height = 5, width = 14, unit = 'in', dpi = 600)
  }

  ts <- seq(2020, 2050, 10)
  for(param_i in unique(df_plot$param)){
    plot <- df_plot %>%
      filter(param %in% param_i, year %in% ts)
    plot_ag(df_plot = plot,
            save_name = paste0(
              getwd(),
              '/outputs/IDBNexus/Charts/Argentina_gcam5p1_MIROC-ESM-CHEM_rcp6p0_PaperReview/',
              param_i,
              '_Argentina_rollmean.pdf'
            ),
            ts = ts)
  }

  for(param_i in unique(df_plot_diff$param)){
    plot <- df_plot_diff %>%
      filter(param %in% param_i, year %in% ts)
    plot_ag(df_plot = plot,
            save_name = paste0(
              getwd(),
              '/outputs/IDBNexus/Charts/Argentina_gcam5p1_MIROC-ESM-CHEM_rcp6p0_PaperReview/',
              param_i,
              '_Argentina_Diff_rollmean.pdf'
            ),
            ts = ts)
  }

}
