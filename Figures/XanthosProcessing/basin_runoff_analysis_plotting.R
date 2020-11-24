# Purpose: Plots Climate Impacts on Water Availability by Basin from Xanthos

library("dplyr")
library("tidyr")
library("ggplot2")
library("readr")
library("zoo")
library(scales)
library(stats)
library(magrittr)
library(data.table)
library(gcamdata)

source('E:/NEXO-UA/Xanthos/example/xanthos_postprocessing_fns_MZ.R')
setwd('E:/NEXO-UA/Results/downscaling/water/xanthos')
rcp_colors <- c("rcp8p5" = "#736F6E",
                "rcp6p0" = "#C0C0C0",
                "rcp4p5" = "#98AFC7",
                "rcp2p6" = "#6698FF",
                "historical" = 'black',
                'historical mean' = 'black')
gcm_colors <- c("NorESM1-M" = "#736F6E",
                "MIROC-ESM-CHEM" = "#C0C0C0",
                "IPSL-CM5A-LR" = "#98AFC7",
                "HadGEM2-ES" = "#6698FF",
                "GFDL-ESM2M" = "#153E7E",
                "watch+wfdei" = 'black',
                'historical' = 'black',
                'historical mean' = 'black')
figures_basepath <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output/figures/basin'
results_basepath <- 'E:/NEXO-UA/Xanthos/example/output'
xanthos_config_names <- c('clim_impacts')
gcm_names <- c('NorESM1-M', 'MIROC-ESM-CHEM', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'GFDL-ESM2M')
# gcm_names <- c('HadGEM2-ES') # temp MZ
gcm_names_incl_hist <- append(gcm_names, 'watch+wfdei')
rcp_names <- c('rcp2p6', 'rcp4p5', 'rcp6p0', 'rcp8p5')
# rcp_names <- c('rcp6p0') # temp MZ
rcp_names_incl_hist <- append(rcp_names, "historical")
xanthos_var_names <- c('Basin_runoff_km3peryear')
time_scale <- '1950_2099'
# water_basins_plot <- c("La Plata", 'Caribbean Coast', 'Magdalena', 'Orinoco', 'Amazon',
#                   'Colombia - Ecuador Pacific Coast', 'Mar Chiquita', 'La Puna Region', 'Salinas Grandes',
#                   'Pampas Region', 'North Argentina South Atlantic Coast', 'South America Colorado', 'Negro',
#                   'Central Patagonia Highlands', 'South Argentina South Atlantic Coast', 'South Chile Pacific Coast', 'Uruguay - Brazil South Atlantic Coast')
water_basins_plot <- c("La Plata", 'La Puna Region', 'Salinas Grandes', 'Mar Chiquita', 'South America Colorado', 'Pampas Region', 'North Argentina South Atlantic Coast', 'Negro', 'Central Patagonia Highlands', 'South Chile Pacific Coast', 'South Argentina South Atlantic Coast'
)

filter_list_2 <- list("Basin_runoff_km3peryear" = water_basins_plot)
stored_in_dir <- 1  # = 1 if dragged whole xanthos folder off pic; 0 if just dragged file down into results dir on comp
run_name <- c('clim_impacts')
gcam_years <- c(2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065, 2070, 2075,
                2080, 2085, 2090, 2095, 2100)

# Set base input files directory
level2_out_dir <- 'E:/NEXO-UA/Results/impacts/water_avail'
extras_dir <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output'
delta <- adjust_gcm_mean(results_basepath, extras_dir, level2_out_dir, time_scale, stored_in_dir, run_name,
                         xanthos_var_names)
deltas_gcm_all <- delta$deltas_gcm_all
runoff_gcm_all_GCAM <- delta$runoff_gcm_all_GCAM_2
runoff_gcm_all_GCAM_3 <- delta$runoff_gcm_all_GCAM_3

# Read in xanthos runoff data; separate out historical and projected data
add_historical <- 1  # Add historical values to dataframe

df_all_runs_basin <- xanthos_proc(xanthos_var_names, xanthos_config_names, gcm_names, rcp_names, time_scale,
                                  results_basepath)$output
df_all_runs_basin <- xanthos_hist_proc(xanthos_var_names, xanthos_config_names, df_all_runs_basin, stored_in_dir,results_basepath, add_historical)$output
df_all_runs_basin_hist <- df_all_runs_basin %>% filter(rcp == 'historical')
df_all_runs_basin <- df_all_runs_basin %>% filter(rcp != 'historical')
df_all_runs_basin$year <- as.numeric(df_all_runs_basin$year)
# Adjust by delta factors
df_all_runs_basin_hist$year <- as.numeric(df_all_runs_basin_hist$year)
df_all_runs_basin <- df_all_runs_basin %>%
  left_join(deltas_gcm_all, by=c('name', 'gcm', 'rcp', 'year')) %>%
  mutate(value = delta_factor*value) %>%
  select(-delta_factor)
#df_all_runs_basin_hist <- df_all_runs_basin_hist %>%
#  left_join(deltas_gcm_all, by=c('name', 'gcm', 'rcp', 'year')) %>%
#  mutate(value = delta_factor*value) %>%
#  select(-delta_factor)
df_all_runs_basin$year <- as.numeric(df_all_runs_basin$year)

# Compute rolling mean--BASIN
roll_window <- 1  # Establish target window for rolling mean; k=1 is no rolling mean, just returns regular values.
df_2_all_runs_basin <- roll_mean(df_all_runs_basin, xanthos_var_names, xanthos_config_names, gcm_names_incl_hist,
                                 rcp_names_incl_hist, k=roll_window)$output
df_2_all_runs_basin$year <- as.numeric(df_2_all_runs_basin$year)

# Compute the mean value, store it for every gcm/rcp combo. Also compute the percent change in every year relative to
# that mean, and store that
df_2_all_runs_basin['mean_2010'] <- 0 # add mean_2010 column
country_list_full <- unique(df_2_all_runs_basin$name)
for (reg in country_list_full){
  for (gcm1 in gcm_names){
    for (rcp1 in rcp_names){
      mean_val <- (df_2_all_runs_basin %>% filter(name==reg, gcm==gcm1, rcp==rcp1, year==2010))$smoothedY[1]
      df_2_all_runs_basin <- df_2_all_runs_basin %>% mutate(mean_2010 = if_else(name==reg & gcm==gcm1 & rcp==rcp1, mean_val, mean_2010))
    }
  }
}
df_2_all_runs_basin <- df_2_all_runs_basin %>% mutate(clim_imp_perc = 100*(smoothedY-mean_2010)/mean_2010) %>% select(-mean_2010)

# Compute historical average, and insert that into historical dataframe (df_all_runs_basin_hist)
df_all_runs_basin_hist$mean_hist <- 0
for(reg in country_list_full){
  df_all_runs_basin_hist <- df_all_runs_basin_hist %>%
    mutate(mean_hist=if_else(name==reg, mean((df_all_runs_basin_hist %>% filter(name==reg))$value), mean_hist))
  merge_df_hist <- df_all_runs_basin_hist %>% filter(year==2010) %>% select(name, mean_hist)
}

# Merge historical averages with df_2_all_runs_basin to compute percentage changes that are projected to occur relative
# to historical averages.
df_2_all_runs_basin <- df_2_all_runs_basin %>% left_join(merge_df_hist, by=c('name'))
df_2_all_runs_basin <- df_2_all_runs_basin %>% mutate(clim_imp_val=mean_hist+mean_hist*(clim_imp_perc/100)) %>%
  mutate(clim_imp_val=if_else(clim_imp_val<0,0,clim_imp_val))

#**************************************** PLOT START ******************************************#
# Plot basin runoff (individual plot for each basin specified earlier, including all the RCP and GCAM combinations).
# This is not separated out by country, this is the entire basin's runoff.
y_ax_lbl <- expression(Annual~Runoff~(km^3))
input <- df_2_all_runs_basin
roll <- 0
start_yr <- 2010
end_yr <- 2050  # 2100
var_names <- c('Basin_runoff_km3peryear')
region_list <- water_basins_plot
region_single_plot(var_names, region_list, input, figures_basepath, start_yr, end_yr,
                   gcm_names, rcp_names, roll, y_ax_lbl, trendline=0)

# Plot smoothed individual basin runoff plots.
y_ax_lbl <- expression(Annual~Runoff~(km^3))
input <- df_2_all_runs_basin
roll <- 2
start_yr <- 2010
end_yr <- 2050
var_names <- c('Basin_runoff_km3peryear')
region_list <- water_basins_plot
region_single_plot(var_names, region_list, input, figures_basepath, start_yr, end_yr,
                   gcm_names, rcp_names, roll, y_ax_lbl, trendline=0)

# Create faceted plot across GCMs and RCPs for Basin runoff
roll <- 0
start_yr <- 2010
end_yr <- 2050  # 2100
start_yr_hist <- 1970
end_yr_hist <- 2009
xanthos_var_names <- c('Basin_runoff_km3peryear')
filter_list_2 <- list("Basin_runoff_km3peryear" = water_basins_plot)
for(var_1 in xanthos_var_names){
  for(reg in filter_list_2[[var_1]]){
    fig_name <- paste0(figures_basepath, '/', reg, "_", var_1, "_", 'gcm_rcp_facet.png')
    plot_df <- df_2_all_runs_basin %>% filter(name==reg, year>=start_yr, year<=end_yr, var==var_1) %>%
      filter(gcm %in% gcm_names_incl_hist, rcp %in% rcp_names_incl_hist)
    plot_df_hist <- df_all_runs_basin_hist %>% filter(name==reg, year>=start_yr_hist, year<=end_yr_hist, var==var_1)
    facet_grid_plot(plot_df, fig_name, rolling=roll, y_lbl=y_ax_lbl, df_all_runs_hist=plot_df_hist, historical=1)
  }
}

# Create a smooth faceted plot across GCMs and RCPs for Basin runoff
roll <- 2
# Use all other values (except "roll" from above)
for(var_1 in xanthos_var_names){
  for(reg in filter_list_2[[var_1]]){
    fig_name <- paste0(figures_basepath, '/', reg, "_", var_1, "_", 'gcm_rcp_facet_smooth.png')
    plot_df <- df_2_all_runs_basin %>% filter(name==reg, year>=start_yr, year<=end_yr, var==var_1) %>%
      filter(gcm %in% gcm_names_incl_hist, rcp %in% rcp_names_incl_hist)
    plot_df_hist <- df_all_runs_basin_hist %>% filter(name==reg, year>=start_yr_hist, year<=end_yr_hist, var==var_1)
    facet_grid_plot(plot_df, fig_name, rolling=roll, y_lbl=y_ax_lbl, df_all_runs_hist=plot_df_hist, historical=1)
  }
}

# Plot basin runoff (individual plot for each basin specified earlier, where all the GCM and RCP combinations are
# combined on the same plot
y_ax_lbl <- expression(Annual~Runoff~(km^3))
input <- df_2_all_runs_basin %>% filter(year>=2009, year<=2050)
roll <- 0
start_yr <- 1970
end_yr <- 2050
var_names <- c('Basin_runoff_km3peryear')
region_list <- water_basins_plot
region_single_plot(xanthos_var_names, region_list, input, figures_basepath, start_yr, end_yr, gcm_names, rcp_names,
                   roll, y_ax_lbl, trendline=0, combined_lines=1, plot_df_hist=df_all_runs_basin_hist)

# Plot facet of values smoothed separately
start_yr <- 2010
start_yr_mod <- 2010
end_yr <- 2050
start_yr_hist <- 1970
end_yr_hist <- 2010
end_yr_hist_mod <- 2010
color_palette <- c('gcm_colors')
model <- c('pm_abcd_mrtm')
xanthos_var_names <- c('Basin_runoff_km3peryear')
# Add columns that plotting module expects to find
runoff_gcm_all_GCAM <- runoff_gcm_all_GCAM %>% mutate(FillPalette = color_palette,
                                                                var=xanthos_var_names, mod=model)
runoff_gcm_all_GCAM_3 <- runoff_gcm_all_GCAM_3 %>% mutate(FillPalette = color_palette,
                                                      var=xanthos_var_names, mod=model)
# Break into a historical portion and a future portion, so they can be plotted separately
runoff_gcm_all_GCAM_hist <- runoff_gcm_all_GCAM %>% filter(year>=start_yr_hist, year<=end_yr_hist_mod) %>%
  mutate(rcp='historical mean', gcm='historical mean')
runoff_gcm_all_GCAM_fut <- runoff_gcm_all_GCAM_3 %>% filter(year>=start_yr_mod, year<=end_yr)
# Produce faceted plot
roll <- 0
filter_list_2 <- list("Basin_runoff_km3peryear" = water_basins_plot)
for(var_1 in xanthos_var_names){
  for(reg in filter_list_2[[var_1]]){
    fig_name <- paste0(figures_basepath, '/', reg, "_", var_1, "_", 'gcm_rcp_facet_GCAM.png')
    plot_df <- runoff_gcm_all_GCAM_fut %>% filter(name==reg, year>=start_yr, year<=end_yr, var==var_1) %>%
      filter(gcm %in% gcm_names_incl_hist, rcp %in% rcp_names_incl_hist)
    plot_df_hist <- runoff_gcm_all_GCAM_hist %>% filter(name==reg, year>=start_yr_hist, year<=end_yr_hist, var==var_1)
    facet_grid_plot(plot_df, fig_name, rolling=roll, y_lbl=y_ax_lbl, df_all_runs_hist=plot_df_hist, historical=1)
  }
}

# Plot all of above lines that appear in facet, but all on the same plot
y_ax_lbl <- expression(Annual~Runoff~(km^3))
#input <- runoff_gcm_all_GCAM_fut %>% filter(year>=2010, year<=2050) %>% mutate(smoothedY=clim_imp_val)
start_yr_hist <- 1990
end_yr_hist <- 2010
input <- df_2_all_runs_basin %>% filter(year>=2010, year<=2050) %>% mutate(smoothedY=clim_imp_val)
runoff_gcm_all_GCAM_hist <- df_2_all_runs_basin %>% filter(year>=start_yr_hist, year<=end_yr_hist_mod) %>%
  mutate(rcp='historical mean', gcm='historical mean') %>% mutate(smoothedY=mean_hist)
roll <- 2
start_yr <- 2010
end_yr <- 2050
var_names <- c('Basin_runoff_km3peryear')
# region_list <- c('Negro', 'South America Colorado')
region_list <- water_basins_plot
gcm_list <- 'MIROC-ESM-CHEM' # 'GFDL-ESM2M' # IPSL-CM5A-LR' # 'HadGEM2-ES'
rcp_list <- c('rcp6p0')  # 'rcp2p6'  # 'rcp8p5'
region_single_plot(xanthos_var_names, region_list, input, figures_basepath, start_yr, end_yr, gcm_names, rcp_names,
                   roll, y_ax_lbl, trendline=0, combined_lines=1, plot_df_hist=runoff_gcm_all_GCAM_hist,
                   all_same_color = 1, titles = 'No', legend_on=F, gcm_list=gcm_list, rcp_list=rcp_list,
                   fig_type='.pdf')

# Plot percentage reduction in smoothed runoff compared with 2010
y_ax_lbl <- expression(Change~('%')~'from'~2010~runoff)
#y_ax_lbl <- expression(atop(Change~('%')~'from'~2010~runoff,
#                            ~from~2010))

input <- df_2_all_runs_basin %>% filter(year>=2010, year<=2050) %>%
  mutate(clim_imp_perc=if_else(year<2010, 0,clim_imp_perc)) %>%
  mutate(smoothedY=clim_imp_perc)
hist_temp <- df_2_all_runs_basin %>% filter(year<=2010, year>=1990) %>% mutate(clim_imp_perc=0)

roll <- 2
start_yr <- 2010
end_yr <- 2050
start_yr_hist <- 1990
end_yr_hist <- 2010
var_names <- c('Basin_runoff_km3peryear')
region_list <- water_basins_plot
gcm_list <- 'HadGEM2-ES' # 'GFDL-ESM2M' # IPSL-CM5A-LR'
rcp_list <- c('rcp8p5')  # 'rcp2p6'
region_single_plot(var_names, region_list, input, figures_basepath, start_yr, end_yr, gcm_names, rcp_names,
                   roll, y_ax_lbl, trendline=0, combined_lines=1, plot_df_hist=hist_temp,
                   all_same_color = 1, titles = 'No', legend_on=F, plot_hist=FALSE, plot_var='perc_red',
                   xmin=2010, xmax=2050, gcm_list=gcm_list, rcp_list=rcp_list, fig_type='.pdf')

#**************************************** PLOT END ******************************************#


#**************************************** CSV & XML START ******************************************#
# Having produced all plots, now save file as csv, in format that will allow it to be converted into gcam-ready xml
variable <- 'runoff'
level2outdir <- 'E:/NEXO-UA/Results/impacts/water_avail'
renewrsc_max_gcam <- paste0(extras_dir, '/', "L201.RenewRsrcCurves_calib_watergap.csv")
read_csv(renewrsc_max_gcam, skip = 4) %>%
  select(region, renewresource) %>% unique() -> region_basin
df_2_all_runs_basin <- df_2_all_runs_basin
name_exclude_list <- c('Andamaman-Nicobar Islands', 'North Marina Islands-Guam', 'Micronesia')
gcam_xanthos_basin_mapping <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output/gcam_xanthos_basin_mapping.csv'
write_csv_file(df_2_all_runs_basin, gcam_years, gcm_names, rcp_names, level2outdir, variable, region_basin=region_basin,
               gcam_xanthos_basin_mapping=gcam_xanthos_basin_mapping, name_exclude_list=name_exclude_list)


# Read in xanthos output file and region files
gcamdatafolder = "E:/NEXO-UA/GCAM-Workspace/gcamdata"

startYear=2015
basin_id <- rep(seq(from=1, to=235, by=1), each=150, times=20)
df_xml <- df_2_all_runs_basin %>% rename(maxSubResource = clim_imp_val) %>% mutate(basin_id=basin_id)

# Prepare the Data by GCAM Basin Region
gcam_prep <- function(df_xml_gcam, extra_dir){
  if(T){
    # From .input/gcamdata/R/zchunk_L201.water_resources_constrained
    GCAM_region_names <- data.table::fread(paste(extras_dir,"/GCAM_region_names.csv",sep=""),header=TRUE)
    iso_GCAM_regID <- data.table::fread(paste(extras_dir,"/iso_GCAM_regID.csv",sep=""),header=TRUE)
    basin_to_country_mapping <- data.table::fread(paste(extras_dir,"/basin_to_country_mapping.csv",sep=""),header=TRUE)
    basin_ids <- data.table::fread(paste(extras_dir,"/basin_ID.csv",sep=""),header=TRUE)
    water_mapping_R_GLU_B_W_Ws_share <- data.table::fread(paste(extras_dir,"/L103.water_mapping_R_GLU_B_W_Ws_share.csv",sep=""),header=TRUE)
    water_mapping_R_B_W_Ws_share <- data.table::fread(paste(extras_dir,"/L103.water_mapping_R_B_W_Ws_share.csv",sep=""),header=TRUE)
    
    # Basin_to_country_mapping table include only one set of distinct basins
    # that are mapped to a single country with largest basin share.
    # Assign GCAM region name to each basin.
    # Basin with overlapping GCAM regions assign to region with largest basin area.
    basin_to_country_mapping %>%
      rename(iso = ISO) %>%
      mutate(iso = tolower(iso)) %>%
      left_join(iso_GCAM_regID, by = "iso") %>%
      # ^^ non-restrictive join required (NA values generated for unmapped iso)
      # basins without gcam region mapping excluded (right join)
      # Antarctica not assigned
      right_join(GCAM_region_names, by = "GCAM_region_ID") %>%
      rename(basin_id = GCAM_basin_ID,
             basin_name = Basin_name) %>%
      select(GCAM_region_ID, region, basin_id) %>%
      arrange(region) ->
      RegionBasinHome
    
    # identify basins without gcam region mapping (anti_join)
    basin_to_country_mapping %>%
      rename(iso = ISO) %>%
      mutate(iso = tolower(iso)) %>%
      left_join(iso_GCAM_regID, by = "iso") %>%
      #not all iso included in basin mapping
      # ^^ non-restrictive join required (NA values generated for unmapped iso)
      anti_join(GCAM_region_names, by = "GCAM_region_ID") ->
      BasinNoRegion
    
    # create full set of region/basin combinations
    # some basins overlap multiple regions
    # Use left join to ensure only those basins in use by GCAM regions are included
    bind_rows(water_mapping_R_GLU_B_W_Ws_share %>%
                rename(basin_id = GLU),
              water_mapping_R_B_W_Ws_share) %>%
      select(GCAM_region_ID, basin_id, water_type) %>%
      filter(water_type == "water withdrawals") %>%
      distinct() %>%
      left_join(basin_ids, by = "basin_id") %>%
      tibble::as_tibble()%>%
      # ^^ non-restrictive join required (NA values generated for unused basins)
      left_join_error_no_match(GCAM_region_names, by = "GCAM_region_ID") %>%
      mutate(water_type = "water withdrawals",
             resource = paste(basin_name, water_type, sep="_")) %>%
      arrange(region, basin_name) ->
      L201.region_basin
    
    # create unique set of region/basin combination with
    # basin contained by home region (region with largest basin area)
    L201.region_basin %>%
      inner_join(RegionBasinHome, by = c("basin_id","GCAM_region_ID","region")) %>%
      arrange(region, basin_name) ->
      L201.region_basin_home
    
    # Re-format to format for ./input/gcamdata/outputs/L201.GrdRenwRsrcMax_runoff.csv which include
    # region, renewresource, sub.renewable.resource, year, maxSubResource
    
    df <- L201.region_basin_home %>%
      dplyr::left_join(df_xml_gcam %>%
                         dplyr::select(basin_id, year, maxSubResource) %>%
                         dplyr::filter(year %in% c(seq(startYear,2100,by=5)))
                       ,by="basin_id") %>%
      dplyr::mutate(sub.renewable.resource="runoff",
                    basin_name=paste(basin_name,"_water withdrawals",sep=""))%>%
      dplyr::select(id=basin_id,region=region,renewresource=basin_name, sub.renewable.resource, year,maxSubResource);
    
    df %>% filter(is.na(renewresource), year==1975)
    df %>% filter(year==1975) %>% nrow()
    df %>% filter(region=="USA",year==1975)%>%arrange(renewresource)
    df %>% filter(region=="China",year==1975)%>%arrange(renewresource)
    df %>% filter(region=="Southeast Asia",year==1975)%>%arrange(renewresource)
    df %>% filter(renewresource=="Hong-Red River")%>%arrange(renewresource)
    df$year%>%unique()
  }
  
  gcam_file <- paste0('E:/NEXO-UA/Results/impacts/water_avail/xml/', paste('runoff_impacts', gcm_gcam, rcp_gcam, sep = '_'), '.xml')
  create_xml(gcam_file) %>%
    add_xml_data(df, "GrdRenewRsrcMaxNoFillOut")%>%
    run_xml_conversion()
}

# prepare GCAM data by GCM and RCP
for(gcm_gcam in gcm_names){
  for(rcp_gcam in rcp_names){
    df_xml_gcam <- df_xml %>% filter(gcm %in% gcm_gcam, rcp %in% rcp_gcam)
    gcam_prep(df_xml_gcam, extra_dir)
  }
}

#**************************************** CSV & XML END ******************************************#


#**************************************** No Longer In Use ******************************************#
# Save as xml
# Use header: GrdRenewRsrcMaxNoFillOut from .\input\gcamdata\inst\extdata\mi_headers\ModelInterface_headers
create_xml("E:/NEXO-UA/Results/impacts/water_avail/xml/runoff_impacts_HadGEM2-ES_rcp6p0.xml") %>%
  add_xml_data(df, "GrdRenewRsrcMaxNoFillOut")%>%
  run_xml_conversion()

#
#
# # Convert csv files to gcam-ready xml files
# var <- "GrdRenewRsrcMax" #  "AgProdChange"
# csvpath <- '/Users/mengqi/Documents/NEXO-UA/Results/impacts/water_avail' # "data/scenario_agprodchange_gcam513_annual"
# xmlpath <- '/Users/mengqi/Documents/NEXO-UA/Results/impacts/water_avail/xml'  # "data/scenario_agprodchange_gcam513_annual/xml"
# csv2xml(csvpath, xmlpath, var)

# Plots to confirm delta scaling and smoothing worked correctly
df_plot <- df_2_all_runs_basin %>% filter(gcm=='HadGEM2-ES', rcp=='rcp6p0', name=='Magdalena')
p <- ggplot(data=df_plot) + geom_line(mapping = aes(x = year, y = value, colour=gcm))
p <- p + geom_line(color='blue', mapping = aes(x = year, y = smoothedY))
p <- p + geom_line(color='green', mapping = aes(x = year, y = clim_imp_val))
p
