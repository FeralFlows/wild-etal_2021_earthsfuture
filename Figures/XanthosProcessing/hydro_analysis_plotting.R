# Purpose: Plots Climate Impacts on Hydropower from Xanthos

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

rcp_colors <- c("rcp8p5" = "#736F6E",
                "rcp6p0" = "#C0C0C0",
                "rcp4p5" = "#98AFC7",
                "rcp2p6" = "#6698FF",
                "historical" = 'black')
gcm_colors <- c("NorESM1-M" = "#736F6E",
                "MIROC-ESM-CHEM" = "#C0C0C0",
                "IPSL-CM5A-LR" = "#98AFC7",
                "HadGEM2-ES" = "#6698FF",
                "GFDL-ESM2M" = "#153E7E",
                "watch+wfdei" = 'black')

figures_basepath <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output/figures/hydro'
results_basepath <- 'E:/NEXO-UA/Xanthos/example/output'
csv_basepath <- 'E:/NEXO-UA/Results/impacts/hydro'
xanthos_config_names <- c('clim_impacts')
gcm_names <- c('NorESM1-M', 'MIROC-ESM-CHEM', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'GFDL-ESM2M')
gcm_names_incl_hist <- append(gcm_names, 'watch+wfdei')
rcp_names <- c('rcp2p6', 'rcp4p5', 'rcp6p0', 'rcp8p5')
rcp_names_incl_hist <- append(rcp_names, "historical")
xanthos_var_names <- c('actual_hydro_by_gcam_region_EJperyr')
time_scale <- '1950_2099'
country_list_plot <- c('Colombia', 'Argentina', 'Uruguay')
run_name <- c('clim_impacts')
gcam_years <- c(2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065, 2070, 2075,
                2080, 2085, 2090, 2095, 2100)
stored_in_dir <- 1  # = 1 if dragged whole xanthos folder off pic; 0 if just dragged file down into results dir on comp
# Add historical values to dataframe
add_historical <- 1
delta_correction <- 0

# Process GCAM hydro values
# Deal first with all data
gcam_ref_hydro_ALL_basepath <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output/reference_gcam_hydro_all_regions.csv'
gcam_ref_hydro_ALL_traject <- read.csv(gcam_ref_hydro_ALL_basepath) %>%
  gather("year", "value", -scenario, -region, -subsector, -Units) %>% mutate(year=str_replace(year, "X", '')) %>%
  select(-scenario, -Units, -subsector) %>% rename(reference=value, period=year)
gcam_ref_hydro_ALL_traject$period <- as.numeric(gcam_ref_hydro_ALL_traject$period)
country_list_full <- as.character(unique(gcam_ref_hydro_ALL_traject$region))
gcam_ref_hydro_traject_final <- gcam_ref_hydro_ALL_traject[FALSE,]
# Loop through years, and add year row if it doesnt already exist
year_vector <- 1990:2100
gcam_year_list <- gcam_ref_hydro_ALL_traject$period
for (reg in country_list_full){
  append_row <- data.frame('region'=reg, 'period'=9999, 'reference' = NA, 'delta'=1)
  temp <- gcam_ref_hydro_ALL_traject %>% filter(region==reg) %>% mutate('delta'=1)
  for(t in year_vector){
    if(!t %in% unique(gcam_year_list)){
      append_row$period <- t
      temp <- rbind(temp, append_row)
    }
  }
  temp <- temp %>% mutate(reference=na.approx(reference, period))  # interpolate reference between gcam time steps
  reference_2010 <- (gcam_ref_hydro_ALL_traject %>% filter(region==reg, period==2010))$reference[1]
  temp <- temp %>% mutate(delta = reference/reference_2010)  # calculate delta values for all years. Delta factor is for energy expansion (via capacity expansion)
  gcam_ref_hydro_traject_final <- rbind(gcam_ref_hydro_traject_final, temp)
}
# final processing steps
gcam_ref_hydro_traject_final <- gcam_ref_hydro_traject_final %>% rename(name=region, year=period)
filter_list_2 <- list("actual_hydro_by_gcam_region_EJperyr" = country_list_full)

# Process Xanthos run data
df_all_runs_hydro <- xanthos_proc(xanthos_var_names, xanthos_config_names, gcm_names, rcp_names, time_scale, results_basepath,
                                  filter_list = filter_list_2)$output
df_all_runs_hydro <- xanthos_hist_proc(xanthos_var_names, xanthos_config_names, df_all_runs_hydro, stored_in_dir,
                                       results_basepath, add_historical, filter_list=filter_list_2)$output
df_all_runs_hydro_hist <- df_all_runs_hydro %>% filter(rcp == 'historical')
df_all_runs_hydro <- df_all_runs_hydro %>% filter(rcp != 'historical')

df_all_runs_hydro$year <- as.numeric(df_all_runs_hydro$year)
df_all_runs_hydro_hist$year <- as.numeric(df_all_runs_hydro_hist$year)

# Compute rolling mean
# Projected
df_2_all_runs_hydro <- roll_mean(df_all_runs_hydro, xanthos_var_names, xanthos_config_names, gcm_names_incl_hist,
                                 rcp_names_incl_hist, region_list=filter_list_2, loess_span=1)$output
df_2_all_runs_hydro$year <- as.numeric(df_2_all_runs_hydro$year)
# Historical
df_2_all_runs_hydro_hist <- roll_mean(df_all_runs_hydro_hist, xanthos_var_names, xanthos_config_names, gcm_names_incl_hist,
                                 rcp_names_incl_hist, region_list=filter_list_2, loess_span=1)$output
df_2_all_runs_hydro_hist$year <- as.numeric(df_2_all_runs_hydro_hist$year)

# Apply delta factor for installed capacity to the smoothed future hydro gen values that are only affected by runoff.
df_2_all_runs_hydro <- df_2_all_runs_hydro %>% left_join(gcam_ref_hydro_traject_final, by=c('name', 'year'))
df_2_all_runs_hydro$delta[is.na(df_2_all_runs_hydro$delta)] <- 1  # Set NA values in delta column to 1

# Compute the mean value, store it for every gcm/rcp combo. Also compute the percent change in every year relative to
# that mean, and store that
df_2_all_runs_hydro['mean_2010'] <- 0 # add mean_2010 column
for (reg in country_list_full){
  for (gcm1 in gcm_names){
    for (rcp1 in rcp_names){
      mean_val <- (df_2_all_runs_hydro %>% filter(name==reg, gcm==gcm1, rcp==rcp1, year==2010))$smoothedY[1]
      df_2_all_runs_hydro <- df_2_all_runs_hydro %>% mutate(mean_2010 = if_else(name==reg & gcm==gcm1 & rcp==rcp1, mean_val, mean_2010))
    }
  }
}
df_2_all_runs_hydro <- df_2_all_runs_hydro %>% mutate(clim_imp_perc = 100*(smoothedY-mean_2010)/mean_2010) %>% select(-mean_2010)

#
df_2_all_runs_hydro <- df_2_all_runs_hydro %>% mutate(clim_imp_val = reference*(1+clim_imp_perc/100)) %>%
  mutate(clim_imp_val=if_else(clim_imp_val<0,0,clim_imp_val))

# adjust smoothedY and rolling_mean by delta
df_2_all_runs_hydro <- df_2_all_runs_hydro %>% mutate(RollingMeanDelta=rolling_mean*delta) %>%
  mutate(smoothedYDelta=smoothedY*delta)
# Adjust smoothedY so it does or does not reflect a delta correction, depending on user preferences
if(delta_correction==1){
  df_2_all_runs_hydro$smoothedY <- df_2_all_runs_hydro$smoothedYDelta
}


# Plot country hydropower where all the GCM and RCP combinations are
# combined on the same plot
y_ax_lbl <- expression(Annual~Hydropower~(TWh))
input <- df_2_all_runs_hydro %>% filter(year>=2010, year<=2050) %>% mutate(smoothedY=clim_imp_val)
roll <- 2
start_yr <- 2010
end_yr <- 2050
start_yr_hist <- 1970
end_yr_hist <- 2010
var_names <- c('actual_hydro_by_gcam_region_EJperyr')
region_list <- c('Argentina')  # country_list_plot
gcm_list <- 'MIROC-ESM-CHEM' #'HadGEM2-ES'  # 'GFDL-ESM2M', 'IPSL-CM5A-LR'
rcp_list <- c('rcp6p0')  # 'rcp2p6' #'rcp8p5'
region_single_plot(var_names, region_list, input, figures_basepath, start_yr, end_yr, gcm_names, rcp_names,
                   roll, y_ax_lbl, trendline=0, combined_lines=1, plot_df_hist=df_2_all_runs_hydro_hist,
                   all_same_color = 1, titles = 'Yes', legend_on=TRUE, xmin=2010, xmax=2050, plot_reference=TRUE,
                   gcm_list=gcm_list, rcp_list=rcp_list, fig_type='.pdf')

# Plot percentage reduction in smoothed hydropower production compared with 2010
y_ax_lbl <- expression(atop(Change~('%')~'in'~hydropower,
                       ~generation~from~2010))
input <- df_2_all_runs_hydro %>% filter(year>=2010, year<=2050)
input <- input %>% mutate(smoothedY=clim_imp_perc)
roll <- 2
start_yr <- 2010
end_yr <- 2050
start_yr_hist <- 1970
end_yr_hist <- 2010
var_names <- c('actual_hydro_by_gcam_region_EJperyr')
region_list <- c('Argentina')  # country_list_plot
gcm_list <- 'MIROC-ESM-CHEM' #'HadGEM2-ES'  # 'GFDL-ESM2M', 'IPSL-CM5A-LR'
rcp_list <- c('rcp6p0')  # 'rcp2p6' #'rcp8p5'
region_single_plot(var_names, region_list, input, figures_basepath, start_yr, end_yr, gcm_names, rcp_names,
                   roll, y_ax_lbl, trendline=0, combined_lines=1, plot_df_hist=df_2_all_runs_hydro_hist,
                   all_same_color = 1, titles = 'Yes', legend_on=TRUE, plot_hist=FALSE, plot_var='perc_red',
                   xmin=2010, xmax=2050, gcm_list=gcm_list, rcp_list=rcp_list, fig_type='.pdf')

# Having produced all plots, now save file as csv, in format that will allow it to be converted into gcam-ready xml
variable <- 'hydro'
write_csv_file(df_2_all_runs_hydro, gcam_years, gcm_names, rcp_names, csv_basepath, variable)

# Convert csv files to gcam-ready xml files
var <- "StubTechFixOut" #  "AgProdChange"
csvpath <- 'E:/NEXO-UA/Results/impacts/hydro' # "data/scenario_agprodchange_gcam513_annual"
xmlpath <- 'E:/NEXO-UA/Results/impacts/hydro/xml'  # "data/scenario_agprodchange_gcam513_annual/xml"
csv2xml(csvpath, xmlpath, var)





#-----------------------------------------------------------------------------------------------------------------------

# PLOTS WE ARE NOT USING ANY LONGER

# Create faceted plot across GCMs and RCPs for country hydropower production
roll <- 0
start_yr <- 2010 # 2010
end_yr <- 2050 # 2100
start_yr_hist <- 1970
end_yr_hist <- 2009
xanthos_var_names <- c('actual_hydro_by_gcam_region_EJperyr')
filter_list_2 <- list("actual_hydro_by_gcam_region_EJperyr" = country_list_plot)
y_ax_lbl <- expression(Annual~Hydropower~(TWh))
for(var_1 in xanthos_var_names){
  for(reg in filter_list_2[[var_1]]){
    fig_name <- paste0(figures_basepath, '/', reg, "_", var_1, "_", 'gcm_rcp_facet.png')
    plot_df <- df_2_all_runs_hydro %>% filter(name==reg, year>=start_yr, year<=end_yr, var==var_1) %>%
      filter(gcm %in% gcm_names_incl_hist, rcp %in% rcp_names_incl_hist)
    plot_df_hist <- df_all_runs_hydro_hist %>% filter(name==reg, year>=start_yr_hist, year<=end_yr_hist, var==var_1)
    facet_grid_plot(plot_df, fig_name, rolling=roll, y_lbl=y_ax_lbl, df_all_runs_hist=plot_df_hist, historical=1)
  }
}

# Create smooth faceted plot across GCMs and RCPs for country hydropower production
roll <- 2
# Use all other values (except "roll" from above)
for(var_1 in xanthos_var_names){
  for(reg in filter_list_2[[var_1]]){
    fig_name <- paste0(figures_basepath, '/', reg, "_", var_1, "_", 'gcm_rcp_facet_smooth.png')
    plot_df <- df_2_all_runs_hydro %>% filter(name==reg, year>=start_yr, year<=end_yr, var==var_1) %>%
      filter(gcm %in% gcm_names_incl_hist, rcp %in% rcp_names_incl_hist)
    plot_df_hist <- df_all_runs_hydro_hist %>% filter(name==reg, year>=start_yr_hist, year<=end_yr_hist, var==var_1)
    facet_grid_plot(plot_df, fig_name, rolling=roll, y_lbl=y_ax_lbl, df_all_runs_hist=plot_df_hist, historical=1)
  }
}

# Hydropower: Countries
y_ax_lbl <- expression(Annual~Hydropower~(TWh))
input <- df_2_all_runs_hydro
roll <- 0
start_yr <- 2010
end_yr <- 2050
var_names <- c('actual_hydro_by_gcam_region_EJperyr')
region_list <- country_list_plot
region_single_plot(var_names, region_list, input, figures_basepath, start_yr, end_yr,
                   gcm_names, rcp_names, roll, y_ax_lbl, trendline=0)

# Smoothed Hydropower by Country
y_ax_lbl <- expression(Annual~Hydropower~(TWh))
input <- df_2_all_runs_hydro
roll <- 2
start_yr <- 2010
end_yr <- 2050
var_names <- c('actual_hydro_by_gcam_region_EJperyr')
region_list <- country_list_plot
region_single_plot(var_names, region_list, input, figures_basepath, start_yr, end_yr,
                   gcm_names, rcp_names, roll, y_ax_lbl, trendline=0)

