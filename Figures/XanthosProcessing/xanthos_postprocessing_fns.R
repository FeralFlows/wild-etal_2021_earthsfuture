library(tidyverse)
if("devtools" %in% rownames(installed.packages()) == F){install.packages("devtools")}
library(devtools)
if("gcamdata" %in% rownames(installed.packages()) == F){install_github(repo="JGCRI/gcamdata")}
library(dplyr)
library(tidyr)
library(foreach)
library(gcamdata)

#-----------------------------------------------------
# FUNCTIONS

line_plot <- function(plot_df, fig_name, rolling=0, y_lbl=NULL, x_lbl=NULL, y_max=NULL, y_min=NULL, trendline=1,
                      title=TRUE, legend_on=TRUE, x_min=NULL, x_max=NULL){

  # ggplot2 Theme
  z_theme <<- theme_bw() +
    theme(
      text =                element_text(family = NULL, face = "plain",colour = "black", size = 8 ,hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9)
      , axis.text.x =       element_text(size=10)
      , axis.text.y =       element_text(size=10)
      ,axis.title.x =       element_text(vjust = -1, margin=margin(t=1,unit="line"))
      ,axis.title.y =       element_text(angle = 90, vjust = 2, margin=margin(r=1,unit="line"))
      ,legend.key =         element_blank()
      ,legend.key.size =    unit(1.5, 'lines')
      ,legend.text =        element_text(size = 10, colour = "black")
      ,legend.title =       element_text(size = rel(1.2), face = NULL, hjust = 0, colour = "black")
      ,strip.background =   element_rect(fill = NA, colour = "black")
      ,plot.margin =        unit(c(1, 1, 1, 1), "lines")
      ,plot.title=          element_text(face="bold", hjust=0,size=12,margin = margin(b=20))
    )

  breakx_minMaster<-5
  breakx_majMaster<-10
  prettyBreaksyMaster<-5
  '%ni%' <- Negate('%in%')
  plot_df_orig <- plot_df %>% filter(rcp %ni% c('historical', 'historical mean'))  # Eliminate out historical if it exists
  plot_df_hist <- plot_df %>% filter(rcp %in% c('historical', 'historical mean'))  # Store historical values in separate DF to be plotted
  line_colors<-get(plot_df$FillPalette)
  #  p <- ggplot(plot_df_orig, aes(x=year, y=value, colour=gcm))
  p <- ggplot(data=plot_df_orig)

  if(rolling==1){
    p <- ggplot(data=plot_df, mapping = aes(x = year, y = rolling_mean, colour=gcm))
  }else if(rolling==2){
    p <- ggplot(data=plot_df, mapping = aes(x = year, y = smoothedY, colour=gcm))
  }else{
    p <- ggplot(data=plot_df, mapping = aes(x = year, y = value, colour=gcm))
  }
  p <- p + geom_point(size=0) + geom_line(size=0.5)  #, se=FALSE, , linetype=gcm
  if(trendline==1){
    p <- p + geom_smooth(method='lm', colour='blue', size=0.5, se=FALSE)
  }
  # Add historical points onto every plot
  p <- p + xlab(x_lbl) + ylab(y_lbl)

  if(!is.null(x_min)){
    p <- p + scale_x_continuous(limits=c(x_min, x_max))
  }
  if(!is.null(y_min)){
    p<-p + scale_y_continuous(limits=c(y_min - 0.1*y_min,1.1*y_max))
  }

  p<-p + scale_color_manual(values=line_colors)
  p <- p + ggtitle(title)
  if(legend_on==FALSE){
    p <- p + guides(color=legend_on)
  }
  p
  ggsave(fig_name, dpi=900, width=4, height=2.5, units="in")
}

line_plot_hist_proj <- function(plot_df, plot_df_hist, fig_name, gcm_names, rcp_names, rolling=0, y_lbl=NULL,
                                x_lbl=NULL, y_max=NULL, y_min=NULL, trendline=1, all_same_color=1, title=NULL, legend_on=TRUE,
                                plot_var=NULL, plot_hist=TRUE, x_min=NULL, x_max=NULL, plot_reference=NULL,
                                gcm_list=NULL, rcp_list=NULL){

  line_colors<-get(plot_df$FillPalette)

  # ggplot2 Theme
  z_theme <<- theme_bw() +
    theme(
      text =                element_text(family = NULL, face = "plain",colour = "black", size = 10 ,hjust = 0.5,
                                         vjust = 0.5, angle = 0, lineheight = 0.9)
      , axis.text.x =       element_text(size=8)
      , axis.text.y =       element_text(size=8)
      ,axis.title.x =       element_text(vjust = -1, margin=margin(t=1,unit="line"))
      ,axis.title.y =       element_text(angle = 90, vjust = 2, margin=margin(r=1,unit="line"))
      ,legend.key =         element_blank()
      ,legend.key.size =    unit(1.5, 'lines')
      ,legend.text =        element_text(size = 8, colour = "black")
      ,legend.title =       element_text(size = rel(1.2), face = NULL, hjust = 0, colour = "black")
      ,strip.background =   element_rect(fill = NA, colour = "black")
      ,plot.margin =        unit(c(1, 1, 1, 1), "lines")
      ,plot.title=          element_text(face="bold", hjust=0.2, vjust = -4, margin = margin(b=20), size=8)
      #,plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

  breakx_minMaster<-5
  breakx_majMaster<-10
  prettyBreaksyMaster<-5
  '%ni%' <- Negate('%in%')
  plot_df_orig <- plot_df %>% filter(rcp %ni% c('historical', 'historical mean'))  # Eliminate out historical if it exists

  if(!is.null(plot_df_hist)){
    line_colors_hist<-get(plot_df_hist$FillPalette)
    plot_df_hist <- plot_df_hist %>% filter(rcp %in% c('historical', 'historical mean'))  # Store historical values in separate DF to be plotted
  }
  #line_colors<-get(plot_df$FillPalette)

  # First, add historical data if user wants to plot it
  if (plot_hist==TRUE){
    if(rolling==1){
      p <- ggplot(data=plot_df_hist, mapping = aes(x = year, y = rolling_mean, colour=gcm, fill=gcm))
    }else if(rolling==2){
      p <- ggplot(data=plot_df_hist, mapping = aes(x = year, y = smoothedY))
    }else{
      p <- ggplot(data=plot_df_hist, mapping = aes(x = year, y = value, colour=gcm, fill=gcm))
    }
    #p <- p + geom_line(size=0.5, color='black', data=plot_df_hist, mapping = aes(x = year, y = value))
    p <- p + geom_line(size=0.5, color='black')
    p <- p + geom_line()
  }else{
    # Historical data not going to be plotted. Need to create new ggplot since it wasnt created above for hist plotting.
    p <- ggplot()
  }

  if(all_same_color==1){
    color_var = 'grey70' # "#153E7E"
  }else{
    color_var = NULL
  }

  for(gcm1 in gcm_names){
    for(rcp1 in rcp_names){
      filtered_df <- plot_df_orig %>% filter(rcp==rcp1, gcm==gcm1)
      if(rolling==1){
        if(all_same_color==1){
          p <- p + geom_line(size=0.5, color = color_var, data=filtered_df, mapping = aes(x = year, y = rolling_mean))
        }else{
          p <- p + geom_line(size=0.5, data=filtered_df, mapping = aes(x = year, y = rolling_mean,
                                                                                          colour=gcm))
        }
      }else if(rolling==2){
        if(all_same_color==1){
          p <- p + geom_line(size=0.5, color = color_var, data=filtered_df, mapping = aes(x = year, y = smoothedY))
        }else{
          p <- p + geom_line(size=0.5, data=filtered_df, mapping = aes(x = year, y = smoothedY,
                                                                                          colour=gcm))
        }

      }else{
        if(all_same_color==1){
          p <- p + geom_line(size=0.5, color = color_var, data=filtered_df, mapping = aes(x = year, y = value))
        }else{
          p <- p + geom_line(size=0.5, data=filtered_df, mapping = aes(x = year, y = value, colour=gcm))
        }
      }

    }
  }

  # Plot reference scenario
  if(!is.null(plot_reference)){
    p <- p + geom_line(size=0.5, linetype=1, color = 'black', data=filtered_df, mapping = aes(x = year, y = reference))
  }

  # SO far, only adding this colored gcm/rcp plots for smoothedY
  # Plot select subset of gcms. Maximum of two lines, or will generate error
  if(rolling==2){
    ctr <- 0
    color_list <- c('dodgerblue3', '#fc9272')  #  #de2d26
    if(!is.null(gcm_list)){
      for(model in gcm_list){
        for(forc in rcp_list){
          ctr <- ctr + 1
          if(ctr>2){
            Print("error: currently can only plot 2 lines of individual gcms/rcps as one color")
          }
          plot_df_orig_2 <- plot_df_orig %>% filter(gcm == model, rcp==forc)  # , rcp == c('rcp2p6')
          plot_df_orig_2$gcm <- as.character(plot_df_orig_2$gcm)
          plot_df_orig_2$rcp <- as.character(plot_df_orig_2$rcp)
          p <- p + geom_line(size=0.5, linetype=1, color = color_list[ctr],
                             data=plot_df_orig_2, mapping = aes(x = year, y = smoothedY))  # linetype=2
        }
      }
    }
  }

  p <- p + xlab(x_lbl) + ylab(y_lbl)
  if(!is.null(y_min)){
    p<-p + scale_y_continuous(limits=c(y_min - 0.1*abs(y_min), 1.1*y_max))
  }
  if(!is.null(x_min)){
    p<-p + scale_x_continuous(limits=c(x_min, x_max))
  }
  p<-p + scale_color_manual(values=line_colors, name = "Time Scale")
  if(!is.null(plot_df_hist)){
    p<-p + scale_color_manual(values=line_colors_hist)
  }
  if(legend_on==FALSE){
    p <- p + guides(color=legend_on)
  }
  p <- p + ggtitle(title)
  p <- p + z_theme
  p
  ggsave(fig_name, dpi=900, width=2.5, height=2.5, units="in")
}

# Faceted grid plot function
facet_grid_plot <- function(plot_df, fig_name, historical=0, rolling=0, y_lbl=NULL, x_lbl=NULL,
                            df_all_runs_hist=NULL){
  breakx_minMaster<-5
  breakx_majMaster<-10
  prettyBreaksyMaster<-5
  plot_df_orig <- plot_df
  line_colors<-get(plot_df$FillPalette)
  if(rolling==1){
    p <- ggplot(data=plot_df_orig, mapping = aes(x = year, y = rolling_mean, colour=gcm))
  }else if (rolling==2){
    p <- ggplot(data=plot_df_orig, mapping = aes(x = year, y = smoothedY, colour=gcm))
  }else{
      p <- ggplot(data=plot_df_orig, mapping = aes(x = year, y = value, colour=gcm))
  }
  # p <- p + geom_smooth(method=lm, colour='black', size=0.5, se=FALSE)  # Linear best fit line
  p <- p + geom_line(size=0.5) #  Used to be first: geom_point(size=0) +  #, se=FALSE, , linetype=gcm

  p <- p + xlab(x_lbl) + ylab(y_lbl)
  p <- p + scale_color_manual(values=line_colors)
  #p <- p + scale_color_brewer(palette='GnBu')
  #p <- p + scale_size_manual(values=c(100, 100, 100, 100))

  if(historical==1){
    # Recreate a data frame that copies the historical values over muiltiple times and lists them with each permutation
    # of GCMs and RCPs, for purposes of facet plotting
    rcp_list <- unique(plot_df$rcp)
    gcm_list <- unique(plot_df$gcm)
    hist_df <- data.frame()
    for(rcp_1 in rcp_list){
      for (gcm_1 in gcm_list){
        input <- df_all_runs_hist  # historical data frame
        input['rcp'] <- rcp_1  # write over the historical rcp category and replace them with each of the rcps
        input['gcm'] <- gcm_1  # write over the historical gcm category and replace them with each of the gcms
        hist_df <- rbind(hist_df, input)  # create new data frame storing historical data with all gcm/rcp combinations in columns
      }
    }
    p <- p + geom_line(size=0.5, color='black', data=hist_df, mapping = aes(x = year, y = value)) #+ geom_line(size=0.5)
  }

  p <- p + facet_grid(vars(gcm), vars(rcp)) + theme(panel.spacing.x = unit(4, "mm"), legend.position="bottom")
  p
  ggsave(fig_name, dpi=900, width=6.5, height=8, units="in")
}

# imports xanthos outputs and reorganizes the data and sets new column values (e.g., rcp and gcm)
xanthos_proc <- function(xanthos_var_names, xanthos_config_names, gcm_names, rcp_names, time_scale, results_basepath,
                         filter_list=NULL, country_grid_id_filepath=NULL, country_names_id=NULL){
  df_all_runs <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(df_all_runs) <- c('name', 'year', 'value', 'mod', 'var', 'gcm', 'rcp', 'FillPalette')
  for(var in xanthos_var_names){
    for(mod in xanthos_config_names){
      for(gcm in gcm_names){
        for(rcp in rcp_names){
          run_name <- paste0(mod, "_", gcm, "_", rcp)  # , "_", time_scale
          if (stored_in_dir==1){
            xanthos_dir <- paste0(results_basepath, '/', run_name)
          }else{
            xanthos_dir <- paste0(results_basepath)
          }
          xanthos_file <- paste0(var, "_", gcm, "_", rcp, "_", time_scale, '.csv')
          xanthos_output_filepath <- paste0(xanthos_dir, '/', xanthos_file)
          if(dir.exists(xanthos_dir)){
            if(var == 'Basin_runoff_km3peryear'){
              input <- read_csv(xanthos_output_filepath) %>% select(-id) %>% gather(year, value, `1950`:`2099`)
              if(!is.null(filter_list)){
                input <- input %>% filter(name %in% filter_list[[var]])
              }
            }else if(var == 'q_km3peryear'){
              grid_country_map <- read_csv(country_grid_id_filepath)
              country_names_id_tbl <- read_csv(country_names_id)
              input <- read_csv(xanthos_output_filepath) %>% gather(year, value, `1950`:`2099`) %>%
                left_join(grid_country_map, by="id") %>% left_join(country_names_id_tbl, by="ctry_code") %>%
                group_by(region, year) %>% summarize(value=sum(value)) %>% ungroup() %>% rename(name=region) %>%
                filter(name %in% filter_list[[var]])
            }else if(var=='actual_hydro_by_gcam_region_EJperyr'){
              rgn32Names <- read_csv('E:/NEXO-UA/Results/downscaling/water/xanthos/output/Rgn32Names.csv')
              input <- read_csv(xanthos_output_filepath) %>% rename(id=region) %>% left_join(rgn32Names, by='id') %>%
                select(-id) %>% gather(year, value, `1950`:`2099`) %>% rename(name=region) %>%
                filter(name %in% filter_list[[var]]) %>% mutate(value=value*277.78)
          }
            input['mod'] <- mod
            input['var'] <- var
            input['gcm'] <- gcm
            input['rcp'] <- rcp
            input['FillPalette'] <- c('gcm_colors')
            df_all_runs <- rbind(df_all_runs, input)
          }
        }
      }
    }
  }
  return(list("output" = df_all_runs))
}

xanthos_hist_proc <- function(xanthos_var_names, xanthos_config_names, df_all_runs, stored_in_dir, results_basepath,
                              add_historical, filter_list=NULL, country_grid_id_filepath=NULL, country_names_id=NULL){
  for(mod in xanthos_config_names){
      if(add_historical==1){
        for(var in xanthos_var_names){
          run_name <- paste0(mod, "_", 'watch+wfdei', "_", '1970_2010')
          if(stored_in_dir==1){
            xanthos_dir <- paste0(results_basepath, '/', run_name)
          }else{
            xanthos_dir <- paste0(results_basepath)
          }
          xanthos_file <- paste0(var, "_", 'watch+wfdei', "_", '1970_2010', '.csv')
          xanthos_output_filepath <- paste0(xanthos_dir, '/', xanthos_file)
          if(dir.exists(xanthos_dir)){
            if(var=='Basin_runoff_km3peryear'){
              input <- read_csv(xanthos_output_filepath) %>% select(-id) %>% gather(year, value, `1970`:`2010`)
              if(!is.null(filter_list)){
                input <- input %>% filter(name %in% filter_list[[var]])
              }
            }else if(var=='q_km3peryear'){
              grid_country_map <- read_csv(country_grid_id_filepath)
              country_names_id_tbl <- read_csv(country_names_id)
              input <- read_csv(xanthos_output_filepath) %>% gather(year, value, `1970`:`2010`) %>%
                left_join(grid_country_map, by="id") %>% left_join(country_names_id_tbl, by="ctry_code") %>%
                group_by(region, year) %>% summarize(value=sum(value)) %>% ungroup() %>% rename(name=region) %>%
                filter(name %in% filter_list[[var]])
            }else if (var=='actual_hydro_by_gcam_region_EJperyr'){
              rgn32Names <- read_csv('E:/NEXO-UA/Results/downscaling/water/xanthos/output/Rgn32Names.csv')
#              input <- read_csv(xanthos_output_filepath) %>% left_join(rgn32Names, by='X1')
#              select(-X1) %>% gather(year, value, `1950`:`2099`) %>% rename(name=region) %>%
#                filter(name %in% filter_list[[var]])
              input <- read_csv(xanthos_output_filepath) %>% rename(id=region) %>% left_join(rgn32Names, by='id') %>%
                select(-id) %>% gather(year, value, `1970`:`2010`) %>% rename(name=region) %>%
                filter(name %in% filter_list[[var]]) %>% mutate(value=value*277.78)

            }
            input['mod'] <- xanthos_config_names
            input['var'] <- var
            input['gcm'] <- 'watch+wfdei'
            input['rcp'] <- 'historical'
            input['FillPalette'] <- c('gcm_colors')
            df_all_runs <- rbind(df_all_runs, input)
          }
        }
      }
    }
    return(list("output" = df_all_runs))
}

# imports xanthos outputs and reorganizes the data and sets new column values (e.g., rcp and gcm)
agmip_proc <- function(agmip_var_names, agmip_config_names, gcm_names, rcp_names, results_basepath,
                         filter_list, water_basin_abbrevc = NULL, country_grid_id_filepath=NULL, country_names_id=NULL,
                       filter_list_2=NULL, gcm_conversion=NULL){
  df_all_runs <- data.frame(matrix(ncol = 14, nrow = 0))
  colnames(df_all_runs) <- c('region', 'basin', 'year', 'value', 'mod', 'var', 'gcm', 'rcp', 'ssp', 'FillPalette')
  for(var in agmip_var_names){
    for(mod in agmip_config_names){
      for(gcm in gcm_names){
        for(rcp in rcp_names){
          for(ssp in ssp_names){
            run_name <- paste0(rcp, "_", gcm, "_", mod, "_", ssp)
            if (stored_in_dir==1){
              agmip_dir <- paste0(results_basepath, '/', run_name)
            }else{
              agmip_dir <- paste0(results_basepath)
            }
            agmip_file <- paste0(var, "_", run_name, '.csv')
            agmip_output_filepath <- paste0(agmip_dir, '/', agmip_file)
            if(dir.exists(agmip_dir)){
              if(var == 'ag_prodchange'){
                input <- read_csv(agmip_output_filepath, skip=4) %>% filter(region %in% filter_list[[var]]) %>%
                  rename(crop = AgSupplySector) %>% rename(basin=AgSupplySubsector) %>% rename(value=AgProdChange) # %>% mutate_if(is.character, str_replace_all, pattern = '', replacement = '')
                #input$basin <- gsub(input$crop, "", input$basin)  # remove crop from basin category
                input$basin <- sapply(seq_along(input$crop), function(x) gsub(input$crop[x], rep("", nrow(input)), input$basin[x]))  # remove crop from basin category
                input <- input %>% mutate(grass = "grass")
                input <- input %>% mutate(tree = "tree")
                input$basin <- sapply(seq_along(input$grass), function(x) gsub(input$grass[x], rep("", nrow(input)), input$basin[x]))  # remove grass
                input$basin <- sapply(seq_along(input$tree), function(x) gsub(input$tree[x], rep("", nrow(input)), input$basin[x]))  # remove tree
                input$basin <- gsub('_', "", input$basin)  # remove underscore
                # input$basin <- mgsub(input$crop, "", input$basin)  # remove crop from basin category
                # input$AgProductionTechnology <- gsub(paste0(input$crop,"_",input$basin,'_'), "", input$AgProductionTechnology)  # remove crop and basin from AgProdTech column
                input$rfd <- grepl("RFD", input$AgProductionTechnology)
                input$irr <- grepl("IRR", input$AgProductionTechnology)
                input$hi <- grepl("hi", input$AgProductionTechnology)
                input$lo <- grepl("lo", input$AgProductionTechnology)
                input <- input %>% mutate(value = 100*value)
                input_temp1 <- input %>% filter(irr == T, hi == T)
                input_temp2 <- input %>% filter(rfd == T, hi == T)
                input <- rbind(input_temp1, input_temp2) %>% select(-hi, -lo, -grass, -tree)
                for (bas in filter_list_2[[var]]){
                  input$basin <- gsub(bas, water_basin_abbrevc[bas], input$basin)  # remove crop from basin category
                }
              }
              input['mod'] <- mod
              input['var'] <- var
              input['gcm'] <- gcm_conversion[gcm]
              input['rcp'] <- rcp
              input['ssp'] <- ssp
              input['FillPalette'] <- c('gcm_colors')
              df_all_runs <- rbind(df_all_runs, input)
            }
          }
        }
      }
    }
  }
  return(list("output" = df_all_runs))
}


yield_proc <- function(agmip_var_names, yield_2010, filter_list, water_basin_abbrevc = NULL,
                       country_names_id=NULL, filter_list_2=NULL){

  input <- read_csv(yield_2010) %>% filter(region %in% filter_list) %>%
    rename(crop = AgSupplySector) %>% rename(basin=AgSupplySubsector)

  if(agmip_var_names=='crop_yield'){
    # Just pull in historical values, which don't require a year
    input <- input %>% select(-year)
  }
  input$basin <- sapply(seq_along(input$crop), function(x) gsub(input$crop[x], rep("", nrow(input)), input$basin[x]))  # remove crop from basin category
  input <- input %>% mutate(grass = "grass")
  input <- input %>% mutate(tree = "tree")
  input$basin <- sapply(seq_along(input$grass), function(x) gsub(input$grass[x], rep("", nrow(input)), input$basin[x]))  # remove grass
  input$basin <- sapply(seq_along(input$tree), function(x) gsub(input$tree[x], rep("", nrow(input)), input$basin[x]))  # remove tree
  input$basin <- gsub('_', "", input$basin)  # remove underscore
  input$rfd <- grepl("RFD", input$AgProductionTechnology)
  input$irr <- grepl("IRR", input$AgProductionTechnology)
  input$hi <- grepl("hi", input$AgProductionTechnology)
  input$lo <- grepl("lo", input$AgProductionTechnology)
  input_temp1 <- input %>% filter(irr == T, hi == T)
  input_temp2 <- input %>% filter(rfd == T, hi == T)
  input <- rbind(input_temp1, input_temp2) %>% select(-hi, -lo, -grass, -tree)
  for (bas in filter_list_2){
    input$basin <- gsub(bas, water_basin_abbrevc[bas], input$basin)  # remove crop from basin category
  }
  return(list("output" = input))
}



# Compute rolling mean
roll_mean <- function(df_all_runs, xanthos_var_names, xanthos_config_names, gcm_names_incl_hist, rcp_names_incl_hist,
                              region_list=NULL, k=10, loess_span=0.5){
  if(is.null(region_list)){
    region_list <- list()
    for(var1 in xanthos_var_names){
      region_list[[var1]] <- unique(df_all_runs$name)
    }
  }
  df_2_all_runs <- data.frame()
  for(var1 in xanthos_var_names){
    for(mod1 in xanthos_config_names){
      for(gcm1 in gcm_names_incl_hist){
        for(rcp1 in rcp_names_incl_hist){
          for(reg in region_list[[var1]]){
            temp <- df_all_runs %>% filter(gcm==gcm1, rcp==rcp1, var==var1, mod==mod1, name==reg) %>%
              mutate(rolling_mean=rollmean(value, k, fill=NA))
            if(nrow(temp)>0){
              loess_func <- loess(value~year, temp, span=loess_span)
              temp <- temp %>% mutate(smoothedY=predict(loess_func))
            }
            df_2_all_runs <- rbind(df_2_all_runs, temp)
          }
        }
      }
    }
  }
  df_2_all_runs$rolling_mean[is.na(df_2_all_runs$rolling_mean)] <- df_2_all_runs$value[is.na(df_2_all_runs$rolling_mean)]
  return(list("output" = df_2_all_runs))
}

region_single_plot <- function(xanthos_var_names, region_list, df_all_runs, figures_basepath, start_yr, end_yr,
                               gcm_names, rcp_names, roll, y_ax_lbl, trendline=1, combined_lines=0, plot_df_hist=NULL,
                               all_same_color = 1, titles=NULL, legend_on=TRUE, plot_var='', plot_hist=TRUE, xmin=NULL,
                               xmax=NULL, plot_reference=NULL, fig_name_append=NULL, gcm_list=NULL, rcp_list=NULL,
                               ymax=NULL, ymin=NULL, fig_type='.png'){
  for(var_1 in xanthos_var_names){
    for(reg in region_list){
      if(is.null(ymax) & is.null(ymin)){
        if(roll==1){
          ymax_across_gcms <- max((df_all_runs %>% filter(name==reg, var==var_1))$rolling_mean)
          ymin_across_gcms <- min((df_all_runs %>% filter(name==reg, var==var_1))$rolling_mean)
        }else if(roll==2){
          ymax_across_gcms <- max((df_all_runs %>% filter(name==reg, var==var_1))$smoothedY)
          ymin_across_gcms <- min((df_all_runs %>% filter(name==reg, var==var_1))$smoothedY)
        }else{
          ymax_across_gcms <- max((df_all_runs %>% filter(name==reg, var==var_1))$value)
          ymin_across_gcms <- min((df_all_runs %>% filter(name==reg, var==var_1))$value)
        }
      }else{
        # User has specified ymin and ymax
        ymax_across_gcms <- ymax
        ymin_across_gcms <- ymin
      }
      if(is.null(titles)){
        title=NULL
      }else if (titles %in% c('yes', 'Yes', 'y', 'Y')){
        title=reg
      }else{
        title=NULL
      }
      if(combined_lines == 0){
        for(gcm1 in gcm_names){
          fig_name <- paste0(figures_basepath, '/', var_1, "_", reg, "_", gcm1, "_", plot_var, "_", if(roll==1){'rolling_mean'}else if(roll==2){'loess'}else{''}, fig_type)
          plot_df <- df_all_runs %>%
            filter(name==reg, gcm==gcm1, year>=start_yr, year<=end_yr, gcm %in% gcm_names, rcp %in% rcp_names, var==var_1)
          if(nrow(plot_df)>0){
            line_plot(plot_df, fig_name, rolling=roll, y_lbl=y_ax_lbl, y_max=ymax_across_gcms, y_min=ymin_across_gcms,
                      trendline=trendline, title=title, legend_on=legend_on, x_min=xmin, x_max=xmax)
          }
        }
      }else{
        fig_name <- paste0(figures_basepath, '/', var_1, "_", reg, "_", plot_var, "_", fig_name_append, "_combined", if(roll==1){'_rolling_mean'}else if(roll==2){'_loess'}else{''}, fig_type)
        plot_df <- df_all_runs %>%
          filter(name==reg, year>=start_yr, year<=end_yr, gcm %in% gcm_names, rcp %in% rcp_names, var==var_1)
        if(plot_hist==TRUE){
          plot_df_hist_2 <- plot_df_hist %>% filter(name == reg, year<=2010)
        }else{
          plot_df_hist_2 <- NULL
        }
        print(reg)
        if(nrow(plot_df)>0){
          line_plot_hist_proj(plot_df, plot_df_hist_2, fig_name, gcm_names, rcp_names, rolling=roll, y_lbl=y_ax_lbl,
                              y_max=ymax_across_gcms, y_min=ymin_across_gcms, x_min=xmin, x_max=xmax,
                              trendline=trendline, all_same_color=all_same_color, title=title,
                              legend_on=legend_on, plot_var=plot_var, plot_hist=plot_hist, plot_reference=plot_reference,
                              gcm_list=gcm_list, rcp_list=rcp_list)
        }
      }
    }
  }
}

region_single_plot_ag <- function(xanthos_var_names, region_list, basin_list, crop_list, df_all_runs, figures_basepath, start_yr,
                                  end_yr, gcm_names, rcp_names, roll, y_ax_lbl, trendline=1){
  for(var_1 in xanthos_var_names){
    for(reg in region_list){
      for(bas in basin_list){
        for (crp in crop_list){
          for (water_type in c('irr', 'rfd')){
            if(roll==1){
              ymax_across_gcms <- max((df_all_runs %>% filter(region==reg, basin==bas, crop==crp, var==var_1))$rolling_mean)
            }else if(roll==2){
              ymax_across_gcms <- max((df_all_runs %>% filter(region==reg, basin==bas, crop==crp, var==var_1))$smoothedY)
            }else{
              ymax_across_gcms <- max((df_all_runs %>% filter(region==reg, basin==bas, crop==crp, var==var_1))$value)
            }
            for(gcm1 in gcm_names){
              fig_name <- paste0(figures_basepath, '/', var_1, "_", reg, "_", bas, "_", crp, "_", water_type, "_", gcm1,
                                 "_", if(roll==1){'rolling_mean'}else if(roll==2){'loess'}else{''}, fig_type)
              if(water_type=='irr'){
                plot_df <- df_all_runs %>%
                  filter(region==reg, basin==bas, crop==crp, gcm==gcm1, year>=start_yr, year<=end_yr, gcm %in% gcm_names,
                         rcp %in% rcp_names, var==var_1, irr==T)
              }else{
                plot_df <- df_all_runs %>%
                  filter(region==reg, basin==bas, crop==crp, gcm==gcm1, year>=start_yr, year<=end_yr, gcm %in% gcm_names,
                         rcp %in% rcp_names, var==var_1, rfd==T)
            }
              if (nrow(plot_df) > 0){
                line_plot(plot_df, fig_name, rolling=roll, y_lbl=y_ax_lbl, y_max=ymax_across_gcms, trendline=trendline)
              }
            }
          }
        }
      }
    }
  }
}

hydro_perc_change <- function(df_hydro, region_list, gcm_names, rcp_names, start_yr, end_yr, roll){
    for(reg in region_list){
      for(gcm1 in gcm_names){
        for (rcp1 in rcp_names){
          df_hydro %>% filter(reg %in% region_list, gcm==gcm1, rcp==rcp1)
        }
      }
    }
}


adjust_gcm_mean <- function(base_dir, extras_dir, level2_out_dir, time_scale, stored_in_dir, run_name,
                            xanthos_var_names, basins_filter=NULL){

  ## Purpose of file is threefold:

  # 1) Despite the fact that the ISIMIP models have been bias corrected using WATCH data,
  # this does not necessarily mean that the mean annual runoff across historical years (1970-2010) in
  # Xanthos results produced with the GCM data will be the same as the mean of the Xanthos results
  # generated using WATCH data directly. So we are having to do sort of a second bias correction here,
  # to be sure that the gcm values in history (as well as those in the future) get corrected. This is
  # particularly useful if we want to plot historical runoff (xanthos forced with watch) on the same
  # plot as GCM projections, because there is a discontinuity starting in 2010 for some GCMs.

  # 2) After performing this correction, the function then smooths the future projections with a LOESS
  # filter, but importantly does so using the historical data as part of the smoothing set, so that you get a smooth
  # continuity between historical and future data points. The function then sets historical points
  # equal to the historical mean.

  # 3) Finally, the function produces Level 2 GCAM files correctly formatted, and containing maxsubresource
  # values, corresponding to the mean annual values in historical years, but loess-smoothed values
  # in future years. To produce the xml that gcam requires, users can either run these files through
  # the model interface with the appropriate header, or can re-build the data system.

  # Authors: Sean Turner, Thomas Wild, Zarrar Khan
  # Date: April 2019

  # Required input files
  gcam_basins <- paste0(extras_dir, '/', "gcam_basin_id.csv")
  renewrsc_max_gcam <- paste0(extras_dir, '/', "L201.RenewRsrcCurves_calib_watergap.csv")
  L201.GrdRenewRsrcMax_runoff <- paste0(extras_dir, '/', "L201.GrdRenewRsrcMax_runoff.csv")

  reanalysis_data <- paste0(extras_dir, '/', "runoff_max_wfdei_1970_2010.csv")
  # get basin data for joining
  read_csv(gcam_basins) %>%
    select(basin.id, basin.name) %>% rename(GCAM_basin_name=basin.name, GCAM_basin_ID=basin.id) ->
    basin_ids

  read_csv(renewrsc_max_gcam, skip = 4) %>%
    select(region, renewresource) %>% unique() -> region_basin

  # read watch reanalysis data
  read_csv(reanalysis_data) %>%
    gather(year, runoff, -name, -id) %>%
    mutate(year = as.integer(year)) -> runoff_wfdei

  baseline_years <- 1970:2009

  # prepare watch output for GCAM (km3 per year)
  runoff_wfdei %>% group_by(id) %>%
    filter(year %in% baseline_years) %>%
    summarise(runoff = mean(runoff)) %>%
    mutate(year = 1970) %>%
    complete(year = seq(1970, 2100, 5), id) %>%
    group_by(id) %>%
    tidyr::fill(runoff) %>% ungroup() %>%
    arrange(id, year) %>%
    rename(basin.id = id,
           runoff.max = runoff) %>%
    select(basin.id, runoff.max, year) %>%
    mutate(runoff.max = round(runoff.max, 3)) %>%
    write_csv("runoff_max_noCC_wfdei.csv")

  # get wfdei mean values for baseline years (GCM deltas to be applied to these values)
  runoff_wfdei %>%
    filter(year %in% baseline_years) %>%
    group_by(id) %>%
    summarise(runoff = mean(runoff)) %>%
    mutate(gcm = "wfdei") ->
    runoff_mean_wfdei_hist

  bind_rows(
    get_gcm("GFDL-ESM2M", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("GFDL-ESM2M", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("GFDL-ESM2M", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("GFDL-ESM2M", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("HadGEM2-ES", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("HadGEM2-ES", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("HadGEM2-ES", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("HadGEM2-ES", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("IPSL-CM5A-LR", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("IPSL-CM5A-LR", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("IPSL-CM5A-LR", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("IPSL-CM5A-LR", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("MIROC-ESM-CHEM", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("MIROC-ESM-CHEM", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("MIROC-ESM-CHEM", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("MIROC-ESM-CHEM", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("NorESM1-M", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("NorESM1-M", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("NorESM1-M", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("NorESM1-M", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names)
  ) ->
    runoff_gcm_all # temp MZ

  # global comparison
  runoff_gcm_all %>%
    filter(rcp == "8p5") %>%
    # filter(rcp == "6p0") %>% # temp MZ
    select(-rcp) %>%
    bind_rows(
      runoff_wfdei %>% mutate(gcm = "wfdei") %>% select(-name)
    ) %>% group_by(id, gcm) %>% summarise(runoff = mean(runoff)) %>%
    ungroup() %>%
    spread(gcm, runoff) %>%
    ggplot(aes(`IPSL-CM5A-LR`, wfdei)) +
    geom_point() + geom_abline() +
    scale_y_log10() + scale_x_log10()
  labs(y = "Global runoff (km3)",
       title = "GCM global runoff comparison (RCP4.5)")


  # Determine the mean annual runoff in historical years (baseline yers) in the GCM runs, so they can be compared to the
  # watch data
  runoff_gcm_all %>%
    filter(year %in% baseline_years) %>%
    group_by(gcm, rcp, id) %>%
    summarise(mean_runoff = mean(runoff)) %>% ungroup() ->
    runoff_gcm_baseline_means

  runoff_gcm_all %>%
    left_join(runoff_gcm_baseline_means,
              by = c("id", "gcm", "rcp")) %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm) %>% rename(hist_mean=runoff), by=c('id')) %>%
    #mutate(delta_factor = runoff / mean_runoff) %>%
    mutate(delta_factor = mean_runoff / hist_mean) %>%
    select(id, year, gcm, rcp, delta_factor) ->
    deltas_gcm_all
  # filter(id == 36) %>%
  # ggplot(aes(year,delta_factor, colour = gcm)) +
  # geom_line() + facet_wrap(~rcp)


  # compare GCM hist means against WATCH means for same period, each plotted point represents a different basin.
  bind_rows(
    runoff_mean_wfdei_hist,
    runoff_gcm_baseline_means %>%
      filter(rcp == "2p6") %>%
      # filter(rcp == "6p0") %>% # temp MZ
      select(-rcp) %>% rename(runoff = mean_runoff)
  ) %>%
    spread(gcm, runoff) %>%
    ggplot(aes(y = wfdei)) +
    geom_point(aes(x = `IPSL-CM5A-LR`)) + geom_abline() +
    # ^^ switch gcm name in above line to view different model
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10')


  # Apply the delta factor to correct all xanthos runoff values produced with GCMs, so they reflect WATCH mean value in
  # historical years (1970-2010)
  runoff_gcm_all %>%
    left_join(deltas_gcm_all,
              by = c("id", "gcm", "rcp", "year")) %>%
    mutate(runoff_adj = delta_factor * runoff) %>%
    select(id, year, runoff_adj, gcm, rcp) ->
    runoff_gcm_all_adj
  # check adjustment--make sure that the delta correction worked by plotting one scenario.
  runoff_gcm_all_adj %>%
    filter(rcp == "4p5") %>% select(-rcp) %>%
    # filter(rcp == "6p0") %>% select(-rcp) %>% # temp MZ
    bind_rows(
      runoff_wfdei %>% mutate(gcm = "wfdei") %>% select(-name)
    ) %>% group_by(year, gcm) %>% summarise(runoff = sum(runoff_adj)) %>%
    ungroup() %>%
    ggplot(aes(year, runoff, colour = gcm)) +
    geom_line() + expand_limits(y = 0) +
    labs(y = "Global runoff (km3)",
         title = "GCM global runoff comparison (RCP4.5)")
         # title = "GCM global runoff comparison (RCP6.0)") # temp MZ

  # apply smoothing
  # get baseline period
  # Given all values from 1970-2010 are equal to the mean for all basins, you end up with only very
  # slightly smoothed values during this time, that flow nicely into the future periods.
  runoff_gcm_all_adj %>%
    left_join(runoff_mean_wfdei_hist %>% rename(baseline_mean = runoff) %>%
                select(-gcm),
              by = c("id")) %>%
    filter(year >= min(baseline_years)) %>%
    mutate(runoff_adj_basemean = if_else(year %in% baseline_years,
                                         baseline_mean, runoff_adj)) %>%
    group_by(id, gcm, rcp) %>%
    nest() %>%
    mutate(model = data %>% map(~loess(runoff_adj_basemean ~ year, data = .))) %>%
    mutate(Pred = map2(model, data, predict)) %>%
    unnest(Pred, data) %>%
    select(id, gcm, rcp, year, Pred) %>%
    rename(runoff_km3perYr = Pred) %>%
    mutate(runoff_km3perYr = if_else(runoff_km3perYr < 0, 0, runoff_km3perYr)) ->
    runoff_gcm_all_adj_smooth


  # runoff_gcm_all_adj_smooth %>%
  #   write_csv("runoff_isimip_gcam_basins.csv")


  runoff_gcm_all_adj_smooth %>%
    group_by(year, gcm, rcp) %>% summarise(runoff = sum(runoff_km3perYr)) %>%
    ungroup() %>%
    ggplot(aes(year, runoff, colour = gcm)) + geom_line() +
    facet_wrap(~rcp) + expand_limits(y = 0)

  # fix 1975 - 2010 to mean hist, so they dont apppear smoothed.
  runoff_gcm_all_adj_smooth %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm),
              by = "id") %>% rename(runoff_hist = runoff) %>%
    mutate(runoff_km3perYr = if_else(year <= 2010,
                                     runoff_hist, runoff_km3perYr)) %>%
    select(-runoff_hist) -> runoff_gcm_all_GCAM_2

  runoff_gcm_all_adj_smooth %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm),
              by = "id") %>% rename(runoff_hist = runoff) %>%
    mutate(runoff_km3perYr = if_else(year <= 2009,
                                     runoff_hist, runoff_km3perYr)) %>%
    select(-runoff_hist) -> runoff_gcm_all_GCAM_3


  # fix 1975 - 2010 to mean hist, so they dont apppear smoothed.
  runoff_gcm_all_adj_smooth %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm),
              by = "id") %>% rename(runoff_hist = runoff) %>%
    mutate(runoff_km3perYr = if_else(year <= 2010,
                                     runoff_hist, runoff_km3perYr)) %>%
    select(-runoff_hist) %>%
    filter(year %in% c(1975, 1990, seq(2005, 2095, 5), 2099)) %>%
    mutate(year = if_else(year == 2099, 2100, as.double(year))) ->
    runoff_gcm_all_GCAM


  runoff_gcm_all_GCAM %>%
    filter(id == 15) %>%
    group_by(year, gcm, rcp) %>% summarise(runoff = sum(runoff_km3perYr)) %>%
    ungroup() %>%
    ggplot(aes(year, runoff, colour = gcm)) + geom_line() +
    facet_wrap(~rcp) + expand_limits(y = 0)

  GCAM_yrs <- runoff_gcm_all_GCAM %>% .$year %>% unique()

  # prepare L2 gcam files for ISI-MIP scenarios
  # Insert these into GCAM and rebuild the data system
  gcm <- "GFDL-ESM2M"
  rcp <- "2p6"
  # gcm <- 'HadGEM2-ES' # temp MZ
  # rcp <- "6p0" # temp MZ



  # Create Level 2 csv files
  gcam_years <- c(2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065, 2070, 2075,
                  2080, 2085, 2090, 2095, 2100)
  gcm_names <- c('NorESM1-M', 'MIROC-ESM-CHEM', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'GFDL-ESM2M')
  # gcm_names <- c('HadGEM2-ES') # temp MZ
  rcp_names <- c('rcp2p6', 'rcp4p5', 'rcp6p0', 'rcp8p5')
  # rcp_names <- c('rcp6p0') # temp MZ
  variable <- 'runoff'
  #write_csv_file(runoff_gcm_all_GCAM, gcam_years, gcm_names, rcp_names, csv_basepath, variable,
  #              basin_ids=basin_ids, region_basin=region_basin, renewrsc_max_gcam=renewrsc_max_gcam, level2_out_dir=level2_out_dir,
  #              L201.GrdRenewRsrcMax_runoff=L201.GrdRenewRsrcMax_runoff)

  # Modify deltas_gcm_all to include basins so it can be  used in separate plotting module that organizes by basin.
  read_csv(gcam_basins) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(deltas_gcm_all, by='id') %>%
    select(-id) ->deltas_gcm_all
  if(!is.null(basins_filter)){
    deltas_gcm_all <- deltas_gcm_all %>% filter(name %in% basins_filter)
  }
  read_csv(gcam_basins) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(runoff_gcm_all_adj_smooth, by='id') %>%
    select(-id) %>% rename(value=runoff_km3perYr) -> runoff_gcm_all_adj_smooth
  if(!is.null(basins_filter)){
    runoff_gcm_all_adj_smooth <- runoff_gcm_all_adj_smooth %>% filter(name %in% basins_filter)
  }
  read_csv(gcam_basins) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(runoff_gcm_all_GCAM_2, by='id') %>%
    select(-id) %>% rename(value=runoff_km3perYr) ->runoff_gcm_all_GCAM_2
  if(!is.null(basins_filter)){
    runoff_gcm_all_GCAM_2 <- runoff_gcm_all_GCAM_2 %>% filter(name %in% basins_filter)
  }
  read_csv(gcam_basins) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(runoff_gcm_all_GCAM_3, by='id') %>%
    select(-id) %>% rename(value=runoff_km3perYr) ->runoff_gcm_all_GCAM_3
  if(!is.null(basins_filter)){
    runoff_gcm_all_GCAM_3 <- runoff_gcm_all_GCAM_3 %>% filter(name %in% basins_filter)
  }

  return(list('deltas_gcm_all' = deltas_gcm_all, 'runoff_gcm_all_adj_smooth' = runoff_gcm_all_adj_smooth,
              'runoff_gcm_all_GCAM_2'=runoff_gcm_all_GCAM_2, 'runoff_gcm_all_GCAM_3'=runoff_gcm_all_GCAM_3))
}

adjust_gcm_mean_country <- function(base_dir, extras_dir, level2_out_dir, country_filter, time_scale, stored_in_dir,
                                     run_name, xanthos_var_names){

  # same as adjust_gcm_mean, except with some modifications to make it work for aggregate country runoff.

  # Required input files
  gcam_countries <- paste0(extras_dir, '/', "gcam_country_id.csv")

  reanalysis_data <- paste0(extras_dir, '/', "country_runoff_max_wfdei_1970_2010.csv")
  # get basin data for joining
  read_csv(gcam_countries) %>%
    select(country.id, country.name) %>% rename(GCAM_country_name=country.name, GCAM_basin_ID=country.id) ->
    basin_ids

  # read watch reanalysis data
  read_csv(reanalysis_data) %>%
    gather(year, runoff, -name, -id) %>%
    mutate(year = as.integer(year)) -> runoff_wfdei

  baseline_years <- 1970:2009

  # prepare watch output for GCAM (km3 per year)
  runoff_wfdei %>% group_by(id) %>%
    filter(year %in% baseline_years) %>%
    summarise(runoff = mean(runoff)) %>%
    mutate(year = 1970) %>%
    complete(year = seq(1970, 2100, 5), id) %>%
    group_by(id) %>%
    tidyr::fill(runoff) %>% ungroup() %>%
    arrange(id, year) %>%
    rename(basin.id = id,
           runoff.max = runoff) %>%
    select(basin.id, runoff.max, year) %>%
    mutate(runoff.max = round(runoff.max, 3)) %>%
    write_csv("country_runoff_max_noCC_wfdei.csv")

  # get wfdei mean values for baseline years (GCM deltas to be applied to these values)
  runoff_wfdei %>%
    filter(year %in% baseline_years) %>%
    group_by(id) %>%
    summarise(runoff = mean(runoff)) %>%
    mutate(gcm = "wfdei") ->
    runoff_mean_wfdei_hist

  # gcm = "GFDL-ESM2M"; rcp = "2p6"
  # read in historical GCM values

  bind_rows(
    get_gcm("GFDL-ESM2M", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("GFDL-ESM2M", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("GFDL-ESM2M", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("GFDL-ESM2M", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("HadGEM2-ES", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("HadGEM2-ES", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("HadGEM2-ES", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("HadGEM2-ES", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("IPSL-CM5A-LR", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("IPSL-CM5A-LR", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("IPSL-CM5A-LR", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("IPSL-CM5A-LR", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("MIROC-ESM-CHEM", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("MIROC-ESM-CHEM", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("MIROC-ESM-CHEM", "6p0", base_dir, stored_in_dir, run_name, time_scale), xanthos_var_names, get_gcm("MIROC-ESM-CHEM", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("NorESM1-M", "2p6", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("NorESM1-M", "4p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names),
    get_gcm("NorESM1-M", "6p0", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names), get_gcm("NorESM1-M", "8p5", base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names)
  ) ->
    runoff_gcm_all

  # Determine the mean annual runoff in historical years (baseline yers) in the GCM runs, so they can be compared to the
  # watch data
  runoff_gcm_all %>%
    filter(year %in% baseline_years) %>%
    group_by(gcm, rcp, id) %>%
    summarise(mean_runoff = mean(runoff)) %>% ungroup() ->
    runoff_gcm_baseline_means

  runoff_gcm_all %>%
    left_join(runoff_gcm_baseline_means,
              by = c("id", "gcm", "rcp")) %>%
    mutate(delta_factor = runoff / mean_runoff) %>%
    select(id, year, gcm, rcp, delta_factor) ->
    deltas_gcm_all


  # Apply the delta factor to correct all xanthos runoff values produced with GCMs, so they reflect WATCH mean value in
  # historical years (1970-2010)
  deltas_gcm_all %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm),
              by = c("id")) %>%
    mutate(runoff_adj = delta_factor * runoff) %>%
    select(id, year, runoff_adj, gcm, rcp) ->
    runoff_gcm_all_adj


  # apply smoothing
  # get baseline period
  # Given all values from 1970-2010 are equal to the mean for all basins, you end up with only very
  # slightly smoothed values during this time, that flow nicely into the future periods.
  runoff_gcm_all_adj %>%
    left_join(runoff_mean_wfdei_hist %>% rename(baseline_mean = runoff) %>%
                select(-gcm),
              by = c("id")) %>%
    filter(year >= min(baseline_years)) %>%
    mutate(runoff_adj_basemean = if_else(year %in% baseline_years,
                                         baseline_mean, runoff_adj)) %>%
    group_by(id, gcm, rcp) %>%
    nest() %>%
    mutate(model = data %>% map(~loess(runoff_adj_basemean ~ year, data = .))) %>%
    mutate(Pred = map2(model, data, predict)) %>%
    unnest(Pred, data) %>%
    select(id, gcm, rcp, year, Pred) %>%
    rename(runoff_km3perYr = Pred) %>%
    mutate(runoff_km3perYr = if_else(runoff_km3perYr < 0, 0, runoff_km3perYr)) ->
    runoff_gcm_all_adj_smooth


  runoff_gcm_all_adj_smooth %>%
    group_by(year, gcm, rcp) %>% summarise(runoff = sum(runoff_km3perYr)) %>%
    ungroup() %>%
    ggplot(aes(year, runoff, colour = gcm)) + geom_line() +
    facet_wrap(~rcp) + expand_limits(y = 0)

  # fix 1975 - 2010 to mean hist, so they dont apppear smoothed.
  runoff_gcm_all_adj_smooth %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm),
              by = "id") %>% rename(runoff_hist = runoff) %>%
    mutate(runoff_km3perYr = if_else(year <= 2010,
                                     runoff_hist, runoff_km3perYr)) %>%
    select(-runoff_hist) -> runoff_gcm_all_GCAM_2

  runoff_gcm_all_adj_smooth %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm),
              by = "id") %>% rename(runoff_hist = runoff) %>%
    mutate(runoff_km3perYr = if_else(year <= 2009,
                                     runoff_hist, runoff_km3perYr)) %>%
    select(-runoff_hist) -> runoff_gcm_all_GCAM_3


  # fix 1975 - 2010 to mean hist, so they dont appear smoothed.
  runoff_gcm_all_adj_smooth %>%
    left_join(runoff_mean_wfdei_hist %>% select(-gcm),
              by = "id") %>% rename(runoff_hist = runoff) %>%
    mutate(runoff_km3perYr = if_else(year <= 2010,
                                     runoff_hist, runoff_km3perYr)) %>%
    select(-runoff_hist) %>%
    filter(year %in% c(1975, 1990, seq(2005, 2095, 5), 2099)) %>%
    mutate(year = if_else(year == 2099, 2100, as.double(year))) ->
    runoff_gcm_all_GCAM

  # Modify deltas_gcm_all to include basins so it can be  used in separate plotting module that organizes by basin.
  read_csv(gcam_countries) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(deltas_gcm_all, by='id') %>%
    select(-id) %>% filter(name %in% basins_filter) -> deltas_gcm_all
  read_csv(gcam_countries) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(runoff_gcm_all_adj_smooth, by='id') %>%
    select(-id) %>% filter(name %in% basins_filter) %>% rename(value=runoff_km3perYr) ->runoff_gcm_all_adj_smooth
  read_csv(gcam_countries) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(runoff_gcm_all_GCAM_2, by='id') %>%
    select(-id) %>% filter(name %in% basins_filter) %>% rename(value=runoff_km3perYr) ->runoff_gcm_all_GCAM_2
  read_csv(gcam_countries) %>% rename(id=basin.id) %>% rename(name=basin.name) %>% left_join(runoff_gcm_all_GCAM_3, by='id') %>%
    select(-id) %>% filter(name %in% basins_filter) %>% rename(value=runoff_km3perYr) ->runoff_gcm_all_GCAM_3
  return(list('deltas_gcm_all' = deltas_gcm_all, 'runoff_gcm_all_adj_smooth' = runoff_gcm_all_adj_smooth,
              'runoff_gcm_all_GCAM_2'=runoff_gcm_all_GCAM_2, 'runoff_gcm_all_GCAM_3'=runoff_gcm_all_GCAM_3))
}

write_csv_file <- function(input, gcam_years, gcms, rcps, csv_basepath, variable,
                           basin_ids=NULL, region_basin=FALSE, renewrsc_max_gcam=FALSE, level2_out_dir=FALSE,
                           L201.GrdRenewRsrcMax_runoff=FALSE, gcam_xanthos_basin_mapping=NULL, name_exclude_list=NULL){
  if(!is.null(gcam_xanthos_basin_mapping)){
    mapping_file <- read.csv(gcam_xanthos_basin_mapping)
  }
  iso_GCAM_reg_ID <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output/iso_GCAM_regID.csv'
  GCAM_basin_country <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output/GCAMBasin_country.csv'
  gcam_region_names <- 'E:/NEXO-UA/Results/downscaling/water/xanthos/output/GCAM_region_names.csv'
  # Process data
  for(clim_mod in gcms){
    for (forc in rcps){
      if(variable=='hydro'){
        export_df <- input %>% select(name, year, var, gcm, rcp, clim_imp_val) %>%
          rename(region=name, fixedOutput=clim_imp_val) %>% # period=year,
          filter(year %in% gcam_years, gcm == clim_mod, rcp==forc) %>% select(-var, -gcm, -rcp)  # period
        export_df$supplysector <- 'electricity'
        export_df$subsector <- 'hydro'
        col_name <- c('stub-technology')
        export_df[[col_name]] <- 'hydro'
        col_name <- c('share.weight.year')  # 'share-weight'
        export_df[[col_name]] <- 0
        export_df <- export_df %>% mutate(share.weight.year=year)
        export_df$subs.share.weight <- 0
        export_df$tech.share.weight <- 0
        # save file
        fileName <- paste0("hydro_impacts", "_", clim_mod, "_", forc, ".csv")
        write.table('INPUT_TABLE,,,,,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE,
                    sep = ',', quote=FALSE)
        write.table('Variable ID,,,,,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE,
                    append=TRUE, quote=FALSE)
        write.table('StubTechFixOut,,,,,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE,
                    append = TRUE, quote=FALSE)
        write.table(',,,,,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE, append = TRUE,
                    quote=FALSE)
        write.table(export_df, paste0(csv_basepath, '/', fileName), append=TRUE, col.names=TRUE, row.names = FALSE,
                    sep=',', quote=FALSE)
      }else if(variable=='runoff'){
        export_df <- input %>% filter(gcm == clim_mod, rcp == forc, year %in% gcam_years) %>%
          rename(year.fillout=year) %>%
          mutate(sub.renewable.resource='runoff')
        if(!is.null(gcam_xanthos_basin_mapping)){
          export_df <- export_df %>%
            left_join(mapping_file, by=c('name')) %>%
            select(-name) %>%
            rename(name = GCAM.basin.name)
        }
        f1 <- read.csv(iso_GCAM_reg_ID, skip=6)
        f2 <- read.csv(GCAM_basin_country)
        f3 <- f1 %>% left_join(f2, by="country") %>% select(-iso, -region_GCAM3)
        f4 <- read.csv(gcam_region_names, skip=6)
        f5 <- f3 %>% left_join(f4, by=c('GCAM_region_ID')) %>% select(-country, GCAM_region_ID)
        export_df <- export_df %>%
          left_join(f5, by=c('name')) %>%
          mutate(renewresource = paste0(name, "_water withdrawals")) %>%
#          left_join(region_basin, by=c('renewresource')) %>%
          rename(maxSubResource=clim_imp_val) %>%
#          mutate(renewresource=str_replace_all(renewresource, '-water withdrawals', '_water withdrawals')) %>%
          filter(region!="NA") %>%
          filter(!name %in% name_exclude_list) %>%
          select(region, renewresource, sub.renewable.resource, year.fillout, maxSubResource)

        renewrsc_max_gcam <- paste0(extras_dir, '/', "L201.RenewRsrcCurves_calib_watergap.csv")

        # save file
        fileName <- paste0("runoff_impacts", "_", clim_mod, "_", forc, ".csv")
        write.table('INPUT_TABLE,,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE,
                    sep = ',', quote=FALSE)
        write.table('Variable ID,,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE,
                    append=TRUE, quote=FALSE)
        write.table('GrdRenewRsrcMax,,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE,
                    append = TRUE, quote=FALSE)
        write.table(',,,,,', file=paste0(csv_basepath, '/', fileName), row.names = FALSE, col.names=FALSE, append = TRUE,
                    quote=FALSE)
        write.table(export_df, paste0(csv_basepath, '/', fileName), append=TRUE, col.names=TRUE, row.names = FALSE,
                    sep=',', quote=FALSE)


      }else if(variable=='agProd'){
        print("no code yet")
      }
    }
  }

}

csv2xml <- function(csvpath, xmlpath, gcam_variable){
  options("gcamdata.use_java"=TRUE)
  filelist <- list.files(path=csvpath, full.names=TRUE, recursive=FALSE)
  filelist <- filelist[(!grepl(".xml", filelist))]
  invisible(foreach(f = filelist) %do% {
    filename1 <- substr(f, (nchar(csvpath)+2), (nchar(f)-4))
    tibble::as.tibble(read.csv(f, skip = 4, stringsAsFactors = F)) -> x
    gcamdata::create_xml(paste0(xmlpath, "/", filename1, ".xml"), mi_header='E:/NEXO-UA/Results/downscaling/water/xanthos/output/ModelInterface_headers.txt') %>%
      gcamdata::add_xml_data(x, gcam_variable) %>%
      gcamdata::run_xml_conversion()
  })
}


get_gcm <- function(gcm, rcp, base_dir, stored_in_dir, run_name, time_scale, xanthos_var_names){
  # run_name_2 <- paste0(run_name, "_", gcm, "_", 'rcp', rcp, "_", time_scale)
  run_name_2 <- paste0(run_name, "_", gcm, "_", 'rcp', rcp) # temp MZ
  if (stored_in_dir==1){
    xanthos_dir <- paste0(base_dir, '/', run_name_2)
  }else{
    xanthos_dir <- paste0(base_dir)
  }
  xanthos_file <- paste0(xanthos_var_names, "_", gcm, "_", 'rcp', rcp, '_', time_scale, '.csv')
  xanthos_output_filepath <- paste0(xanthos_dir, '/', xanthos_file)
  read_csv(xanthos_output_filepath) %>%
    gather(year, runoff, -name, -id) %>%
    select(-name) %>%
    mutate(gcm = gcm, rcp = rcp, year = as.integer(year)) %>%
    mutate(rcp = paste0('rcp', rcp))
}
