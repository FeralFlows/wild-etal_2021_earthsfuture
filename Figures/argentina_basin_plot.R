# This is to plot the basins within Argentina

if("devtools" %in% rownames(installed.packages()) == F){install.packages("devtools")}
library(devtools)
if("metis" %in% rownames(installed.packages()) == F){devtools::install_github(repo="JGCRI/metis")}
library(metis)
if('sp' %in% rownames(installed.packages()) == F){install.packages('sp')}
library(sp)

# Create Basins and hydrologic watersheds within Argentina (better solution)
m_basin <- metis::mapGCAMBasins
m_country <- metis::mapCountries

m_argentina <- m_country[m_country@data$subRegion %in% c("Argentina"),]; sp::plot(m_argentina)

# Basins within Argentina
m_argentina_basin <- sp::spTransform(m_argentina, raster::crs(m_basin))
m_argentina_basin <- raster::crop(m_basin, m_argentina_basin)
m_argentina_basin@data <-  droplevels(m_argentina_basin@data)
sp::plot(m_argentina_basin)

# Plot Argentina Maps at basin and hydroshed level
metis.map(m_argentina_basin, fileName = 'Argentina_basin', labels = T, labelsSize = 0.6)
