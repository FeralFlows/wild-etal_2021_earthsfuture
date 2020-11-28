##################################################################
# Plot CO2 emission under Reference and Policy Scenarios         #
# Author: Mengqi Zhao                                      #
# Email: mengqiz@umd.edu                                         #
# Last Update: 2020-10-26                                        #  
##################################################################

library(ggplot2)

# Argentina: Plot Reference case CO2 emissions and Policy Co2 emissions CO2 emissions.
CO2_Emissions_RefPolicy_Argentina <- read.csv("E:/NEXO-UA/Results/metis/Figure_EmissionsGoalsLine.csv", skip=1) # for Argentina
CO2_Emissions_RefPolicy_ROW <- read.csv("E:/NEXO-UA/Results/metis/Figure_EmissionsGoalsLine_ROW.csv", skip=1) # for rest of the world

emission_plot <- function(CO2_Emissions_RefPolicy, region_name, y_lim){
  CO2_Emissions_RefPolicy <- CO2_Emissions_RefPolicy %>%
    gather(Scenario, value, Reference:Policy)
  p <- ggplot(data = CO2_Emissions_RefPolicy %>% filter(!Scenario=='NDC'), mapping = aes(x = year, y = value, color = Scenario, fill=Scenario))
  p <- p + scale_color_manual(values=c("#1bab55", "black"))
  p <- p + ylim(y_lim)
  p <- p + ylab(expression(CO[2]~Emissions~(10^6~tons)))
  p <- p + xlab('Year')
  p <- p + geom_line(size=1)
  p <- p + theme(text = element_text(size=12), axis.text.x = element_text(size=12))
  p <- p + theme_bw()
  dirOutputs <- 'E:/NEXO-UA/Results/metis/outputs'
  fname <- paste('Figure_emissCO2_RefPolicyCompare', region_name, sep = '_')
  pdfpng <- 'pdf'
  figWidth <- 4
  figHeight <- 3
  p
  metis.printPdfPng(figure=p,
                    dir=dirOutputs,
                    filename=fname,
                    figWidth=figWidth,
                    figHeight=figHeight,
                    pdfpng=pdfpng)
}


emission_plot(CO2_Emissions_RefPolicy_Argentina, 'Argentina', c(0,300))
emission_plot(CO2_Emissions_RefPolicy_ROW, 'ROW', c(0,60000))
