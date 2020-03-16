#################################################
## Paper Figure SPECTRUM ########################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

#################################################

## PLOT

for(var in c("TEMP","PREC","ISOT", "ITPC")){
  pdf(file = paste0("Plots/COMPARE_Plot_7_Spectra_",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH*2/3, height = PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  #png(file = paste0("Plots/ECHAM_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor("#08306b", alpha = 0.5), 
        ylim = c(0.000005,10000), xlim = c(1/300, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#08306b", lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 2)
  mtext("Power spectral sensity", side = 2, line= 2)
  LLines(ANALYSIS_HadCM3$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = "#4292c6", lw = 2)
  LLines(LogSmooth(ANALYSIS_HadCM3$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#4292c6", lw = 2)
  #LPlot(LogSmooth(ANALYSIS_HadCM3$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#4292c6")
  #text(0.2, 8e2, "ECHAM5 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
       labels = c(1/0.002, 300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))
  axis(side = 2, at = c(1e-3, 1e-1, 1e1, 1e3), 
       labels = c(1e-3, 0.1, 10, 1000))
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor("#91002B", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = "#91002B", lw = 2)
  
  LLines(ANALYSIS_HadCM3$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor("salmon", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS_HadCM3$SPECTRA$MEAN_SPEC[[paste0("SIM_ds_",var)]]$spec), col = "salmon", lw = 2)
  
  #text(0.02, 0.01, "ECHAM5 down-sampled", col = "#91002B")
  #text(0.02, 0.001, "10y filter", col = "#91002B")
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor("black", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = "black", lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  
  if(var == "ISOT"){
    legend("bottomleft", legend = c("ECHAM5 (full)", "HadCM3 (full)", "ECHAM5 (down-sampled)", "HadCM3 (down-sampled)", "Records"), 
           col = c("#08306b","#4292c6","#91002B","salmon","black"), lwd = c(1,1,1,1,2), lty = c(1,1,1,1,1), bty = "n") 
    text(0.5, 500, TeX("$\\delta^{18}O$ in prec"), adj = 1)
  }
  if(var == "ITPC"){
    legend("bottomleft", legend = c("ECHAM5 (full)", "HadCM3 (full)", "ECHAM5 (down-sampled)", "HadCM3 (down-sampled)", "Records"), 
           col = c("#08306b","#4292c6","#91002B","salmon","black"), lwd = c(1,1,1,1,2), lty = c(1,1,1,1,1), bty = "n") 
    text(0.5, 500, TeX("prec-weighted $\\delta^{18}O$"), adj = 1)
  }
  
  
  dev.off()
  
}

# LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = adjustcolor("#08306b", alpha = 0.5), 
#       ylim = c(0.000005,10000), xlim = c(1/300, 0.5),
#       ylab = "",
#       xaxt = 'n',
#       yaxt = "n",
#       xlab = "", lwd = 2)
# LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#08306b", lw = 3)#,
# #main = TeX("Mean Spectra from cave locations (res>8)"))
# mtext("Period (y)", side = 1, line= 2)
# mtext("Power spectral sensity", side = 2, line= 2)
# LLines(ANALYSIS_HadCM3$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = "#4292c6", lw = 2)
# LLines(LogSmooth(ANALYSIS_HadCM3$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec), col = "#4292c6", lw = 2)
# 
# LLines(MeanSpectrum(list(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY14,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY21,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY33,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY48,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY51,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY76,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY85,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY93,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY94,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY95,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY97))$spec, ylim = c(0.001,1000))
# 
# #LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_",var)]]$spec, col = "green", lw = 3)
# 
# LLines(MeanSpectrum(list(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY111,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY113,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY117,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY118,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY123,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY128,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY137,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY147,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY165,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY172,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY178,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY179,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY185,
#                         ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY187))$spec, col = "red")
# 
# LLines(MeanSpectrum(list(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY201,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY202,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY212,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY222,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY226,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY238,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY240,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY242,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY253,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY278,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY286,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY289,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY294,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY296))$spec, col = "darkgreen")
# 
# LLines(MeanSpectrum(list(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY305,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY326,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY329,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY330,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY335,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY351,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY358,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY361,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY378,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY390,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY395,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY399))$spec, col = "yellow")
# 
# LLines(MeanSpectrum(list(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY420,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY422,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY430,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY435,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY442,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY443,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY447,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY448,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY461,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY464,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY496,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY498))$spec, col = "orange")
# 
# LLines(MeanSpectrum(list(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY506,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY514,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY522,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY528,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY538,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY539,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY541,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY544,
#                          # ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY546,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY547,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY560,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY564))$spec, col = "orange")
# 
# 
# 
# LLines(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY548, col = "blue")
# 
# LLines(MeanSpectrum(list(ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY613,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY620,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY621,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY623,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY672,
#                          ANALYSIS$SPECTRA$SIM_full$a$ITPC$ENTITY673))$spec, col = "blue")
# 
#git comment