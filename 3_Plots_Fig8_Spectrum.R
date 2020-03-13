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

for(var in c("ISOT", "ITPC")){
  pdf(file = paste0("Plots/ECHAM7k_7_Spectra_",var,".pdf"), width = 8*2/3, height = 5*1.2*2/3)
  #png(file = paste0("Plots/ECHAM7k_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$full$spec, col = adjustcolor("#08306b", alpha = 0.5), 
        ylim = c(0.000005,1000), 
        xlim = c(1/3000, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$full$spec), col = "#08306b", lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 2)
  mtext("Power spectral sensity", side = 2, line= 2)
  #text(0.2, 8e2, "ECHAM5 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.0005,0.001,0.002,0.0033333333, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), 
       labels = c(1/0.0005, 1/0.001, 1/0.002, 300, 1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5))
  axis(side = 2, at = c(1e-3, 1e-1, 1e1, 1e3), 
       labels = c(1e-3, 0.1, 10, 1000))
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$ds$spec, col = adjustcolor("#91002B", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$ds$spec), col = "#91002B", lw = 2)
  #text(0.02, 0.01, "ECHAM5 down-sampled", col = "#91002B")
  #text(0.02, 0.001, "10y filter", col = "#91002B")
  
  LLines(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec, col = adjustcolor("black", alpha = 0.5), lw = 2)
  LLines(LogSmooth(ANALYSIS$SPECTRA$MEAN_SPEC$Record$spec), col = "black", lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  
  if(var == "ISOT"){
    legend("bottomleft", legend = c("ECHAM5 (full)", "ECHAM5 (down-sampled)", "Records"), 
           col = c("#08306b","#91002B","black"), lwd = c(1,1,2), lty = c(1,1,1), bty = "n") 
    text(0.5, 500, TeX("ECHAM5: $\\delta^{18}O$ in prec"), adj = 1)
  }
  if(var == "ITPC"){
    legend("bottomleft", legend = c("ECHAM5 (full)", " ...  3y filter", " ... 12y filter", "ECHAM5 (down-sampled)", " ... 4y filter", "Records"), 
           col = c("#08306b","#4292c6","#c6dbef","#91002B","#91002B","black"), lwd = c(1,1,1,1,1,2), lty = c(1,1,1,1,3,1), bty = "n") 
    text(0.5, 500, TeX("ECHAM5: prec-weighted $\\delta^{18}O$"), adj = 1)
  }
  
  
  dev.off()
  
}

