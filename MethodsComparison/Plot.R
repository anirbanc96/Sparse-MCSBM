library(readr)
library(tidyverse)
ValsAMP <- read.csv("~/Dropbox (Penn)/Finite degree Multi-Layer SBM/IMA Revision Codes/MethodsComparison/AMP/ValsCompAMP.csv")[,-1]
ValsBP <- read.csv("~/Dropbox (Penn)/Finite degree Multi-Layer SBM/IMA Revision Codes/MethodsComparison/BP/ValsCompRatio.csv")[,-1]
ValsDCMASE <- read.csv("~/Dropbox (Penn)/Finite degree Multi-Layer SBM/IMA Revision Codes/MethodsComparison/DCMASE/ValsCompDCMASE.csv")[-1,-1]

names(ValsDCMASE) <- names(ValsAMP)

thresh <- seq(0, 4.5, 0.5)

for (i in 1:3){

  ratio.data <- rbind(ValsAMP[i,], ValsBP[i,], ValsDCMASE[i,])
  ratio.data[is.na(ratio.data)] <- 0
  ratio.data <- as.data.frame(t(rbind(thresh, ratio.data)))
  colnames(ratio.data) <- c("thershold", "AMP", "BP", "DC-MASE")
  
  comp.plot <- ratio.data %>%
    ggplot(aes(x = thershold, y = AMP)) + 
    geom_vline(xintercept = 1, linetype = "dotted", size = 0.5) +
    geom_line(aes(y = AMP, color = "AMP")) +
    geom_point((aes(y = AMP, color = "AMP"))) +
    geom_line(aes(y = BP, color = "BP")) +
    geom_point(aes(y = BP, color = "BP")) +
    geom_line(aes(y = `DC-MASE`, color = "DC-MASE")) +
    geom_point(aes(y = `DC-MASE`, color = "DC-MASE")) +
    labs(title = "Empirical Overlap", x = "Combined Signal Strength",
         y = "Overlap", color = "Methods") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.22, 0.73),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  file_name = paste("comapre_", i, ".pdf", sep="")
  
  ggsave(plot = comp.plot, filename = file_name, device = "pdf",
         width = 4, height = 3)
}

