# Run experiments over N
setwd("~/Downloads/Sparse-MCSBM-main/Over N")
source("Body.R")
source("Plot.R")

rm(list = ls())

# Run experiments over M
setwd("~/Downloads/Sparse-MCSBM-main/Over M")
source("Body.R")
source("Plot.R")

rm(list = ls())

# Run experiments over ratio
setwd("~/Downloads/Sparse-MCSBM-main/OverRatio")
source("Body.R")
source("RatioPlot.R")

rm(list = ls())

# Run experiments for comparison
setwd("~/Downloads/Sparse-MCSBM-main/MethodsComparison/BP")
source("BodyComp.R")
rm(list = ls())
setwd("~/Downloads/Sparse-MCSBM-main/MethodsComparison/AMP")
source("AMP.R")
rm(list = ls())
setwd("~/Downloads/Sparse-MCSBM-main/MethodsComparison/DCMASE")
source("Overlap.R")
rm(list = ls())
setwd("~/Downloads/Sparse-MCSBM-main/MethodsComparison")
source("Plot.R")