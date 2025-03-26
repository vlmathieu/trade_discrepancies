# for (pkg in list(snakemake@params[["r_packages"]])){
# #   if (!require(pkg)) {
# #     install.packages(pkg, repos = "http://cran.us.r-project.org")
# #   }
#   library(pkg)
# }

library("ggplot2")
library("ggrepel")
library("ggpubr")
library("hrbrthemes")
library("scales")
library("patchwork")
# # library("poweRlaw")
# library("maps")
# library("geosphere")
# # library("CoordinateCleaner")
library("dplyr")
# library("reshape2")
# # library("FAOSTAT")
# library("utils")

ip <- as.data.frame(installed.packages()[, c(1, 3:4)])
ip <- ip[is.na(ip$Priority), 1:2, drop = FALSE]
print(ip)
