r_packages <- c("ggplot2",
                "ggrepel",
                "ggpubr",
                "hrbrthemes",
                "scales",
                "patchwork",
                "poweRlaw",
                "maps",
                "geosphere",
                "CoordinateCleaner",
                "dplyr",
                "reshape2",
                "FAOSTAT",
                "utils")

for (pkg in r_packages){
  if (!require(pkg)) {
    install.packages(pkg)
    library(pkg)
  }
}

for (pkg in r_packages){
  print(pkg)
}
