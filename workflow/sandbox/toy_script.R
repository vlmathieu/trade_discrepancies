r_packages <- c("ggplot2",
                "ggh4x",
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
  }
}

library("ggplot2")
library("dplyr")
library("hrbrthemes")
library("reshape2")
library("ggrepel")

library("ggh4x")
library("ggrepel")
library("ggpubr")
library("scales")
library("patchwork")
library("grid")
library("tidyverse")
library("ggstream")
library("showtext")
library("ggtext")

data <- read.csv(
  file = "/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/degree_concentration.csv", # nolint
  header = TRUE,
  sep = ";")
head(data)
data$trade_degree_imp
summary(data)

pal <- c("#3D85F7", "#C32E5A")

lst <- gsub("_exp", "", names(data)[grep("exp", names(data))])
metrics <- lst[lst != "nb_edge"]

myplots <- list()

for (metric in metrics) {

  # Select columns for plot
  columns <- c("period", paste0(metric, c("_imp", "_exp")))

  # Define order for corresponding colors
  order <- paste0(metric, c("_imp", "_exp"))

  # Define name for y axis
  if ("mean" %in% metric){
    name <- "Mean"
  } else if ("var" %in% metric) {
    name <- "Variance"
  } else if ("skew" %in% metric) {
    name <- "Skewness"
  } else if ("kurt" %in% metric) {
    name <- "Kurtosis"
  }

  # Produce basic plot
  plot_connectivity <- data %>%
    filter(cmd == 12) %>%
    # Select relevant columns
    select(all_of(columns)) %>%
    # Convert wide to long dataframe
    melt(id.vars = "period",
         variable.name = "trade_flow",
         value.name = "metric") %>%
    arrange(period) %>%
    # Convert trader_type to factor and order its values
    mutate(trade_flow = factor(trade_flow, levels = order)) %>%

    # Create ggplot line plot
    ggplot(aes(x = period,
               y = metric,
               color = trade_flow,
               label = trade_flow)) +
    geom_point(alpha = 0.6) +
    # geom_line(alpha = 0.6) +
    geom_smooth(method = NULL,
                se = FALSE,
                linewidth = 1)

  # Collect final y values of geom_smooth to plot legend label
  y_imp <- filter(ggplot_build(plot_connectivity)$data[[2]],
                  x == max(data$period) & grepl("_imp", label))$y
  y_exp <- filter(ggplot_build(plot_connectivity)$data[[2]],
                  x == max(data$period) & grepl("_exp", label))$y

  # Collect y_min and y_max for setting axis limits
  y_min <- min(min(ggplot_build(plot_connectivity)$data[[2]]$y),
               min(select(filter(data, cmd == 12), contains(metric))))
  y_max <- max(max(ggplot_build(plot_connectivity)$data[[2]]$y),
               round(max(select(filter(data, cmd == 12), contains(metric)))))

  # Complete plot with annotation and theme
  plot_connectivity <- plot_connectivity +
    # End of chart labels
    annotate("text",
             x = max(data$period) + 0.4,
             y = y_imp,
             label = "Imports",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = pal[1]) +

    annotate("text",
             x = max(data$period) + 0.4,
             y = y_exp,
             label = "Exports",
             hjust = 0,
             size = 3,
             lineheight = .8,
             fontface = "bold",
             color = pal[2]) +

    # Set up scale colors, breaks, and limits + themes
    scale_color_manual(values = pal) +
    scale_x_continuous(
      breaks = c(1996, 2000, 2005, 2010, 2015, 2020, 2022),
      labels = c("1996", "2000", "2005", "2010", "2015", "2020", "2022")
    ) +
    scale_y_continuous(limits = c(trunc(y_min - 1),
                                  round(y_max + 1)),
                       expand = c(0, 0)) +
    labs(x = "Year",
         y = paste0(name,
                    " of the nb. of trading partners")) +
    coord_cartesian(clip = "off") +
    theme_ipsum(axis_title_size = 11) +
    theme(legend.position = "none",
          plot.margin = margin(20, 20, 20, 20))

  # Store plot in list
  myplots[[match(metric, metrics)]] <- plot_connectivity
}

# Produce composite plot using patchwork
final_plot <- myplots[[1]] + myplots[[2]] + myplots[[3]] + myplots[[4]] +
  plot_layout(ncol = 2)

ggsave(
  "connectivity.png",
  plot = final_plot,
  device = "png",
  path = "/Users/valentinmathieu/Desktop/",
  scale = 1,
  width = 2480 * 1.5,
  height = 1240 * 2,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = "white"
)
