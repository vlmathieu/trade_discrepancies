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
  file = "/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/contributor_profiles.csv", # nolint
  header = TRUE,
  sep = ";")
head(data)
to_match <- paste0("primary_value", c("_imp", "_exp"))
size_col <- colnames(data)[grep(paste(to_match, collapse = "|"), colnames(data))]

pal <- c("#3D85F7", "#C32E5A")
columns <- c("country",
             "nb_edge_imp",
             "nb_edge_exp",
             "primary_value_imp",
             "primary_value_exp")

plot_profile <- data %>%
  filter(cmd == 12, period == 2001) %>%
  select(all_of(columns)) %>%
  mutate(across(all_of(size_col), ~ .x / 1000000)) %>%
  arrange(desc(size_col[grep("_imp", size_col)])) %>%
  setNames(gsub("primary_value", "size", names(.))) %>%
  mutate(country = factor(country)) %>%

  ggplot(aes(x = nb_edge_exp,
             y = nb_edge_imp,
             label = country)) +

  geom_segment(aes(x = 0,
                   y = 0,
                   xend = ceiling(max(data$nb_edge_exp) / 10) * 10,
                   yend = ceiling(max(data$nb_edge_imp) / 10) * 10),
               color = "grey",
               linetype = "dotted") +

  geom_point(aes(size = size_exp), alpha = 0.6, color = pal[2]) +
  geom_point(aes(size = size_imp), alpha = 0.6, color = pal[1]) +

  geom_text_repel(size = 2,
                  min.segment.length = 1.2,
                  position = position_nudge_repel(x = 3, y = 3),
                  fontface = "bold") +
  geom_label(label = 2022, x = 114, y = 4, size = 5) +

  scale_size(range = c(1, 20),
             breaks = c(1000, 2000, 5000, 10000),
             limits = c(0, 20000),
             name = "Traded value (million US$)") +
  scale_x_continuous(name = "Number of export partners",
                     limits = c(0, ceiling(max(data$nb_edge_exp) / 10) * 10),
                     breaks = seq(0,
                                  ceiling(max(data$nb_edge_exp) / 10) * 10,
                                  by = 20),
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Number of import partners",
                     limits = c(0, 120),
                     breaks = seq(20,
                                  ceiling(max(data$nb_edge_exp) / 10) * 10,
                                  by = 20),
                     expand = c(0, 0)) +
  theme_ipsum(axis_title_size = 11)

ggsave(
  "contributor_profiles.png",
  plot = plot_profile,
  device = "png",
  path = "/Users/valentinmathieu/Desktop/",
  scale = 1,
  width = 2480 * 1.1,
  height = 1240 * 1.5,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = "white"
)
