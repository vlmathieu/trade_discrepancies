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
  file = "/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/market_concentration.csv", # nolint
  header = TRUE,
  sep = ";")
head(data)
data[data$period == 2022 & data$cmd == 12, ]
select(filter(data, period == 2022, cmd == 12), traded_value_exp)

pal <- c("#3D85F7", "#C32E5A")
order <- c("hhi_imp", "hhi_exp")

plot_mkt <- data %>%
  filter(cmd == 12) %>%
  select(c("period", "hhi_imp", "hhi_exp")) %>%
  melt(id.vars = "period",
       variable.name = "trader_type",
       value.name = "hhi") %>%
  arrange(period) %>%
  mutate(trader_type = factor(trader_type, levels = order)) %>%
  ggplot(aes(x = period,
             y = hhi,
             color = trader_type,
             label = trader_type)) +
  geom_line() +

  # End of chart labels
  annotate("text", 
           x = max(data$period) + 0.4,
           y = data[data$period == max(data$period) &
                      data$cmd == prod, ]$hhi_imp,
           label = "Imports",
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = pal[1]) +

  annotate("text",
           x = max(data$period) + 0.4,
           y = data[data$period == max(data$period) &
                      data$cmd == prod, ]$hhi_exp,
           label = "Exports",
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = pal[2]) +

  # Segments
  geom_segment(aes(x = min(data$period),
                   y = 0.25,
                   xend = max(data$period),
                   yend = 0.25),
               color = "black",
               linetype = "dotted") +
  annotate("text",
           x = max(data$period) + 0.4,
           y = 0.25,
           label = "HHI = 0.25",
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = min(data$period),
                   y = 0.15,
                   xend = max(data$period),
                   yend = 0.15),
               color = "black",
               linetype = "dotted") +
  annotate("text",
           x = max(data$period) + 0.4,
           y = 0.15,
           label = "HHI = 0.15",
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = min(data$period),
                   y = 0.01,
                   xend = max(data$period),
                   yend = 0.01),
               color = "black",
               linetype = "dotted") +
  annotate("text",
           x = max(data$period) + 0.4,
           y = 0.01,
           label = "HHI = 0.01",
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  scale_color_manual(values = pal) +
  scale_x_continuous(
    breaks = c(1996, 2000, 2005, 2010, 2015, 2020, 2022),
    labels = c("1996", "2000", "2005", "2010", "2015", "2020", "2022")
  ) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
  labs(x = "Year",
       y = "Herfindahl-Hirschman Index") +
  coord_cartesian(clip = "off") +
  theme_ipsum(axis_title_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(40, 80, 20, 20))

ggsave(
  "mkt_concentration.png",
  plot = plot_mkt,
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

