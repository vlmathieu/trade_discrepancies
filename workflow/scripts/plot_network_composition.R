library("ggplot2")
library("hrbrthemes")
library("dplyr")
library("reshape2")

data <- read.csv(file = snakemake@input[[1]],
                 header = TRUE,
                 sep = ";")

pal <- c("#3D85F7", "#A49393", "#C32E5A")
order <- c("nb_pure_imp", "nb_mixed", "nb_pure_exp")

plot_compo <- data %>%
  filter(cmd == 12) %>%
  select(c("period", "nb_pure_exp", "nb_pure_imp", "nb_mixed")) %>%
  melt(id.vars = "period",
       variable.name = "trader_type",
       value.name = "nb_country") %>%
  arrange(period) %>%
  mutate(trader_type = factor(trader_type, levels = order)) %>%
  ggplot(aes(x = period,
             y = nb_country,
             fill = trader_type,
             label = trader_type)) +
  geom_area(alpha = 0.75) +

  # End of chart labels
  annotate("text", x = 2022.4,
           y = data[data$period == 2022 &
                      data$cmd == 12, ]$tot_nb_nodes -
             data[data$period == 2022 &
                    data$cmd == 12, ]$nb_pure_imp / 2,
           label = paste(as.character(data[data$period == 2022 &
                                             data$cmd == 12, ]$nb_pure_imp),
                         "pure importers",
                         sep = " "),
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = pal[1]) +

  annotate("text", x = 2022.4,
           y = (data[data$period == 2022 &
                       data$cmd == 12, ]$nb_mixed +
                  data[data$period == 2022 &
                         data$cmd == 12, ]$nb_pure_exp) / 2,
           label = paste(as.character(data[data$period == 2022 &
                                             data$cmd == 12, ]$nb_mixed),
                         "mixed countries",
                         sep = " "),
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = pal[2]) +

  annotate("text", x = 2022.4,
           data[data$period == 2022 &
                  data$cmd == 12, ]$nb_pure_exp / 2,
           label = paste(as.character(data[data$period == 2022 &
                                             data$cmd == 12, ]$nb_pure_exp),
                         "pure exporters",
                         sep = " "),
           hjust = 0,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = pal[3]) +

  # Segments
  geom_segment(aes(x = 1996,
                   y = 0,
                   xend = 1996,
                   yend = data[data$period == 1996 &
                                 data$cmd == 12, ]$tot_nb_nodes),
               color = "black") +
  geom_point(aes(x = 1996,
                 y = data[data$period == 1996 &
                            data$cmd == 12, ]$tot_nb_nodes),
             color = "black") +
  annotate("text", x = 1996,
           y = data[data$period == 1996 &
                      data$cmd == 12, ]$tot_nb_nodes + 5,
           label = as.character(data[data$period == 1996 &
                                       data$cmd == 12, ]$tot_nb_nodes),
           hjust = 0.5,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = 2000,
                   y = 0,
                   xend = 2000,
                   yend = data[data$period == 2000 &
                                 data$cmd == 12, ]$tot_nb_nodes),
               color = "black") +
  geom_point(aes(x = 2000,
                 y = data[data$period == 2000 &
                            data$cmd == 12, ]$tot_nb_nodes),
             color = "black") +
  annotate("text", x = 2000,
           y = data[data$period == 2000 &
                      data$cmd == 12, ]$tot_nb_nodes + 5,
           label = as.character(data[data$period == 2000 &
                                       data$cmd == 12, ]$tot_nb_nodes),
           hjust = 0.5,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = 2005,
                   y = 0,
                   xend = 2005,
                   yend = data[data$period == 2005 &
                                 data$cmd == 12, ]$tot_nb_nodes),
               color = "black") +
  geom_point(aes(x = 2005,
                 y = data[data$period == 2005 &
                            data$cmd == 12, ]$tot_nb_nodes),
             color = "black") +
  annotate("text", x = 2005,
           y = data[data$period == 2005 &
                      data$cmd == 12, ]$tot_nb_nodes + 5,
           label = as.character(data[data$period == 2005 &
                                       data$cmd == 12, ]$tot_nb_nodes),
           hjust = 0.5,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = 2010,
                   y = 0,
                   xend = 2010,
                   yend = data[data$period == 2010 &
                                 data$cmd == 12, ]$tot_nb_nodes),
               color = "black") +
  geom_point(aes(x = 2010,
                 y = data[data$period == 2010 &
                            data$cmd == 12, ]$tot_nb_nodes),
             color = "black") +
  annotate("text", x = 2010,
           y = data[data$period == 2010 &
                      data$cmd == 12, ]$tot_nb_nodes + 5,
           label = as.character(data[data$period == 2010 &
                                       data$cmd == 12, ]$tot_nb_nodes),
           hjust = 0.5,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = 2015,
                   y = 0,
                   xend = 2015,
                   yend = data[data$period == 2015 &
                                 data$cmd == 12, ]$tot_nb_nodes),
               color = "black") +
  geom_point(aes(x = 2015,
                 y = data[data$period == 2015 &
                            data$cmd == 12, ]$tot_nb_nodes),
             color = "black") +
  annotate("text", x = 2015,
           y = data[data$period == 2015 &
                      data$cmd == 12, ]$tot_nb_nodes + 5,
           label = as.character(data[data$period == 2015 &
                                       data$cmd == 12, ]$tot_nb_nodes),
           hjust = 0.5,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = 2020,
                   y = 0,
                   xend = 2020,
                   yend = data[data$period == 2020 &
                                 data$cmd == 12, ]$tot_nb_nodes),
               color = "black") +
  geom_point(aes(x = 2020,
                 y = data[data$period == 2020 &
                            data$cmd == 12, ]$tot_nb_nodes),
             color = "black") +
  annotate("text", x = 2020,
           y = data[data$period == 2020 &
                      data$cmd == 12, ]$tot_nb_nodes + 5,
           label = as.character(data[data$period == 2020 &
                                       data$cmd == 12, ]$tot_nb_nodes),
           hjust = 0.5,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  geom_segment(aes(x = 2022,
                   y = 0,
                   xend = 2022,
                   yend = data[data$period == 2022 &
                                 data$cmd == 12, ]$tot_nb_nodes),
               color = "black") +
  geom_point(aes(x = 2022,
                 y = data[data$period == 2022 &
                            data$cmd == 12, ]$tot_nb_nodes),
             color = "black") +
  annotate("text", x = 2022,
           y = data[data$period == 2022 &
                      data$cmd == 12, ]$tot_nb_nodes + 5,
           label = as.character(data[data$period == 2022 &
                                       data$cmd == 12, ]$tot_nb_nodes),
           hjust = 0.5,
           size = 3,
           lineheight = .8,
           fontface = "bold",
           color = "black") +

  scale_fill_manual(values = pal) +
  scale_x_continuous(
    breaks = c(1996, 2000, 2005, 2010, 2015, 2020, 2022),
    labels = c("1996", "2000", "2005", "2010", "2015", "2020", "2022")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Year",
       y = "Number of trading countries") +
  coord_cartesian(clip = "off") +
  theme_ipsum(axis_title_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(40, 80, 20, 20))

ggsave(
  basename(snakemake@output[[1]]),
  plot = plot_compo,
  device = "png",
  path = dirname(snakemake@output[[1]]),
  scale = 1,
  width = 2480 * 1.1,
  height = 1240 * 1.5,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = "white"
)
