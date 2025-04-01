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
library("ggh4x")
library("ggrepel")
library("ggpubr")
library("hrbrthemes")
library("scales")
library("patchwork")
library("dplyr")
library("grid")

data <- read.csv(
  file = "/Users/valentinmathieu/Desktop/wd/trade_discrepancies/results/processed_data/network_analysis/output/network_metrics.csv",
  header = TRUE,
  sep = ";")
head(data)

head(data[data$omission == "" & data$product == 12 & data$period == 1999, ])

df_exp <- data[(data$product == 12 &
                  data$trader_type == "exporter" &
                  data$omission == ""),
               c("period",
                 "tot_nb_nodes",
                 "nb_nodes",
                 "traded_value")]
df_exp <- df_exp[order(df_exp$period), ]
df_exp$traded_value <- df_exp$traded_value / 1000000000
df_exp$node_share <- df_exp$nb_nodes / df_exp$tot_nb_nodes * 100
df_imp <- data[(data$product == 12 &
                  data$trader_type == "importer" &
                  data$omission == ""),
               c("period",
                 "tot_nb_nodes",
                 "nb_nodes",
                 "traded_value")]
df_imp <- df_imp[order(df_imp$period), ]
df_imp$traded_value <- df_imp$traded_value / 1000000000
df_imp$node_share <- df_imp$nb_nodes / df_imp$tot_nb_nodes * 100

coeff <- 0.115 # Coefficient to scale second y-axis.

plot_size <- ggplot(df_exp,
                    aes(x = period,
                        group = 1)) +
  geom_point(data = df_exp,
             aes(y = nb_nodes),
             color = "red",
             shape = 16,
             size = 2.2) +
  geom_point(data = df_imp,
             aes(y = nb_nodes),
             color = "cyan4",
             shape = 15,
             size = 2.2) +
  geom_line(aes(y = nb_nodes), color = "red") +
  geom_line(data = df_imp,
            aes(x = period,
                y = nb_nodes),
            linetype = "longdash",
            color = "cyan4") +
  geom_point(data = df_exp,
             aes(y = tot_nb_nodes),
             color = "black",
             alpha = 1,
             shape = 18,
             size = 2.2) +
  geom_line(aes(y = tot_nb_nodes),
            color = "black",
            alpha = 0.5) +
  geom_area(data = df_exp,
            aes(x = period,
                y = traded_value / coeff),
            fill = "red",
            color = "red",
            alpha = 0.4) +
  geom_area(data = df_imp,
            aes(x = period,
                y = traded_value / coeff),
            fill = "cyan4",
            color = "cyan4",
            linetype = "longdash",
            alpha = 0.4) +
  scale_y_continuous(
    # Features of the first axis
    name = "Nb. of nodes",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~. * coeff,
                        name = "Total amount of traded value (Billion US$)"),
  ) +
  labs(x = "Year") +
  theme_ipsum(axis_title_size = 11)

ggsave(
  "size.png",
  plot = plot_size,
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

df <- data[(data$product == 12 &
              data$omission == ""),
           c("period",
             "trader_type",
             "tot_nb_nodes",
             "nb_pure_exp",
             "nb_pure_imp",
             "nb_mixed")]
df <- df[order(df$period), ]
df$node_share <- df$nb_nodes / df$tot_nb_nodes * 100
df_exp <- df[df$trader_type == "exporter",
             c("period", "trader_type", "node_share")]
df_imp <- df[df$trader_type == "importer",
             c("period", "trader_type", "node_share")]
df_share <- cbind(exp = df_exp, imp = df_imp)
df_share <- df_share[, c("exp.period", "exp.node_share", "imp.node_share")]
names(df_share)[names(df_share) == "exp.period"] <- "period"

plot_nb_nodes <- ggplot(df,
                        aes(x = period,
                            group = trader_type,
                            color = trader_type,
                            shape = trader_type)) +
  geom_point(aes(y = nb_nodes)) +
  geom_line(aes(y = nb_nodes)) +
  geom_line(aes(y = tot_nb_nodes), color = "black") +
  ylim(110, 215) +
  # Colors for the lines
  scale_color_manual(values = c("#C32E5A", "#3D85F7", "black")) +
  labs(x = "Year",
       y = "Number of trading countries",
       color = "Type of trading country",
       shape = "Type of trading country") +
  coord_cartesian(xlim = c(1996, 2022)) +
  theme_ipsum(axis_title_size = 11) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

plot_share_node <- ggplot(df_share,
                          aes(x = period)) +
  geom_line(aes(y = exp.node_share, color = "exporter")) +
  geom_line(aes(y = imp.node_share, color = "importer")) +
  stat_difference(aes(ymin = exp.node_share,
                      ymax = imp.node_share),
                  alpha = 0.3) +
  ylim(65, 100) +
  # Colors for the lines
  scale_color_manual(values = c("#C32E5A", "#3D85F7")) +
  # Colors for the fill. They are lighter versions of the line colors.
  # The third one is required because lines are sometimes equal
  scale_fill_manual(
    values = c(
      colorspace::lighten("#3D85F7"),
      colorspace::lighten("#C32E5A"),
      "grey60"
    ),
    labels = c("More importers", "More exporters", "same")
  ) +
  labs(x = "Year",
       y = "Share of trading countries in (%)",
       color = "Type of trading country") +
  coord_cartesian(xlim = c(1996, 2022)) +
  theme_ipsum(axis_title_size = 11)

# plot_compo <- plot_nb_nodes / plot_share_node

# grid.newpage()
# plot_compo <- grid.draw(rbind(ggplotGrob(plot_nb_nodes),
#                               ggplotGrob(plot_share_node),
#                               size = "last"))

plot_compo <- plot_nb_nodes + plot_spacer() + 
  plot_share_node + plot_layout(ncol = 1, heights = c(10, -4, 10))

ggsave(
  "compo.png",
  plot = plot_compo,
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

