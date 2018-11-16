library(ggplot2)

ns_ws_color_values <- c("dodgerblue", "darkorange2")
volcano_colors <- scale_color_manual(values = c("grey50", "seagreen3"), guide = FALSE)
ns_ws_colors <- rep(ns_ws_color_values, 11)
# from brewer.pal(name = 'PuBu', n = 9)
heatmap_colors <- c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF",
                    "#3690C0", "#0570B0", "#045A8D", "#023858")
