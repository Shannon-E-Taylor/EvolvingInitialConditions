


library('ggplot2')
library('data.table')
library('tidyr')
library('dplyr')
library('patchwork')
library('ggpubr')

source('config.R')





data = read.csv('../data/all_overviews.csv')


PSM_len_synonyms <- c("psm_length", "PSM_lnegth", "psm", "PSM_lenght", 
                      "PSM_length                "
)

embryo_len_synonums <- c("embyro_length", "embryo_length        ",  
                         "embryo_length                 ", 
                         "embryo_length "
                         
)


data <- data %>% 
  mutate(
    Text = if_else(
      Text %in% PSM_len_synonyms,
      "PSM_length",
      Text
    )
  )%>% 
  mutate(
    Text = if_else(
      Text %in% embryo_len_synonums,
      "embryo_length",
      Text
    )
  ) 




# Extract lengths
lengths <- data[data$Text %in% c(
  'embryo_length', 'PSM_length', 'PSM_width', 'somite_length',
  'species', 'somite_stage', 
  'somite_width', 'paraxial-mesoderm', 'head', 'hw'),
  c('image_id', 'Text', 'length', 'species', 'somite_count')]

overview_measurements <- pivot_wider(lengths, 
                                     names_from = 'Text', 
                                     values_from = 'length', values_fn = mean)

overview_measurements$species <- factor(overview_measurements$species, 
                                        levels = c('A. calliptera', 'R. chillingali', 'M. zebra'))

overview_measurements <- overview_measurements[overview_measurements$species %in% c('A. calliptera', 'R. chillingali'), ]

overview_measurements <- overview_measurements[overview_measurements$somite_count > 0, ]

plt_params <- list(
  scale_color_manual(values = cichlid_colours), 
  scale_fill_manual(values = cichlid_colours), 
  theme_bw(),
  theme(legend.position = 'none'), 
  scale_x_discrete(labels=c("AC","RC"))# , 
  # stat_compare_means(method = "t.test", vjust = 0.4)
)





embryo_len <- ggplot(overview_measurements, 
                  aes(x = somite_count, y = embryo_length, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylim(0, 2600) + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ylab('Embryo length (\u00B5m)') + 
  xlab('Somite count') + 
  ggtitle('B. Embryo length') +   theme(legend.position = 'none') + 
  stat_poly_line(se = FALSE)



psm_len <- ggplot(overview_measurements, 
                  aes(x = somite_count, y = PSM_length, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylim(0, 900) + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ylab('PSM length (\u00B5m)') + 
  xlab('Somite count') + 
  ggtitle('D. PSM length') + 
  stat_poly_line(
    data = overview_measurements[
      overview_measurements$somite_count < 15,], 
    se = FALSE) +   theme(legend.position = 'none') + 
  stat_poly_line(data = overview_measurements[overview_measurements$somite_count > 15,], se = FALSE)


psm_len.prop <- ggplot(overview_measurements, 
                  aes(x = somite_count, y = PSM_length/embryo_length*100, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylim(0, 60) + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ylab('PSM length (%EL)') + 
  xlab('Somite count') +   theme(legend.position = 'none') + 
  ggtitle('E. Relative PSM length') + 
  stat_poly_line(
    se = FALSE) 

head_len.prop <- ggplot(overview_measurements, 
                       aes(x = somite_count, y = head/embryo_length*100, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylim(0, 60) + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ylab('Head length (%EL)') + 
  xlab('Somite count') + 
  ggtitle('F. Relative head length') + 
  theme(legend.position = 'none') + 
  stat_poly_line(
    se = FALSE) 




overview_measurements$head_prop = (overview_measurements$head / overview_measurements$embryo_length) * 100
overview_measurements$psm_prop = 100 - (overview_measurements$PSM_length / overview_measurements$embryo_length) * 100


overview_measurements_m <- overview_measurements |>
  select(somite_count, species, head_prop, psm_prop) |>
  tidyr::pivot_longer(
    cols = c(head_prop, psm_prop),
    names_to = "variable",
    values_to = "value"
  )


proportions <- ggplot(overview_measurements_m, 
                      aes(x = somite_count, y = value, color = species, 
                          shape = variable, # group = interaction(variable), 
                          )
                      ) +  
  geom_smooth(method = 'loess', se = FALSE, # color = 'black', 
              linewidth = 1) + 
  geom_point(size = 1.5, stroke=NA, alpha = 1) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) + 
  # scale_x_continuous(limits = c(0, 100), expand = c(0, 2)) + 
  scale_colour_manual(values = cichlid_colours) +  # , labels = c('AC', 'RC', 'MZ'))+
  ylab("AP position (%)") + 
  xlab("Somite count") + 
  coord_flip() + 
  scale_x_reverse(limits= c(42, 0), expand = c(0,0)) + 
  scale_shape_manual(values = c("head_prop" = 16, "psm_prop" = 17), 
                     labels = c("head_prop" = "Head boundary", "psm_prop" = "PSM boundary")) +
  theme(legend.position = 'bottom') + 
  ggtitle('G. Embryo proportions over time')  + 
  theme(legend.margin=margin(-10, 0, 0, 0)) + 
  labs(shape = 'Boundary') + 
  guides(
    color = guide_legend(nrow = 1),   # Put species legend entries in one row
    shape = guide_legend(nrow = 1)    # Put shape legend entries in one row
  ) +  
  theme(legend.position = 'bottom', 
        legend.box = "vertical")  # Stack legends vertically


design <- '
ABBB
CDDD
EEFF
GGGG
'


w = 0.6
lw = 8

plt_params_2 <- list(
  scale_color_manual(values = cichlid_colours), 
  scale_fill_manual(values = cichlid_colours), 
  theme_bw(),
  theme(legend.position = 'bottom') 
  # stat_compare_means(method = "t.test", vjust = 0.4)
)



el <- ggplot(overview_measurements[overview_measurements$somite_count > 1 & 
                                     overview_measurements$somite_count <=5,], 
             aes(x = somite_count, y = embryo_length, 
                 color = species, fill= species, group = species)) + 
  # geom_boxplot(alpha = 0.3, outliers = FALSE, notch = FALSE) +
  # geom_jitter(width = 0.1, height = 0, size = 1) +
  geom_point(position = position_dodge(width = w), 
             alpha = 0.3) + 
  plt_params_2 + 
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 2, 
    # alpha = 0.7,
    position = position_dodge(width = w), 
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3,
    # shape = 0,
    # color = 'red', 
    position = position_dodge(width = w), 
    show.legend=FALSE
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.3,
    linewidth = 1, alpha = 0.7, 
    position = position_dodge(width = w), 
  ) +
  # stat_compare_means(method = 't.test') + 
  scale_x_continuous(
    breaks = 2:5,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme(
    panel.grid.major.x = element_line(colour = "grey95",linewidth = lw),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_line(colour = "white"),
    panel.grid.minor.y = element_line(colour = "white"),
    legend.position = 'none'
  ) + 
  ylab('Embryo length (\u00B5m)') + 
  xlab('Somite count') + 
  ggtitle('A. Embryo length (early)')



w = 0.6
lw = 8

plt_params_2 <- list(
  scale_color_manual(values = cichlid_colours), 
  scale_fill_manual(values = cichlid_colours), 
  theme_bw(),
  theme(legend.position = 'bottom') 
  # stat_compare_means(method = "t.test", vjust = 0.4)
)

psm <- ggplot(overview_measurements[overview_measurements$somite_count > 1 & 
                                      overview_measurements$somite_count <=5,], 
              aes(x = somite_count, y = PSM_length, 
                  color = species, fill= species, group = species)) + 
  # geom_boxplot(alpha = 0.3, outliers = FALSE, notch = FALSE) +
  # geom_jitter(width = 0.1, height = 0, size = 1) +
  geom_point(position = position_dodge(width = 0.6), alpha = 0.3) + 
  plt_params_2 + 
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 2, 
    # alpha = 0.7,
    position = position_dodge(width = w), 
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3,
    # shape = 0,
    # color = 'red', 
    position = position_dodge(width = w), 
    show.legend=FALSE
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.3,
    linewidth = 1, alpha = 0.7, 
    position = position_dodge(width = w), 
  ) +
  scale_x_continuous(
    breaks = 2:5,
    minor_breaks = NULL,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme(
    panel.grid.major.x = element_line(colour = "grey95",linewidth = lw),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_line(colour = "white"),
    panel.grid.minor.y = element_line(colour = "white"),
    legend.position = 'none' 
  ) + 
  ylab('PSM length (\u00B5m)') + 
  xlab('Somite count') + 
  ggtitle('C. PSM length (early)')


el + embryo_len +  psm + psm_len + psm_len.prop + head_len.prop + proportions + 
  plot_layout(design = design, axes="collect_y")

ggsave('../graphs/FigS3.png', dpi = 600, width = 8, height = 10)
