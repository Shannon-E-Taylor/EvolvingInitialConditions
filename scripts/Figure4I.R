
library('factoextra')

library('ggplot2')
library('data.table')
library('tidyr')
library('dplyr')
library('stringr')
library('patchwork')
library('broom')
library('latex2exp')
library('ggpmisc')

library('ggpubr')


# Import variables from config file 
source("config.R") 



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

# overview_measurements <- overview_measurements[overview_measurements$species %in% c('A. calliptera', 'R. chillingali'), ]


overview_measurements <- overview_measurements %>% 
  filter(
    species%in%c('A. calliptera', 'R. chillingali')& 
      somite_count > 1 &
      somite_count <=15
  )

overview_measurements$class <- '2-5ss'
overview_measurements[overview_measurements$somite_count > 5,]$class = '6-15ss'

overview_measurements$class <- factor(
  overview_measurements$class, levels = c('2-5ss', '6-15ss'))



plt_params <- list(
  scale_color_manual(values = cichlid_colours), 
  scale_fill_manual(values = cichlid_colours), 
  theme_bw(),
  theme(legend.position = 'none'), 
  scale_x_discrete(labels=c("AC","RC"))# , 
  # stat_compare_means(method = "t.test", vjust = 0.4)
)

plot_species_boxplot <- function(data, yvar, ylab) {
  
  yvar <- rlang::ensym(yvar)
  
  fc <- data %>%
    group_by(class, species) %>%
    summarise(mean = mean(!!yvar, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = species, values_from = mean) %>%
    mutate(
      label = sprintf("%.2f-fold",
                      `R. chillingali` / `A. calliptera`),
      x = 1.5,
      y = min(data[[rlang::as_name(yvar)]], na.rm = TRUE) * 1.05
    )
  
  ggplot(data,
         aes(species, !!yvar,
             colour = species,
             fill = species)) +
    geom_boxplot(alpha = .3, outlier.shape = NA) +
    geom_jitter(width = .1, size = 1) +
    stat_summary(fun = mean, geom = "point",
                 colour = "red", shape = 3) +
    stat_compare_means(method = 't.test', size = 3)  + 
    geom_text(data = fc,
              aes(x = x, y = y, label = label),
              size = 3, 
              inherit.aes = FALSE) +
    facet_wrap(~class) +
    plt_params +
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.1))
    ) + 
    labs(y = ylab, x = NULL)
}


psm <- plot_species_boxplot(overview_measurements, 'PSM_length', 'PSM length (\u00B5m)')

el <- plot_species_boxplot(overview_measurements, 'embryo_length', 'Embryo length (\u00B5m)')
psm + el


ggsave('../graphs/Fig4I.png', width = 6, height = 3, dpi =600)



