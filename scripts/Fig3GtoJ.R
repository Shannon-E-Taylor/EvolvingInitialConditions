
library('ggplot2')
library('data.table')
library('tidyr')
library('dplyr')
library('patchwork')
library('ggpubr')

source('config.R')
# contrast colour #7A9E9F

# Read in, and proccess, the data
volumes <- read.csv('../data/volumes.csv')
volumes$somite_count <- as.integer(volumes$somite_count)

# Remove data lacking somite count and incorrect species 
volumes <- volumes[volumes$somite_count > 5 ,]#
volumes <- volumes[volumes$species %in% c('A. calliptera', 'R. chillingali'), ]

# somite volume 0 means I didn't measure it 
# (because somite is not complete in image)
volumes[volumes$mean_volume_of_somites == 0, ]$mean_volume_of_somites = NA

# and get rid of mis-annotated data
to_remove = c('Merian_xx', # doesn't actually have psm - mislabelled 
              'Rickshaw_1', 'Rickshaw_2') # staging of these is wrong. 
volumes = volumes[volumes$incomplete != TRUE,]
volumes = volumes[!volumes$embryo_id %in% c(to_remove), ]

volumes <- assign_class(
  volumes, 
  AC_midseg_threshold = AC_midseg_threshold, 
  RC_midseg_threshold = RC_midseg_threshold, 
  AC_lateseg_threshold = AC_lateseg_threshold, 
  RC_lateseg_threshold = RC_lateseg_threshold
)

volumes$class <- factor(
  volumes$class,
  levels = c("onset", setdiff(levels(volumes$class), "onset"))
)
volumes$class[volumes$somite_count <= 5] <- "onset"

# Read in, and process, the cell counts data 
cell_counts <- read.csv('../data/cell_numbers.csv')
cell_counts$somite_count <- as.integer(cell_counts$somite_count)

# Remove data lacking somite count and incorrect species 
cell_counts <- cell_counts[cell_counts$somite_count > 5 ,]#
cell_counts <- cell_counts[cell_counts$species %in% c('A. calliptera', 'R. chillingali'), ]

somite_const = 15


mean_chilli_vol = mean(volumes[volumes['somite_count'] < somite_const & volumes['species'] == 'R. chillingali', ]$psm_volume, na.rm = TRUE)
mean_calli_vol = mean(volumes[volumes['somite_count'] < somite_const & volumes['species'] == 'A. calliptera', ]$psm_volume, na.rm = TRUE)

mean_chilli_count = mean(cell_counts[cell_counts$somite_count < somite_const & cell_counts$species == 'R. chillingali', ]$cells_in_psm, na.rm = TRUE)
mean_calli_count = mean(cell_counts[cell_counts$somite_count < somite_const & cell_counts$species == 'A. calliptera', ]$cells_in_psm, na.rm = TRUE)



plt_params <- list(
  scale_color_manual(values = cichlid_colours), 
  scale_fill_manual(values = cichlid_colours), 
  theme_bw(),
  theme(legend.position = 'none'), 
  scale_x_discrete(labels=c("AC","RC"))# , 
  # stat_compare_means(method = "t.test", vjust = 0.4)
)


# Helper to calculate fold changes
fold_change_labels <- function(data, value, group = species, facet = class,
                               ref = "A. calliptera",
                               test = "R. chillingali",
                               y_mult = 1.05) {
  
  value <- rlang::ensym(value)
  group <- rlang::ensym(group)
  facet <- rlang::ensym(facet)
  
  data %>%
    group_by(!!facet, !!group) %>%
    summarise(mean = mean(!!value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!group, values_from = mean) %>%
    mutate(
      fold = .data[[test]] / .data[[ref]],
      x = 1.5,
      y = max(data[[rlang::as_name(value)]], na.rm = TRUE) * y_mult,
      label = sprintf("%.2f-fold", fold)
    )
}

pos_const = 10


stat.lm <- lm(
  psm_volume~somite_count+species, 
  data = volumes[volumes$class =='early', ])
anova_res <- anova(stat.lm)
eta_res <- eta_squared(stat.lm)

p_species <- anova_res["species", "Pr(>F)"]
eta_species <- eta_res$Eta2[eta_res$Parameter == "species"]

label <- paste0(
  "p = ", format.pval(p_species, digits = 2),
  "\nη² = ", round(eta_species, 2)
)


volumes %>% 
  # filter(class=='early') %>%
  # group_by(date, species)  %>%
  # exclude clutches where we have less than two entries, 
  # and where there is less than one hour difference between ages 
  # as that will make the regression unreliable 
  filter(!is.na(psm_volume)) %>%
  count(species)


vol <- ggplot(volumes[volumes$class =='early', ], 
              aes(x = species, y = psm_volume/1000000, 
                  color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE, notch = FALSE) + 
  geom_jitter(width = 0.1, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 9/pos_const, size = 3, 
           label = eval(paste(c(round(mean_chilli_vol / mean_calli_vol, 2)), 'fold\nchange'))) + 
  ylab('PSM volume (10\u2076 \u00B5m\u00b3)') +
  scale_y_continuous(limits = c(0, 9)) + 
  ggtitle('D. PSM volume') +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = label,
    vjust = 1.5,
    size = 3, 
  )


stat.lm <- lm(
  cells_in_psm~somite_count+species, 
  data = cell_counts[cell_counts$somite_count < 15, ])
anova_res <- anova(stat.lm)
eta_res <- eta_squared(stat.lm)

p_species <- anova_res["species", "Pr(>F)"]
eta_species <- eta_res$Eta2[eta_res$Parameter == "species"]

label <- paste0(
  "p = ", format.pval(p_species, digits = 2),
  "\nη² = ", round(eta_species, 2)
)

cell_counts %>% 
  filter(somite_count < 15) %>%
  # group_by(date, species)  %>%
  # exclude clutches where we have less than two entries, 
  # and where there is less than one hour difference between ages 
  # as that will make the regression unreliable 
  filter(!is.na(cells_in_psm)) %>%
  count(species)


cells <- ggplot(cell_counts[cell_counts$somite_count < 15, ], 
                aes(x = species, y = cells_in_psm, color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE, notch = FALSE) + 
  geom_jitter(width = 0.1, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 10000/pos_const, size = 3, 
           label = eval(paste(c(round(mean_chilli_count / mean_calli_count, 2)), 'fold\nchange'))) + 
  ylab('Cells in PSM') + 
  ylim(0, 10000) + 
  ggtitle('G. Cells in PSM')+
  annotate(
    "text",
    x = 1,
    y = Inf,
    label = label,
    vjust = 1.5,
    size = 3, 
  )




a_guess = 100000
r_guess = 0

get_exponential_coefficents <- function(df, variable){
  
  exponential_model = nls(
    as.formula(paste(variable, "~ a*exp(r*somite_count)")), 
    data = df, 
    start = list(a = a_guess, r = r_guess) 
  )
  # extract parameters. 
  exponential_coef = coef(exponential_model)
  a = exponential_coef[1]
  r = exponential_coef[2]
  
  return(c(a, r))
  
}

out <- get_exponential_coefficents(
  volumes[volumes$species == 'R. chillingali' & 
            volumes$class == 'early', ], 
  variable = 'psm_volume'
)
a_RC_early_psm = out[[1]]
r_RC_early_psm = out[[2]]

out = get_exponential_coefficents(
  volumes[volumes$species == 'R. chillingali' & 
            volumes$class != 'early', ], 
  variable = 'psm_volume'
)
a_RC_midlate_psm = out[[1]]
r_RC_midlate_psm = out[[2]]


out = get_exponential_coefficents(
  volumes[volumes$species == 'A. calliptera' & 
            volumes$class == 'early', ], 
  variable = 'psm_volume'
)
a_AC_early_psm = out[[1]]
r_AC_early_psm = out[[2]]

out = get_exponential_coefficents(
  volumes[volumes$species == 'A. calliptera' & 
            volumes$class != 'early', ], 
  variable = 'psm_volume'
)


a_AC_midlate_psm = out[[1]]
r_AC_midlate_psm = out[[2]]

exp_function <- function(x, a, r) {
  a * exp(r * x) / 1000000
}




volumes_lm <- lm(
  log(psm_volume) ~ somite_count + species,
  data = volumes[volumes$somite_count > 15, ]
)

anova_res <- anova(volumes_lm)
eta_res <- eta_squared(volumes_lm)

p_somite <- anova_res["somite_count", "Pr(>F)"]
p_species <- anova_res["species", "Pr(>F)"]

eta_somite <- eta_res$Eta2[eta_res$Parameter == "somite_count"]
eta_species <- eta_res$Eta2[eta_res$Parameter == "species"]

label <- paste0(
  "Somite count:  p = ", format.pval(p_somite, digits = 2),
  ", η² = ", round(eta_somite, 2),
  "\nSpecies:  p = ", format.pval(p_species, digits = 2),
  ", η² = ", round(eta_species, 2), ''
)





psm_vol <- ggplot(volumes, 
                  aes(x = somite_count, y = psm_volume/1000000, 
                      color = species), 
) + 
  geom_point(size = 1) + 
  theme_bw() + 
  scale_color_manual(values = cichlid_colours) +
  geom_function(fun = exp_function,
                args = list(a = a_RC_midlate_psm, r = r_RC_midlate_psm),
                color = cichlid_colours[[2]], 
                linewidth = 1,
                xlim = c(17, 40)) +
  geom_function(fun = exp_function,
                args = list(a = a_AC_midlate_psm, r = r_AC_midlate_psm),
                color = cichlid_colours[[1]], 
                linewidth = 1,
                xlim = c(17, 32)) +
  # stat_poly_line(data= volumes[volumes$somite_count < somite_const, ], se = FALSE) + 
  ylab('PSM volume (10\u2076 \u00B5m\u00b3)') +
  xlab('Somite count') + 
  ggtitle('E. PSM volume') + 
  scale_y_continuous(limits = c(0, 9)) + 
  xlim(0, 40)+ 
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label,
    hjust = 1.1,
    vjust = 1.5, 
    size = 3
  ) +
  theme(legend.position = 'none')

out = get_exponential_coefficents(
  volumes[#volumes$species == 'R. chillingali' & 
    volumes$class == 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_RC_som_early = out[[1]]
r_RC_som_early = out[[2]]

out = get_exponential_coefficents(
  volumes[#volumes$species == 'A. calliptera' &
    volumes$class == 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_AC_som_early = out[[1]]
r_AC_som_early = out[[2]]

out = get_exponential_coefficents(
  volumes[#volumes$species == 'R. chillingali' & 
    volumes$class != 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_RC_som_midlate = out[[1]]
r_RC_som_midlate = out[[2]]

out = get_exponential_coefficents(
  volumes[#volumes$species == 'A. calliptera' &
    volumes$class != 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_AC_som_midlate = out[[1]]
r_AC_som_midlate = out[[2]]

som_vol <- ggplot(volumes, 
                  aes(x = somite_count, 
                      y = mean_volume_of_somites/1000000, 
                      color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  scale_color_manual(values = cichlid_colours) +
  geom_function(fun = exp_function, 
                args = list(a = a_AC_som_midlate, r = r_AC_som_midlate), 
                color = cichlid_colours[[1]],
                linewidth = 1,
                xlim = c(16, 30)
  ) + 
  geom_function(fun = exp_function, 
                args = list(a = a_RC_som_midlate, r = r_RC_som_midlate), 
                color = cichlid_colours[[2]],
                linewidth = 1,
                xlim = c(16, 40)) + 
  ylab('Nascent somite \nvolume (10\u2076 \u00B5m\u00b3)') + 
  xlab('Somite count') + 
  xlim(0, 40) + 
  ylim(0, .25)+ 
  ggtitle('F. Somite volume') + 
  theme(legend.position = 'none')


ymax = 1000





volumes_lm <- lm(
  cells_in_psm ~ somite_count * species,
  data = cell_counts[cell_counts$somite_count > 15, ]
)

anova_res <- anova(volumes_lm)
eta_res <- eta_squared(volumes_lm)

p_somite <- anova_res["somite_count", "Pr(>F)"]
p_species <- anova_res["species", "Pr(>F)"]

eta_somite <- eta_res$Eta2[eta_res$Parameter == "somite_count"]
eta_species <- eta_res$Eta2[eta_res$Parameter == "species"]

label <- paste0(
  "Somite count:  p = ", format.pval(p_somite, digits = 2),
  ", η² = ", round(eta_somite, 2),
  "\nSpecies:  p = ", format.pval(p_species, digits = 2),
  ", η² = ", round(eta_species, 2), ' '
)




psm_cells <- ggplot(cell_counts, 
                    aes(x = somite_count, y = cells_in_psm, color = species), 
) + 
  geom_point(size = 1) + 
  theme_bw() + 
  scale_color_manual(values = cichlid_colours) +
  ylab('Cells in PSM') + 
  xlab('Somite count') + 
  ggtitle('H. Cells in PSM') + 
  xlim(0, 40) + 
  ylim(0, 10000) + 
  theme(legend.position = 'bottom') +   annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label,
    hjust = 1.1,
    vjust = 1.5, 
    size = 3
  ) +
  stat_poly_line(data = cell_counts[cell_counts$somite_count > somite_const,], se = FALSE)


design = '
ABBCDD
'



vol + ggtitle('G. PSM volume') + psm_vol +  ggtitle('H. PSM volume') + 
  cells + ggtitle('I. Cells in PSM') + psm_cells + ggtitle('J. Cells in PSM') +
  plot_layout(design = design, # guides = 'collect', 
              axes = 'collect_y') & 
  theme(legend.position='none')

ggsave('../graphs/Fig3GtoJ.png', width = 10, height = 3, dpi =600)


