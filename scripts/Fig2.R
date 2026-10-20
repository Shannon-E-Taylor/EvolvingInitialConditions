library('ggplot2')
library('ggpmisc')
library('effectsize')
library('patchwork')
library('magrittr')
library('dplyr')
library('broom')

data <- read.csv('../data/somite_counts.csv')
cichlid_colours = c('#C3BE46FF', '#4D4D4DFF',  '#00CCFF') # calliptera, Ramphos, zebras


psm_data = read.csv('../data/computed_psm_dimensions.csv')

psm_data$somite_count = as.integer(psm_data$somite_count)

psm_data = psm_data[
  (psm_data$species %in% c('A. calliptera', 'R. chillingali')) & (psm_data$somite_count > 5), 
]


#####
# stats 
#####

som_lm <- lm(som ~ somite_count + species, data = psm_data 
)

anova_res <- anova(som_lm)
eta_res <- eta_squared(som_lm)


p_somite <- anova_res["somite_count", "Pr(>F)"]
p_species <- anova_res["species", "Pr(>F)"]

eta_somite <- eta_res$Eta2[eta_res$Parameter == "somite_count"]
eta_species <- eta_res$Eta2[eta_res$Parameter == "species"]

label <- paste0(
  "Somite count:  p = ", format.pval(p_somite, digits = 2),
  ", η² = ", round(eta_somite, 2),
  "\nSpecies:  p = ", format.pval(p_species, digits = 2),
  ",   η² = ", round(eta_species, 2), ' '
)




psm_data %>% 
  # group_by(date, species)  %>%
  # exclude clutches where we have less than two entries, 
  # and where there is less than one hour difference between ages 
  # as that will make the regression unreliable 
  filter(!is.na(som)) %>%
  count(species)




som_len <- psm_data %>% 
  filter(species%in%c('A. calliptera', 'R. chillingali')) %>% 
  ggplot(
              aes(x = somite_count, y = som, 
                  colour = species)
) +
  theme_bw() + 
  xlim(0, 40) + 
  xlab('Somite Stage') + 
  scale_colour_manual(values = cichlid_colours) + 
  scale_y_continuous(expand = expansion(mult = 0.1)) + 
  ylab('Somite length (\u00B5m)') +
    stat_poly_line(se = FALSE) + 
    # stat_poly_eq(use_label(c('eq', 'R2')),
    #              label.x = "right",
    #              label.y = c(0.95, 0.85))  +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label,
    hjust = 1.1,
    vjust = 1.5
  ) +
  geom_point(size = 1) + theme(legend.position = 'bottom') + 
  ggtitle('C. Nascent somite length')




#

#########
# Somite volumnes 

source('config.R')
# contrast colour #7A9E9F

# Read in, and proccess, the data
volumes <- read.csv('../data/volumes.csv')
volumes$somite_count <- as.integer(volumes$somite_count)

# Remove data lacking somite count and other species 
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

exp_function <- function(x, a, r) {
  a * exp(r * x) / 1000000
}

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



volumes_lm <- lm(
  log(mean_volume_of_somites) ~ somite_count + species,
  data = volumes[volumes$somite_count > 15, ]
)

anova_res <- anova(volumes_lm)
eta_res <- eta_squared(volumes_lm)

p_somite <- anova_res["somite_count", "Pr(>F)"]
p_species <- anova_res["species", "Pr(>F)"]

eta_somite <- eta_res$Eta2[eta_res$Parameter == "somite_count"]
eta_species <- eta_res$Eta2[eta_res$Parameter == "species"]

label <- paste0(
  "Somite count: p = ", format.pval(p_somite, digits = 2),
  ",  η² = ", round(eta_somite, 2),
  "\nSpecies: p = ", format.pval(p_species, digits = 2),
  ",       η² = ", round(eta_species, 2), ' '
)


volumes %>% 
  # group_by(date, species)  %>%
  # exclude clutches where we have less than two entries, 
  # and where there is less than one hour difference between ages 
  # as that will make the regression unreliable 
  filter(!is.na(mean_volume_of_somites)) %>%
  count(species)


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
  annotate(
                  "text",
                  x = Inf,
                  y = Inf,
                  label = label,
                  hjust = 1.1,
                  vjust = 1.5
                ) +
  ylab('Nascent somite volume (10\u2076 \u00B5m\u00b3)') + 
  xlim(0, 40) + 
  xlab('Somite Stage') + 
  ylim(0, .25)+ 
  ggtitle('D. Somite volume') + 
  theme(legend.position = 'bottom')





### CLOCK RATE 


library(broom)

data <- read.csv('../data/somite_counts.csv')


datesfemales_to_avoid = c("13/03/2023", "24/03/2024")

data <- data %>% 
  filter(
  !date%in%datesfemales_to_avoid & 
    somites > 0
  )

# data <- data[!data$date %in% datesfemales_to_avoid, ]

data_sub <- data %>%
  # exclude A. calliptera older than 70 hours from data
  # as they have finished segmentation
  filter(
    (species != "Astatotilapia calliptera") | (min_age < 70)
  ) %>% 
  # exclude clutches where we have less than two entries, 
  # and where there is less than one hour difference between ages 
  # as that will make the regression unreliable 
  group_by(date, species) %>%
  filter(n() > 2,  (max(min_age) - min(min_age)) >= 2) %>% 
  ungroup() 

data_sub$age_mins = data_sub$min_age * 60





df_lm <- data_sub %>% 
  group_by(date, species)  %>%
  do({
    model <- lm(somites ~ age_mins, data = .)
    tidy_model <- tidy(model)
    tidy_model$n <- nrow(.)  # Add number of data points used
    tidy_model
  }) 

data_sub %>% 
  distinct(date, species) %>%
  count(species)

keep_dates = unique(df_lm$date)


data_sub %>% 
  count(species)


df_mean <- df_lm %>%
  filter(term != '(Intercept)') %>% 
  group_by(species) %>%
  summarise(mean_value = mean(estimate, na.rm = TRUE), 
            median_value = median(estimate, na.rm = TRUE))

label_df <- df_mean %>%
  mutate(
    label = paste0(
      "Mean: ", round(1 / mean_value, 1),
      "\nMedian: ", round(1 / median_value, 1)
    ),
    y_pos = Inf
  )


df_mean$slope = 1/df_mean$median_value

hjust = 0.4

rate <- ggplot(df_lm[df_lm$term!='(Intercept)', ], 
               aes(x = species, y = 1 / estimate, 
                   color = species, 
                   shape = species,
                   fill = species#, size = n/2
               )) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE, width = 0.5) + 
  geom_jitter(height=0, width = 0.1) + 
  scale_colour_manual(values = cichlid_colours) + 
  scale_fill_manual(values = cichlid_colours) + 
  theme_bw() + 
  # ylab('Somites produced per hour') + 
  ylab('Rate of somite production (mins)') + 
  ylim(0, 150) + 
  xlab('') + 
  theme(legend.position = 'none') + 
  ggtitle('A. Indistinguishable rate') + 
  # stat_summary(fun='mean', size = 1.5, position=position_nudge(x = -hjust, y = 0)) + 
  scale_x_discrete(labels=c('AC', 'RC')) + 
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "t.test", 
                     label.y = 140, 
                     label.x = 1.5)+ 
  stat_summary(fun = mean, geom = "point",
               colour = "red", shape = 3) + 
  geom_text(
    data = label_df,
    aes(
      x = species,
      y = -Inf,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = -1,
    color = "black",
    size = 3
  )

# compare our data to Marconi et al. 2023 
# rate + geom_point(
#   data = data.frame(species = c("Rhamphochromis sp. 'chilingali'", 
#                                 "Astatotilapia calliptera"), 
#                     Value = c(68, 74)),
#   aes(x = as.numeric(factor(species)) + 0.1, y = Value),
#   hjust = 1.5,
#   color = "blue",
#   size = 3
# )



slope_AC = df_mean[df_mean$species=='Astatotilapia calliptera',]$slope / 60
slope_RC = df_mean[df_mean$species=="Rhamphochromis sp. 'chilingali'",]$slope / 60

data$slope = slope_AC
data[data$species=="Rhamphochromis sp. 'chilingali'", ]$slope = slope_RC


lm_AC <- lm(min_age - slope*somites ~ 1, 
            data = subset(data, species=='Astatotilapia calliptera'))


lm_RC <- lm(min_age - slope*somites ~ 1, 
            data = subset(data, species=="Rhamphochromis sp. 'chilingali'"))

# lm_AC <- lm(min_age ~ somites * species, data = data_sub)

new <- data.frame(somites = c(1, 33, 1, 39), 
                  species = c('Astatotilapia calliptera', 'Astatotilapia calliptera', 
                              "Rhamphochromis sp. 'chilingali'", 
                              "Rhamphochromis sp. 'chilingali'"))

# manually computed slopes and intercepts for each species
coef_table <- data.frame(
  species = c('Astatotilapia calliptera', "Rhamphochromis sp. 'chilingali'"),
  slope   = c( slope_AC , slope_RC ),
  intercept = c( coef(lm_AC)[1] , coef(lm_RC)[1] )
)


predict_df <- data.frame(
  somites = seq(5, 40, by = 5)
)


predict_df$AC <- coef_table$intercept[coef_table$species == "Astatotilapia calliptera"] +
  coef_table$slope[coef_table$species == "Astatotilapia calliptera"] *
  predict_df$somites

predict_df$RC <- coef_table$intercept[coef_table$species == "Rhamphochromis sp. 'chilingali'"] +
  coef_table$slope[coef_table$species == "Rhamphochromis sp. 'chilingali'"] *
  predict_df$somites


# merge coefficients into your data
new2 <- merge(new, coef_table, by = "species")

# compute min_age = slope * somites + intercept
new2$age <- new2$slope * new2$somites + new2$intercept

df_range <- new2 %>%
  group_by(species) %>%
  summarise(
    somite_min = min(somites),
    somite_max = max(somites),
    age_min = age[which.min(somites)],
    age_max = age[which.max(somites)]
  )

df_range %>%
  mutate(predicted_slope = (age_max - age_min) / (somite_max - somite_min))


plt <- ggplot(data, 
              aes(x = min_age, y = somites, 
                  shape = species, 
                  color = species)) + 
  geom_hline(yintercept=39, 
             color = cichlid_colours[2], 
             linetype = 'longdash', 
             linewidth = 1, show.legend = FALSE,
             alpha = 0.5) + 
  geom_hline(yintercept=33, 
             color = cichlid_colours[1], 
             linetype = 'longdash', show.legend = FALSE,
             linewidth = 1, 
             alpha = 0.5) + 
  geom_point(size = 2) + 
  geom_segment(
    data = df_range,
    aes(x = age_min, xend = age_max,
        y = somite_min, yend = somite_max, color = species),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = cichlid_colours) + 
  # stat_poly_eq(data = data_sub, 
  #              use_label("eq", "R2"), 
  #              label.x = "right", label.y = "bottom", 
  #              method = 'lm') + 
  theme_bw() + 
  coord_cartesian(clip = "off") +  # allow drawing outside plot area
  ylab('Somite number') + 
  xlab('Time (hpf)') + 
  ggtitle('B. Increased duration of segmentation') + 
  theme(legend.position = 'bottom', 
        legend.margin=margin(-10, 0, 0, 0), 
        plot.margin = margin(20, 20, 20, 20)
  ) + geom_segment(data = df_range, aes(x = age_min, xend = age_max,
                                        y = somite_max + 2, 
                                        color = species), 
                   linewidth = 4, lineend = "square", 
                   alpha = 0.5) + 
  geom_label(data = df_range, 
             aes(
               x = age_min, y = somite_max + 4,
               label = paste0(round(age_min, 0), 'hpf'), 
               color = species), 
             show.legend = FALSE, 
             size = 3
  ) + 
  geom_label(data = df_range, 
             aes(
               x = age_max, y = somite_max + 4,
               label = paste0(round(age_max, 0), 'hpf'), 
               color = species), 
             show.legend = FALSE, 
             size = 3
  )



#####
# SAVE 



design <- '
ABBB
CCDD
'

rate + plt +som_len + som_vol + 
  plot_layout(design = design, 
              # guides = 'collect'
              )  




ggsave('../graphs/Fig2.png',
       width = 10, height = 6, dpi =600)











data_sub %>% 
  ggplot(
    aes(x = min_age, y = somites,
        color = date, shape = date)) +
  geom_point() +
  scale_shape_manual(values=1:21, guide = 'none') +
  stat_poly_line(# data = data_sub, 
    se = TRUE, linewidth = 0.1) +
  # scale_colour_manual(values = cichlid_colours) +
  # stat_poly_eq(data = data_sub, use_label("eq", "R2")) +
  facet_grid(cols = vars(species)) + 
  theme_bw() + 
  labs(color = '') + 
  ylab('Somite number') + 
  xlab('Time (hpf)') + 
  ggtitle('Regressions on individual clutches') + 
  theme(legend.position = 'none') 

ggsave('../graphs/FigS2.png', width = 10, height = 5, dpi = 600)


  