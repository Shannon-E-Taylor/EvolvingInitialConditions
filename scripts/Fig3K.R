
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


#############
# MORPHOMETRICS # 
#################



psm_data = read.csv('../data/computed_psm_dimensions.csv')

psm_data$somite_count = as.integer(psm_data$somite_count)

psm_data = psm_data[
  (psm_data$species %in% c('A. calliptera', 'R. chillingali')) & (psm_data$somite_count > 0), 
]

drops = c('somite_count_y', 'species_y', 'somite_count_x', 'species_x', 
          # 'overview_id', 
          'X', 'xsize', 'zsize'#, #'psm_id'
)
psm_data = psm_data[, !names(psm_data) %in% drops]

psm_data = psm_data[!is.na(psm_data$species), ]



# method and params in config 
psm_data = assign_class(psm_data, 
                        AC_midseg_threshold, 
                        RC_midseg_threshold, 
                        AC_lateseg_threshold, 
                        RC_lateseg_threshold
)



cichlid_colours = c("A. calliptera"="#C3BE46FF", "R. chillingali"="#4D4D4DFF")



plot_absolute_stage <- function(data, y, species_diff){
  
  plt <- ggplot(data, 
                aes(x = somite_count, y = get(y), 
                    colour = species)
                ) +
    theme_bw() + 
    xlim(0, 40) + 
    xlab('Somite Stage') + 
    theme(legend.position = "bottom") + 
    scale_colour_manual(values = cichlid_colours) + 
    scale_y_continuous(expand = expansion(mult = 0.1)) 
  
  if (species_diff){
    plt <- plt + 
      stat_poly_line(se = FALSE) + 
      stat_poly_eq(use_label(c('eq', 'R2')), 
                   label.x = "right", 
                   label.y = c(0.95, 0.85))  + 
      geom_point()
  }
  
  else {
    plt <- plt + 
      stat_poly_line(se = FALSE, color = 'black', show.legend=FALSE) + 
      stat_poly_eq(color = 'black', 
                   use_label(c('eq', 'R2'))) + 
      geom_point(show.legend=FALSE)
  }
  
  
  return(plt)
}

psm <- plot_absolute_stage(psm_data, 'psm', species_diff = TRUE) + 
  ylab('PSM length (\u00B5m)')  

noto <- plot_absolute_stage(psm_data, 'noto', species_diff = TRUE) + 
  ylab('notocord length (\u00B5m)') 

mpz <- plot_absolute_stage(psm_data, 'mpz', species_diff = TRUE) + 
  ylab('MPZ length (\u00B5m)')

width <- plot_absolute_stage(psm_data, 'psm_width', species_diff = FALSE) + 
  ylab('PSM width (\u00B5m)') + theme(legend.position='none')

depth <- plot_absolute_stage(psm_data, 'psm_depth', species_diff = FALSE) + 
  ylab('PSM depth (\u00B5m)') + theme(legend.position='none')

som <- plot_absolute_stage(psm_data, 'som', species_diff = TRUE) + 
  ylab('Somite length (\u00B5m)') 

len.width <- plot_absolute_stage(psm_data, 'len_width', species_diff = FALSE) + 
  ylab('PSM length : width') 

len.depth <- plot_absolute_stage(psm_data, 'len_depth', species_diff = FALSE) + 
  ylab('PSM length : depth')

width.depth <- plot_absolute_stage(psm_data, 'width_depth', species_diff = FALSE) + 
  ylab('PSM width : depth') 

n <- plot_absolute_stage(psm_data, 'n', species_diff = FALSE) + 
  ylab('Notocord width (\u00B5m)')

som.psm <- plot_absolute_stage(psm_data, 'psm_som', TRUE) + 
  ylab('PSM : somite length') 

psm + noto + mpz + 
  width + depth + som +  
  len.width + len.depth + width.depth + 
  n + som.psm + 
  plot_layout(ncol = 3, guides = 'collect') & 
  theme(legend.position = 'bottom')


ggsave('../graphs/FigS4.png', width = 10, height = 10, dpi = 300)

#########
# PCA 
#########



psm_data_for_pca = psm_data
psm_data_for_pca = na.omit(psm_data_for_pca)

# psm_data_for_pca <- psm_data_for_pca[psm_data_for_pca$somite_count > 5, ]

psm_data_names =  names(psm_data_for_pca)
keeps = which(!psm_data_names %in% c('species', 'somite_count', 
                                     'class', 'propsom', 'overview_id', 
                                     'psm_id', 'max_stage'))



numerical_data = psm_data_for_pca[keeps]

res.pca <- prcomp(numerical_data, scale = TRUE)



# fviz_eig(res.pca, addlabels = TRUE)


pca_values = data.frame(predict(res.pca))

pca_values$species = psm_data_for_pca$species
pca_values$somite_count = psm_data_for_pca$somite_count

pca_values = assign_class(pca_values, 
                          AC_midseg_threshold, 
                          RC_midseg_threshold, 
                          AC_lateseg_threshold, 
                          RC_lateseg_threshold
)

pca_values$class <- as.character(pca_values$class)
pca_values$class[pca_values$somite_count <= 5] <- "onset"
pca_values$class <- factor(pca_values$class, levels = c("onset", "early", "mid", "late"))


# Calculate the hulls for each group
hull_cyl <- pca_values %>%
  group_by(species, class) %>%
  slice(chull(PC1, PC2))


# extract explained variance for each PC
PC1_perc = round(get_eigenvalue(res.pca)[1,2], digits=0)
PC2_perc = round(get_eigenvalue(res.pca)[2,2], digits=0)
PC3_perc = round(get_eigenvalue(res.pca)[3,2], digits=0)


pca <- ggplot() + 
  # Plot light grey points for all data
  geom_point(data = select(pca_values, c('PC1', 'PC2', 'species')), 
             aes(x = PC1, y = PC2, shape = species), 
             color = 'darkgrey', 
             alpha = 1, size = 1) + 
  # Draw polygons around groups (hull) with desired aesthetics
  geom_polygon(data = hull_cyl, 
               aes(x = PC1, y = PC2, linetype = species), 
               color = 'darkgrey', alpha = 0) + 
  # Add the colored points on top, with somite_count as color
  geom_point(data = pca_values, 
             aes(x = PC1, y = PC2, 
                 color = species, 
                 shape = species)) + 
  # Color scale for somite_count
  # scale_color_viridis_c() + 
  scale_colour_manual(values = cichlid_colours) + 
  theme_bw() + 
  facet_grid( ~ class) + 
  xlab(sprintf("PC1: %.f%% explained variance", PC1_perc)) + 
  ylab(sprintf("PC2: %.f%% explained variance", PC2_perc)) + 
  theme(legend.position = 'bottom') + 
  ggtitle('K. PCA separates embryos by stage and species')



ggsave('../graphs/Fig3K.png', width = 7, height = 3, dpi = 600)


# used in response to reviewers

# pca_values %>% 
#   count(species)

# 
# 
# pc1 <- ggplot(
#   pca_values, aes(x = somite_count, y = PC1, shape = species, color = species)
# ) + 
#   geom_point()
# 
# pc2 <- ggplot(
#   pca_values, aes(x = somite_count, y = PC2, shape = species, color = species)
# ) + 
#   geom_point()
# 
# pc1 + pc2 +  
#   plot_layout(guides = 'collect') &
#   scale_colour_manual(values = cichlid_colours)&
#   theme_bw() & 
#   theme(legend.position = 'bottom') & 
#   xlab('Somite count')
# 
# 
# ggsave('../graphs/PCA_stage.png', width = 7, height = 3, dpi = 600)
# 
# 

