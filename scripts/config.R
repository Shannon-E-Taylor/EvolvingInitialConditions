library('ggpmisc')

# config.R
cichlid_colours = c('#C3BE46FF', '#4D4D4DFF',  '#00CCFF') # calliptera, Ramphos, zebras
whole_embryo_fname = '../data/rois_phalloidin.csv'
tbox_fname = '../data/tailbud_rois.csv'
metadata_fname = '../data/whole_embryo_scan_metadata.csv'



tag_levels = list(c("A'", "A''", "B'", "B''", "C'"," C''", "D'", "D''"))

absolute_stage_theme <- list(
  stat_poly_line(se = TRUE),  
  theme_bw(), 
  xlim(0, 40), 
  xlab('Somite Stage'), 
  theme(legend.position = "none"), 
  scale_colour_manual(values = cichlid_colours)
  
)

relative_stage_theme <- list(
  stat_poly_line(se = TRUE), 
  theme_bw(), 
  xlim(0, 100), 
  xlab('Percent somitogenesis'), 
  scale_colour_manual(values = cichlid_colours), 
  theme(
    axis.text.y=element_blank()), 
  ylab(NULL)
)


AC_midseg_threshold=15
RC_midseg_threshold=17

AC_lateseg_threshold=25
RC_lateseg_threshold=30

assign_class <- function(data, 
                         AC_midseg_threshold, 
                         RC_midseg_threshold, 
                         AC_lateseg_threshold, 
                         RC_lateseg_threshold
){
  
  
  # set Class (stage) by somite stage. 
  data$class = 'early'
  
  data[
    (data$species == 'A. calliptera') & 
      (data$somite_count > AC_midseg_threshold)
    , ]$class = 'mid'
  
  data[
    (data$species == 'R. chillingali') & 
      (data$somite_count > RC_midseg_threshold)
    , ]$class = 'mid'
  
  
  data[
    (data$species == 'A. calliptera') & 
      (data$somite_count > AC_lateseg_threshold)
    , ]$class = 'late'
  
  data[
    (data$species == 'R. chillingali') & 
      (data$somite_count > RC_lateseg_threshold)
    , ]$class = 'late'
  
  data$class <- factor(data$class, 
                       levels = c('early', 'mid', 'late'))
  
  return(data)
  
}

not_all_na <- function(x) any(!is.na(x))

