library('ggplot2')


update_geom_defaults("text", list(size = 3))


data <- read.csv('../data/pharyngula_stage_somite_counts.csv')
data <- data[data$somite_count > 0, ]


metadata <- read.csv('../data/pharyngula_stage_scans_metadata.csv')
metadata$species <- NULL

data <- merge(metadata, data, by.x = 'image_id', by.y = 'id')

df = data[data$species == 'R. chillingali', ]

mean_somite <- mean(df$somite_count)
sd_somite <- sd(df$somite_count)

# Create the histogram using ggplot2
plt <- ggplot(df, aes(x = somite_count)) +
  geom_histogram(binwidth = 1, fill = "#4D4D4DFF", color = "white") +
  labs(x = "Somite number", y = "Frequency") +
  annotate(
    "text",
    x = mean_somite, y = 2.5, vjust = -1,
    label = paste("Mean: ", round(mean_somite, 1), " ", 
                  "\u00B1", " ", round(sd_somite, 1)),
    color = "black"
  ) + 
  ggtitle("R. chillingali\n somite counts") +
  theme_bw()


# Create a data frame
df = data[data$species == 'A. calliptera', ]

mean_somite <- mean(df$somite_count)
sd_somite <- sd(df$somite_count)

my_y_title <- expression(paste(italic("A. calliptera"), " \n somite number"))

# Create the histogram using ggplot2
plt <- ggplot(df, aes(x = somite_count)) +
  geom_histogram(binwidth = 1, fill = "#C3BE46FF", color = "white") +
  labs(x = "Somite number", y = "Frequency") +
  annotate(
    "text",
    x = mean_somite, y = 2.5, vjust = -1,
    label = paste("Mean: ", round(mean_somite, 1), " ", 
                  "\u00B1", " ", round(sd_somite, 1)),
    color = "black"
  ) + 
  ggtitle("A. calliptera\n somite counts") + theme_bw()


source("config.R")
# 

data<- data[data$species %in% c('A. calliptera', 'R. chillingali'), ]


ac_mean_som <- mean(data[data$species == 'A. calliptera',]$somite_count)
ac_sd_som <- sd(data[data$species == 'A. calliptera',]$somite_count)

rc_mean_som <- mean(data[data$species == "R. chillingali",]$somite_count)
rc_sd_som <- sd(data[data$species == "R. chillingali",]$somite_count)

ac <- paste("A. calliptera: ", round(ac_mean_som, 1), "", 
            "\u00B1", " ", round(ac_sd_som, 1))

rc <- paste("R. chillingali: ", round(rc_mean_som, 1), "", 
            "\u00B1", " ", round(rc_sd_som, 1)) 

counts <- ggplot(data, aes(x = somite_count, fill = species)) + 
  geom_histogram(binwidth = 1, color = "white", position = "stack",) + 
  scale_fill_manual(values = cichlid_colours) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  annotate('text', x=-Inf, y=Inf, vjust=1.2, hjust=-.1, 
           label = paste(ac, rc, sep = '\n'), 
           # color = cichlid_colours[1]
  ) + 
  # annotate('text', x=Inf, y=Inf, vjust=2, hjust=1.1, 
  #          label = paste(rc), 
  #          # color = cichlid_colours[1]
  # ) + 
  ggtitle("Somite counts") + 
  labs(x = "Somite number", y = "Frequency") + 
  theme(legend.margin=margin(-10, 0, 0, 0)) 

ggsave('../graphs/FigS1.png', counts, width = 900, height = 900, units = 'px')

