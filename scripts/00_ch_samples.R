setwd("~/chinook_male_age/scripts")

library(tidyverse); library(ggplot2); library(scatterpie); library(viridis)

# Read in data and tidy population labels.
dat <- read.csv("../data/chRADseq_samples.csv") %>% 
  select(c("fishID", "pop", "year", "age")) %>% 
  mutate(pop = as.factor(str_to_title(gsub("_", " ", pop)))) 

tab <- dat %>% 
  group_by(pop, age, year) %>% 
  tally()

# Summarise sample distributions and write.
write.csv(tab, "../data/sample_numbers.csv", row.names = FALSE)

# Proportions of samples by location and age by sampling year.
# Preparation for using geom_scatterpie.
props <- dat %>% 
  group_by(pop, age, year) %>% 
  tally() %>% 
  pivot_wider(names_from = year, values_from = n) %>%
  # Odd syntax here because column names are numbers/years.
  mutate(total = sum(`2016`, `2017`, `2018`, `2019`, `2020`, `2021`, na.rm = T),
         # For labeling categorical variables because geom_scatterpie requires numerical axes.
         popn  = as.numeric(as.factor(pop))) %>% 
  # Necessary to remove NAs beforehand.
  mutate_if(is.integer, ~replace(., is.na(.), 0)) %>% 
  # This makes the legend more sensible. 
  relocate(`2016`, .before = `2017`)

# Define scalar for modifying piechart sizes and legend.
scz <- 250
  
ggplot() + geom_scatterpie(data = props,
           # 300 is just a scalar for visualization.
           aes(x = popn, y = age, r = total/scz),
                  cols = colnames(props)[3:8]) +
  # Enables categorical labels on continuous axis.
  scale_x_continuous(breaks = c(1,2,3), 
                     labels = c("Chilliwack River Fall",
                                "Puntledge River Fall",
                                "Qualicum River")) +
  labs(x = NULL, y = "Age") + theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_fill_manual(values = alpha(c(viridis_pal(option = "H")(6)), 4/5)) +
  guides(fill = guide_legend(nrow = 1)) +
  # Add rectangle to clarify that legend refers to total sample size by age and location.
  geom_rect(aes(xmin = 6/11, xmax = 9/5, ymax = 5.6, ymin = 4.5),
            colour = "black", fill = "white") +
  # Add scatterpie radius legend inside rectangle defined above.
  # Seems difficult to get radius legend outside of plot area, but 
  # our data give us space for this anyway.
  geom_scatterpie_legend(radius = props$total/scz, x = 1, y = 5,
                         labeller = function(x) scz*x) +
  annotate("text", label = "Number of samples", x = 1.1, y = 5.5)

ggsave("../plots/sample_dist.tiff", dpi = 300, width = 7, height = 10)

