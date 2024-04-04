setwd("~/chinook_male_age/scripts")

library(tidyverse); library(ggplot2)

dat <- read.csv("../data/RADseq_sample_info.csv") %>% 
  select(c("fishID", "pop", "year", "age")) %>% 
  mutate(pop = as.factor(str_to_title(gsub("_", " ", pop)))) 

tab <- dat %>% 
  group_by(pop, age) %>% 
  tally()

write.csv(tab, "../data/sample_numbers.csv", row.names = FALSE)

ggplot(data = tab, aes(x = pop, y = age)) +
  geom_point(aes(size = n)) +
  scale_size(limits = c(0,100), 
             breaks = c(1, 25, 50, 75, 100)) +
  labs(x = NULL, y = "Age") + theme_bw() +
  theme(legend.title = element_blank())

