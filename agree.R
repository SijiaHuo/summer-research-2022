library(lubridate)
library(tidyverse)

#runs the same function from the analysis.R file, but slightly modified 

pdat <- tab %>% #tab %>% 
  filter(result %in% c("positive", "negative") & date>=first_day) %>%
  group_by(testType, date) %>%
  summarize(k = sum(result=="positive"),
            n = n(), patientId = patientId, result = result, .groups = "drop") %>%
  ungroup() %>%
  group_by(patientId) %>%
  filter(length(unique(testType))!=1)

  #filter(as.integer(n_distinct(testType)) != 1) #gets rid of patients who have data on the same test type (no comparisons possible)

filtered_pdat <- pdat %>% 
  arrange(date) %>%
  mutate(lagDate = lag(date), days_elapsed = as.numeric(date - lagDate)) %>%
  mutate(p = k/n)

filtered_pdat2 <- filtered_pdat %>%
  group_by(patientId) %>%
  filter(days_elapsed<=7) %>%
  filter(length(unique(testType))==2)

sample <- unique(filtered_pdat2$patientId)[50:100]

filtered_pdat2 %>%
  filter(patientId %in% sample) %>%
  ggplot(aes(date, p, shape=testType, color=result)) +
  geom_text(aes(label=substr(testType, 1, 1)), size=3.5) +
  facet_wrap(~patientId) +
  theme_bw()

