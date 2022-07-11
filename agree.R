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

filtered_pdat <- filtered_pdat %>% 
  filter(days_elapsed<=14) %>%
  group_by(patientId, testType) %>%
  mutate(agrees = ifelse(sum(testType == 'AntigensSelfTest') == sum(testType=='Molecular'), TRUE, FALSE))

sample <- unique(filtered_pdat$patientId)[1:25]

filtered_pdat %>%
  filter(patientId %in% sample) %>%
  ggplot(aes(date, p, shape=testType, color=result)) +
  geom_text(aes(label=substr(testType, 1, 1)), size=3) +
  facet_wrap(~patientId) +
  theme_bw()

