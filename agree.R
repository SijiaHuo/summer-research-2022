library(lubridate)
library(tidyverse)

#runs the same function from the analysis.R file, but slightly modified 

pdat <- tab %>% #tab %>% 
  filter(date >= result %in% c("positive", "negative") & date>=first_day) %>%
  group_by(testType, date) %>%
  summarize(k = sum(result=="positive"),
            n = n(), patientId = patientId, .groups = "drop") %>% 
  select(patientId, testType, date, n, k) %>%
  ungroup() %>%
  group_by(patientId) %>%
  filter(as.integer(n_distinct(testType)) != 1, ) #gets rid of patients who have data on the same test type (no comparisons possible)

filtered_pdat <- pdat %>% 
  arrange(date) %>%
  mutate(lagDate = lag(date), days_elapsed = as.numeric(date - lagDate)) %>%
  filter(days_elapsed <= 4) %>% # excludes tests done that aren't within 4 days of lag date
  mutate(p = n/k)

filtered_pdat %>% 
  group_by(patientId) %>%
  ggplot(aes(date, p, shape=testType, color=testType)) +
  scale_shape_manual(values=c(3, 16, 17)) +
  geom_point() 
