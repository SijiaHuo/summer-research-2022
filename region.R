library(tidyverse)
library(lubridate)
library(reshape2)
library(data.table)
library('ggplot2')
library('ggpubr')
library(npreg)
library(sjPlot)
# analysis starts ----------------------------------------


## data for model
dat <- tab %>% 
  filter(date >= result %in% c("positive", "negative") & date>=first_day) %>%
  group_by(testType, date, region) %>%
  summarize(k = sum(result=="positive"),
            n = n(), .groups = "drop") %>% 
  select(region, testType, date, n, k) %>%
  pivot_wider(names_from=testType, values_from = c(n, k)) %>% 
  mutate(k7_AntigensSelfTest = zoo::rollsum(k_AntigensSelfTest, 7, align = "right", na.pad = TRUE),
         n7_AntigensSelfTest = zoo::rollsum(n_AntigensSelfTest, 7, align = "right", na.pad = TRUE)) %>%
  mutate(k7_Molecular = zoo::rollsum(k_Molecular, 7, align = "right", na.pad = TRUE),
         n7_Molecular = zoo::rollsum(n_Molecular, 7, align = "right", na.pad = TRUE)) %>%
  mutate(p_AntigensSelfTest = k_AntigensSelfTest / n_AntigensSelfTest,
         p7_AntigensSelfTest = k7_AntigensSelfTest / n7_AntigensSelfTest,
         p_Molecular = k_Molecular / n_Molecular,
         p7_Molecular = k7_Molecular / n7_Molecular)  


#p7 is the 7 day positivity rate that we want to estimate

## notice correlation much stronger with weekly average:

dat %>% 
  mutate(lower = qbinom(0.025, n_AntigensSelfTest, p_AntigensSelfTest)/n_AntigensSelfTest,
         upper = qbinom(0.975, n_AntigensSelfTest, p_AntigensSelfTest)/n_AntigensSelfTest) %>%
  ggplot(aes(date, p7_Molecular)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5) +
  geom_point(aes(y=p_AntigensSelfTest)) +
  facet_wrap(~region)


dat %>% 
  mutate(lower = qbinom(0.025, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest,
         upper = qbinom(0.975, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest) %>%
  ggplot(aes(date, p7_Molecular)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5) +
  geom_point(aes(y=p7_AntigensSelfTest)) +
  facet_wrap(~region)


## Also note that when lack of fit starts, jump in tests:

dat %>% filter(!region %in% c("No reportada")) %>%
  ggplot(aes(date, n_AntigensSelfTest)) + 
  geom_col() +
  geom_line(aes(y=n7_AntigensSelfTest/7))  +
  facet_wrap(~region)



dat %>% filter(!region %in% c("No reportada")) %>%
  group_by(date) %>%
  summarize(across(starts_with(c("n","k")), sum, na.rm=TRUE)) %>%
  mutate(p_AntigensSelfTest = k_AntigensSelfTest / n_AntigensSelfTest,
         p7_AntigensSelfTest = k7_AntigensSelfTest / n7_AntigensSelfTest,
         p_Molecular = k_Molecular / n_Molecular,
         p7_Molecular = k7_Molecular / n7_Molecular)  %>%
  mutate(lower = qbinom(0.025, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest,
         upper = qbinom(0.975, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest) %>%
  ggplot(aes(date, p7_Molecular)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5) +
  geom_point(aes(y=p7_AntigensSelfTest)) 


## does a linear model help?


fit <- glm(cbind(k7_Molecular, n7_Molecular) ~ p7_AntigensSelfTest, 
           family = "binomial", data = dat)

mutate(dat, p_hat = predict(fit, newdata = dat, type = "response")) %>%
  ggplot(aes(date, p7_Molecular)) +
  geom_point() +
  geom_line(aes(y=p_hat))
