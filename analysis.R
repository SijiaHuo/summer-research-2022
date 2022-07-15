library(tidyverse)
library(lubridate)
library(reshape2)
library(data.table)
library('ggplot2')
library('ggpubr')
library(npreg)

library(sjPlot)
library(forecast)
# analysis starts ----------------------------------------


## data for model
dat <- tab %>% 
  filter(date >= result %in% c("positive", "negative") & date>=first_day) %>%
  group_by(testType, date) %>%
  summarize(k = sum(result=="positive"),
            n = n(), .groups = "drop") %>% 
  select(testType, date, n, k) %>%
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
  geom_point(aes(y=p_AntigensSelfTest))


colors = c("Molecular" = "blue", "Home Test" = "red")

dat %>% 
  mutate(lower = qbinom(0.025, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest,
         upper = qbinom(0.975, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest) %>%
  ggplot(aes(date, p7_Molecular, color = "Molecular")) +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = "black"), width = 0.5) +
  geom_point(aes(y=p7_AntigensSelfTest, color = "Home Test")) + 
  ggtitle("Molecular and Home Tests Positivity Rate") +
  scale_color_manual(values = colors) + 
  labs(x = "Date", y = "7 day Molecular and Home Test Positivity Rate", color = "Test Type") +
  theme_bw()


## Also note that when lack of fit starts, jump in tests:

dat %>% ggplot(aes(date, n_AntigensSelfTest)) + geom_col() +geom_line(aes(y=n7_AntigensSelfTest/7))
## does a linear model help?


fit <- glm(cbind(k7_Molecular, n7_Molecular) ~ p7_AntigensSelfTest, 
           family = "binomial", data = dat)

mutate(dat, p_hat = predict(fit, newdata = dat, type = "response")) %>%
  ggplot(aes(date, p7_Molecular)) +
  geom_point() +
  geom_line(aes(y=p_hat), color = "red") +
  geom_line(aes(y=p7_AntigensSelfTest), color = "blue")


## Correlation
dat %>% filter(date>=make_date(2022,3,1) & date<make_date(2022,5,1)) %>%
  ggplot(aes(x = p7_Molecular, y = p7_AntigensSelfTest)) + 
  geom_point() +
  geom_abline() +
  geom_smooth(method='lm', formula= y~x) +
  # stat_cor(method = "pearson", label.x = 0.2, label.y = 0.1) +
  labs(title="Correlation Between Tests from March to April", 
  x = "Molecular 7 day avg", y = "Home test 7 day avg") + theme_bw() + xlim(0,0.4) + ylim(0,0.4)

dat %>% filter(date>=make_date(2021,12,1) & date<make_date(2022,3,1)) %>%
  ggplot(aes(x = p7_Molecular, y = p7_AntigensSelfTest)) + 
  geom_point() +
  geom_abline() +
  geom_smooth(method='lm', formula= y~x) +
  # stat_cor(method = "pearson", label.x = 0.05, label.y = 0.3) +
  labs(title="Correlation Between Tests from December to February", 
       x = "Molecular 7 day avg", y = "Home test 7 day avg") + theme_bw() + xlim(0,0.4) + ylim(0,0.4)

dat %>% filter(date>=make_date(2022,5,1) & date<make_date(2022,7,1)) %>%
  ggplot(aes(x = p7_Molecular, y = p7_AntigensSelfTest)) + 
  geom_point() +
  geom_abline() +
  geom_smooth(method='lm', formula= y~x) +
  # stat_cor(method = "pearson", label.x = 0.4, label.y = 0.3) +
  labs(title="Correlation Between Tests from May to June", 
       x = "Molecular 7 day avg", y = "Home test 7 day avg") + theme_bw() + xlim(0.2,0.6) + ylim(0.2,0.6)


dat[is.na(dat)] <- 0

deaths %>% group_by(date) %>% filter(date>=make_date(2021,12,15)) %>%
  summarize(total = n()) %>%
  ggplot(aes(date, total)) +
  geom_line()


#-------------------------------------------------



## Correlation
dat %>%  
  ggplot(aes(x = p7_Molecular, y = p7_AntigensSelfTest)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.75)


# Region 


#region = c("Metro", "Bayamón", "Arecibo", "Ponce", "Caguas", "Fajardo", "Mayagüez")
population = c(Metro = 682054,
               Bayamón = 537123,
               Arecibo = 394774,
               Ponce = 474603,
               Caguas = 529505,
               Fajardo = 116148,
               Mayagüez = 459487)

regionalPopulation = data.frame(region = names(population), population  = unname(population))

molecularRegion = tab %>% group_by(region) %>% filter(testType == "Molecular") %>%
  summarize(n = n())
  

HomeTestRegion = tab %>% group_by(region) %>% filter(testType != "Molecular") %>%
  summarize(n = n())

molecularRegion = full_join(molecularRegion, regionalPopulation)

HomeTestRegion = full_join(HomeTestRegion, regionalPopulation)

molecularRegion = molecularRegion %>% group_by(region) %>%
  mutate(rate = (n/population) * 1000)

HomeTestRegion = HomeTestRegion %>% group_by(region) %>%
  mutate(rate = (n/population) * 1000)


x1 = deaths %>% group_by(date) %>% filter(date>=make_date(2021,12,15)) %>%
  summarize(n = sum(CO_CLASIFICACION=="CONFIRMADO")) %>%
  mutate(wday = lubridate::wday(date,label=TRUE,abbr=FALSE)) %>%
  ggplot(aes(date, n, color = as.factor(wday))) +
  geom_point() + theme_bw() + labs(title = "Deaths by COVID-19 From December 2021 to July 2022"
                                   , x = "Date", y = "Number of deaths", color = "Days")

x2 = hosp %>% group_by(FE_REPORTE) %>%
  mutate(FE_REPORTE = as.Date(FE_REPORTE)) %>%
  mutate(wday = lubridate::wday(date,label=TRUE,abbr=FALSE)) %>%
  ggplot(aes(FE_REPORTE, CAMAS_ADULTOS_COVID, color = as.factor(wday))) +
  geom_point() + theme_bw() + labs(title = "Adult Hospitalization by COVID-19 from December 2021 to
                                   July 2022", x = "Date", y = "Hospitalizations", color = "Days")

x3 = dat %>% group_by(date) %>%
  mutate(wday = lubridate::wday(date,label=TRUE,abbr=FALSE)) %>%
  ggplot(aes(date, p_Molecular, color = as.factor(wday))) +
  geom_point() + theme_bw() + labs(title = "Positivity Rate Molecular Tests", x = "Date", y = "Positivity rate", color = "Days")

library(cowplot)

plot_grid(x1, x2, x3)


#write.csv(tab, "/Users/michaelterrefortes/Documents/CovidProject/tab.csv" )

