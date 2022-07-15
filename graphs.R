library(lubridate)

# Positivity rate per day for each type of test
colors = c("Sunday" = 1, "Monday" = 2, "Tuesday" = 3, "Wednesday" = 4,
           "Thursday" = 5, "Friday" = 6, "Saturday" = 7)

tab %>% group_by(testType, date) %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(n=n(), pos = sum(result=="positive")) %>%
  mutate(wday = lubridate::wday(date,label=TRUE,abbr=FALSE)) %>%
  ggplot(aes(date, pos/n, color = as.factor(wday))) +
  geom_point() +
  facet_wrap(~ testType, ncol = 2) +
  labs(title = "Positivity Rate by Test Type", x = "Date", y = "Positive rate",
       color = "Days") + theme_bw()



dates = as.Date("2022-07-08")
lubridate::wday(dates, label = TRUE, abbr = FALSE)


tab %>% group_by(testType, date) %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(n=n(), pos = sum(result=="positive")) %>%
  mutate(wday = wday(date)) %>%
  ggplot(aes(date, pos/n, color = testType)) +
  geom_point() +
  geom_smooth()


# Bar plot of home tests reported per day
tab %>% group_by(date) %>%
  filter(testType != "Molecular") %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(n=n()) %>%
  ggplot(aes(date, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas por día") +
  xlab("Fecha")

# Bar plot of number of accumulated home tests per day
tab %>% arrange(date) %>%
  filter(testType != "Molecular") %>%
  filter(date>=make_date(2021,12,1)) %>%
  group_by(date) %>% 
  summarize(n=n()) %>%
  mutate(n=cumsum(n)) %>%
  ggplot(aes(date, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas acumuladas") +
  xlab("Fecha")

# home test daily positive rate
tab %>% group_by(date) %>%
  filter(testType != "Molecular") %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(n=n(), p = sum(result == "positive")) %>%
  mutate(rate = p/n) %>%
  ggplot(aes(date, rate)) +
  geom_point()  +
  geom_smooth(method="loess", span = 1/10) +
  ylab("Tasa de positividad diaria") +
  xlab("Fecha")

# Bar plot positive cases per day of home tests
tab %>% group_by(date, testType) %>% filter(date>=make_date(2021,12,1)) %>%
  summarize(p = sum(result == "positive")) %>%
  mutate(rate = p) %>%
  ggplot(aes(date, log2(p))) +
  geom_col()  +
  ylab("Positivos por día") +
  xlab("Fecha") +
  facet_wrap(~ testType)

# bar plot number of reported home tests per week
tab %>% group_by(week = round_date(date, unit="week")) %>% 
  filter(testType != "Molecular") %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(n=n()) %>%
  ggplot(aes(week, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas por semana") +
  xlab("Fecha")


tab %>% group_by(testType, week = round_date(date, unit="week")) %>% 
  filter(date>=make_date(2021,12,15)) %>%
  summarize(n=n(), pos = sum(result=="positive")) %>%
  #mutate(wday = wday(date)) %>%
  ggplot(aes(week, pos/n)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ testType, nrow = 2) 


# Boxplot per test

tab %>% group_by(date, testType) %>%  filter(date>=make_date(2021,12,1)) %>%
  summarize(pos=sum(result=="positive"), n=n()) %>%
  ggplot(mapping = aes(x = testType, y = pos/n)) +
  geom_violin() +
  geom_boxplot(width = 0.1) 


# Add dates to table 
dat = as.data.table(dat) %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(weekday = weekdays(date)) %>% 
  dcast(date + n_AntigensSelfTest + n_Molecular + k_AntigensSelfTest + k_Molecular + 
          k7_AntigensSelfTest + n7_AntigensSelfTest + k7_Molecular + n7_Molecular +
          p_AntigensSelfTest + p7_AntigensSelfTest + p_Molecular + p7_Molecular 
        ~ weekday, fun.aggregate = length)


deaths %>% group_by(date) %>% filter(date>=make_date(2021,12,1)) %>%
  summarize(Female = sum(CO_SEXO=="F"), Man = sum(CO_SEXO=="M"), total = n()) %>%
  ggplot() +
  geom_point(aes(date, Female), color = "red") +
  geom_point(aes(date, Man), color = "blue")






