# Positivity rate per day for each type of test
tab %>% group_by(testType, date) %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(n=n(), pos = sum(result=="positive")) %>%
  mutate(wday = wday(date)) %>%
  ggplot(aes(date, pos/n, color = as.factor(wday))) +
  geom_point() +
  facet_wrap(~ testType, nrow = 2) 


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
