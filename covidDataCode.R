library(tidyverse)
library(lubridate)
library(splines)
installed.packages("zoo")
library(zoo)
library(forecast) 
library(ggplot2)
library(ggpubr)
install.package("caTools")
library(caTools)

# initializations

cases_url <- "https://bioportal.salud.pr.gov/api/administration/reports/orders/basic"
# age_levels <-  paste(seq(0, 125, 5), "to", seq(4, 129, 5))
# the_years <- seq(2020, year(today()))
# imputation_delay  <- 2
# first_day <- make_date(2020, 3, 12)

# moving average ----------------------------------------------------------

ma7 <- function(d, y, k = 7) 
  tibble(date = d, moving_avg = as.numeric(stats::filter(y, rep(1/k, k), side = 1)))

sum7 <- function(d, y, k = 7) 
  tibble(date = d, moving_sum = as.numeric(stats::filter(y, rep(1, k), side = 1)))

# -- Fixed values
pr_pop <- 3285874 ## population of puerto rico

icu_beds <- 229 #if available beds is missing change to this

first_day <- make_date(2020, 3, 12)

last_complete_day <- today() - 1

the_years <- seq(2020, year(today()))

age_levels <-  c("0 to 9", "10 to 19", "20 to 29", "30 to 39", "40 to 49", "50 to 59", "60 to 69", 
                 "70 to 79", "80 to 89", "90 to 99", "100 to 109", "110 to 119", "120 to 129")

age_levels <-  paste(seq(0, 125, 5), "-", seq(4, 129, 5))

imputation_delay  <- 2

alpha <- 0.05

get_bioportal <- function(url){
  jsonlite::fromJSON(
    rawToChar(
      httr::GET(url, httr::content_type('application/json'),
                httr::add_headers('Accept-Enconding'="br"))$content)
  )
}

cases_url_molecular <-  paste0(cases_url,"?testType=Molecular&createdAtStartDate=2021-12-01T04:00:00Z&createdAtEndDate=2022-06-27T04:00:00Z")

caseras_url_antigens <-  paste0(cases_url,"?testType=AntigensSelfTest")

test_types <- c("Molecular", "AntigensSelfTest")


tab1 <- get_bioportal(cases_url_molecular)
tab2 <- get_bioportal(caseras_url_antigens)

tab <- rbind(tab1, tab2)

tab <- tab %>%  
  as_tibble() %>%
  mutate(collectedDate = ymd_hms(collectedDate, tz = "America/Puerto_Rico"),
         reportedDate = ymd_hms(reportedDate, tz = "America/Puerto_Rico"),
         orderCreatedAt = ymd_hms(orderCreatedAt, tz = "America/Puerto_Rico"),
         resultCreatedAt = ymd_hms(resultCreatedAt, tz = "America/Puerto_Rico"),
         ageRange       = na_if(ageRange, "N/A"),
         ageRange       = factor(ageRange, levels = age_levels),
         region = na_if(region, "N/A"),
         region = ifelse(region == "Bayamon", "Bayamón", region),
         region = ifelse(region == "Mayaguez", "Mayagüez", region),
         region = replace_na(region, "No reportada"),
         region = factor(region),
         result = tolower(result),
         result = case_when( 
           (grepl("covid", result) | grepl("sars-cov-2", result)) &
             grepl("positive", result) ~ "positive",
           grepl("influenza", result) ~ "other",
           grepl("positive", result) ~ "positive",
           grepl("negative", result) ~ "negative",
           result == "not detected" ~ "negative",
           TRUE ~ "other")) %>%
  arrange(reportedDate, collectedDate, patientId) 

## fixing bad dates
tab <- tab %>% 
  mutate(date = if_else(collectedDate > reportedDate, reportedDate, collectedDate)) %>% ## if collectedDate is in the future make it reportedDate
  mutate(date = if_else(is.na(collectedDate), reportedDate - days(imputation_delay),  collectedDate)) %>%
  mutate(date = if_else(!year(date) %in% the_years, reportedDate - days(imputation_delay),  date)) %>%  
  mutate(date = as_date(date)) %>%
  filter(year(date) %in% the_years & date <= today()) %>%
  arrange(date, reportedDate)

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

httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
url <- "https://covid19datos.salud.gov.pr/estadisticas_v2/download/data/sistemas_salud/completo"
hosp <- read.csv(text = rawToChar(httr::content(httr::GET(url)))) %>% 
  mutate(date = as_date(FE_REPORTE)) %>%
  filter(date >= first_day) %>%
  arrange(date) 

age_starts <- c(0, 10, 15, 20, 30, 40, 65, 75)
age_ends <- c(9, 14, 19, 29, 39, 64, 74, Inf)


url <- "https://covid19datos.salud.gov.pr/estadisticas_v2/download/data/defunciones/completo"
deaths <-  read.csv(text = rawToChar(httr::content(httr::GET(url)))) %>% 
  rename(ageRange = TX_GRUPO_EDAD, date = FE_MUERTE) %>%
  mutate(date = as_date(ymd_hms(date, tz = "America/Puerto_Rico"))) %>%
  mutate(age_start = as.numeric(str_extract(ageRange, "^\\d+")), 
         age_end = as.numeric(str_extract(ageRange, "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange, levels = age_levels)) 

# Rollmean graph
tab %>% group_by(testType) %>%
  group_by(testType, date) %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(pos = sum(result == 'positive'), n = n()) %>%
  mutate(seven_day_avg = rollmean(pos/n, 7, align = 'left', fill = 0)) %>%
  ggplot(aes(date, y = seven_day_avg)) + 
  geom_line(color = "blue", size = .7) +
  geom_point() +
  facet_wrap(~ testType)

tab3 <- tab %>% group_by(date, testType) %>% filter(testType == "Molecular") %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(pos = sum(result == "positive"), n =n())

tab4 <- tab %>% group_by(date, testType) %>% filter(testType != "Molecular") %>% 
  summarize(pos = sum(result == "positive"), n =n())


tab5 = tab3 %>% group_by(date) %>%
  summarize(avgM = pos / n)

tab6 = tab4 %>% group_by(date) %>%
  summarize(avgA = pos / n)

rm(tab3, tab4)

tabAvg = full_join(tab5, tab6)

rm(tab5, tab6)

tabAvg[is.na(tabAvg)] <- 0 

# Correlation between positivity rates
tabAvg %>% filter(date>=make_date(2021,12,1)) %>% 
  ggplot(aes(x = avgM, y = avgA)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.75)



## Binomial
# pr_pop <- 3285874
# 
# y1 = tab %>% filter(testType=="Molecular")
# 
# binaryData <- as.integer(as.data.frame(y1)$result == "positive")
# 
# positive_cases <- nrow(y1[y1$result=='positive',])
# n <- nrow(y1)
# 
# probability <- positive_cases / n
# 
# y <- rbinom(positive_cases, n, probability) / n
# hist(y)
# 
# ggplot(data=y) +
#   geom_histogram(aes(x=y))

## Binomial molecular

posM = tab %>% filter(testType=="Molecular" & date>=make_date(2021,12,15)) 
binaryData1 = as.integer(as.data.frame(posM)$result == "positive")

positive_cases1 <- nrow(posM[posM$result=='positive',])
n1 <- nrow(posM)

probability1 <- positive_cases1 / n

y <- rbinom(positive_cases1, n1, probability1) / n1
p1 = hist(y)
plot(p1)



#Binomial hometest

posA = tab %>% filter(testType!="Molecular")
binaryData2 = as.integer(as.data.frame(posA)$result == "positive")


positive_cases2 <- nrow(posA[posA$result=='positive',])
n2 <- nrow(posA)

probability2 <- positive_cases2 / n2

y <- rbinom(positive_cases2, n2, probability2) / n2

p2 = hist(y)

plot(p2)


# plot( p1, col=rgb(0,0,1,1/4))  # first histogram
# plot( p2, col=rgb(1,0,0,1/4), add=T) 


tab %>% group_by(date, testType) %>%  filter(date>=make_date(2021,12,1)) %>%
  summarize(pos=sum(result=="positive"), n=n()) %>%
  ggplot(mapping = aes(x = testType, y = pos/n)) +
  geom_violin() +
  geom_boxplot(width = 0.1) 
 


# n1 <- tab %>% group_by(date, testType) %>% filter(date>=make_date(2021,12,1)) %>%
#   summarize(pos=sum(result=="positive"), n=n()) 
# 
# v = n1$pos[n1$testType=="Molecular"]/n1$n[n1$testType=="Molecular"]
# v1 = n1$pos[n1$testType!="Molecular"]/n1$n[n1$testType!="Molecular"]

# Standard deviation
sd(binaryData1)
sd(binaryData2)


# Variance 
var(binaryData1)
var(binaryData2)

# Function standard error
se <- function(x) sqrt(var(x) / length(x))

se(binaryData1)
se(binaryData2)

# tt = tab %>% group_by(date) %>% filter(testType == "Molecular") %>%
#   filter(date>=make_date(2021,12,1)) %>%
#   summarize(n=n(), pos = sum(result=="positive"))
# 
# t1 = tab %>% group_by(date) %>% filter(testType != "Molecular") %>%
#   filter(date>=make_date(2021,12,1)) %>%
#   summarize(n=n(), pos = sum(result=="positive"))

# ?geom_ribbon
# 
# # Generate data
# huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
# h <- ggplot(t1, aes(date))
# 
# h + geom_ribbon(aes(ymin=0, ymax=pos/n))
# h + geom_area(aes(y = pos/n))
# 
# # Orientation cannot be deduced by mapping, so must be given explicitly for
# # flipped orientation
# h + geom_area(aes(x = pos/n, y = date), orientation = "y")
# 
# # Add aesthetic mappings
# h +
#   geom_ribbon(aes(ymin = pos/n - e2, ymax = pos/n + e2), fill = "grey70") +
#   geom_line(aes(y = pos/n))

# summary(tab5)

tests = tab %>% group_by(date, testType) %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(total=n(), pos = sum(result=="positive"), neg=sum(result=="negative"))

var(tests$pos[tests$testType=="Molecular"])
var(tests$pos[tests$testType!="Molecular"])

se(tests$pos[tests$testType=="Molecular"])
se(tests$pos[tests$testType!="Molecular"])


h = dnorm(v, mean(v), sd(v))

plot(v,h)

tab6 %>% group_by(date) %>%
  summarize(avgP = ifelse(avgA == 0, 0, avgA - e2)) %>%
  ggplot(aes(date, avgP)) +
  geom_point()

# Logistic regression model

tabM = tab %>% group_by(date) %>% filter(date>=make_date(2021,12,15) & date<make_date(2022,6,27)) %>% 
  filter(testType=="Molecular") %>% summarize(pos=sum(result=="positive"), n=n())
tabA = tab %>% group_by(date) %>% filter(date>=make_date(2021,12,1) & date<make_date(2022,6,27)) %>% 
  filter(testType!="Molecular") %>% summarize(pos=sum(result=="positive"), n=n())

# tabM = select(tabM,pos,n)
# tabA = select(tabA,pos,n)

tabA = tabA %>% group_by(date) %>%
  mutate(avgDaily = pos/n)

train_dataA = tabA %>% filter(date>=make_date(2021,12,15) & date<make_date(2022,3,22))
train_dataM = tabM %>% filter(date>=make_date(2021,12,15) & date<make_date(2022,3,22))

test_dataA = tabA %>% filter(date>=make_date(2022,3,22))
test_dataM = tabM %>% filter(date>=make_date(2022,3,22))


# split1<- sample(c(rep(0, 0.5 * nrow(tabA)), rep(1, 0.5 * nrow(tabA))))
# split1
# 
# train_dataA <- tabA[split1 == 0, ] 
# train_dataM <- tabM[split1 == 0, ] 
# 
# test_dataA <- tabA[split1 == 1, ] 
# test_dataM <- tabM[split1 == 1, ] 

# train_reg <- subset(tabA, split == "TRUE")
# test_reg <- subset(tabM, split == "FALSE")

model <- glm(cbind(pos, n-pos) ~ train_dataA$pos, data = train_dataM, 
family = binomial)

model

summary(model)

probabilities <- predict(model, full_join(test_dataA, test_dataM), type = "response")

tabM %>% filter(date>=make_date(2021,12,1)) %>%
  ggplot(aes(date,pos/n)) +
  geom_point()

probabilities

plot(probabilities) 
plot(tabM$date, tabM$pos/tabM$n)

# probabilities <- ifelse(probabilities >0.5, 1, 0)
# 
# table(test_dataA$avgDaily, probabilities)
# 
# fit <- model$fit
# fit2 <- model$family$linkinv(fit)


tab %>% group_by(region, testType) %>%
  summarize(pos = sum(result=="positive"), n = n()) %>%
  ggplot(aes(region,pos/n)) +
  geom_col() +
  facet_wrap(~ testType)

#GLM

model = glm(cbind(pos, n-pos) ~ tab4$pos, data = tab4, family=binomial)
probabilities <- predict(model, type = "response")

conf_int <- confint(model)

plot(probabilities)

mydata <- tab4 %>%
  dplyr::select_if(is.numeric)


predictors <- colnames(mydata)

mydata <- mydata %>%
  ungroup()%>%
  mutate(logit = log(probabilities / (1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

ggplot(mydata , aes(logit, predictor.value))+
  geom_point(size = 0.3, alpha = 0.5) +
  geom_smooth(method = "loess")
  

tab %>% group_by(date) %>% filter(result!="Molecular") %>%
  summarize(pos = count(result=="positive"))






