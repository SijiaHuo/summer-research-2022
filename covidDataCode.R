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
library(data.table) # for dcast()

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


tab %>% group_by(testType, week = round_date(date, unit="week")) %>% 
  filter(date>=make_date(2021,12,15)) %>%
  summarize(n=n(), pos = sum(result=="positive")) %>%
  #mutate(wday = wday(date)) %>%
  ggplot(aes(week, pos/n)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ testType, nrow = 2) 

weekData = tab %>% group_by(testType, week = round_date(date, unit="week")) %>% 
  filter(date>=make_date(2021,12,15)) %>%
  summarize(n=n(), pos = sum(result=="positive"))


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
  mutate(seven_day_avg = rollmean(pos/n, 7, align = 'right', fill = 0)) %>%
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
  ggplot(aes(x = date, y = avgA)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.75)


tab %>% group_by(date, testType) %>%  filter(date>=make_date(2021,12,1)) %>%
  summarize(pos=sum(result=="positive"), n=n()) %>%
  ggplot(mapping = aes(x = testType, y = pos/n)) +
  geom_violin() +
  geom_boxplot(width = 0.1) 
 

# Standard deviation
sd(binaryData1)
sd(binaryData2)


# Variance 
var(binaryData1)
var(binaryData2)

# Function standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}

se(binaryData1)
se(binaryData2)


tests = tab %>% group_by(date, testType) %>%
  filter(date>=make_date(2021,12,1)) %>%
  summarize(total=n(), pos = sum(result=="positive"), neg=sum(result=="negative"))

var(tests$pos[tests$testType=="Molecular"])
var(tests$pos[tests$testType!="Molecular"])

se(tests$pos[tests$testType=="Molecular"])
se(tests$pos[tests$testType!="Molecular"])


# Logistic regression model daily

tabM = tab %>% group_by(date) %>% filter(date>=make_date(2021,12,15) & date<make_date(2022,6,23)) %>%
  filter(testType=="Molecular") %>% summarize(MolecularPos=sum(result=="positive"), MolecularTotal=n())
tabA = tab %>% group_by(date) %>% filter(date>=make_date(2021,12,1) & date<make_date(2022,6,23)) %>%
  filter(testType!="Molecular") %>% summarize(HTPos=sum(result=="positive"), HTTotal=n())


          
          # tabM = select(tabM,date,pos,n,avg)
          # tabA = select(tabA,pos,n)

tabA = tabA %>% group_by(date) %>%
  mutate(HTAvgDaily = HTPos/HTTotal)

tabM = tabM %>% group_by(date) %>%
  mutate(MolecularAvgDaily = MolecularPos/MolecularTotal)

          # train_dataA = tabA %>% filter(date>=make_date(2021,12,15) & date<make_date(2022,3,22))
          # train_dataM = tabM %>% filter(date>=make_date(2021,12,15) & date<make_date(2022,3,22))
          # 
          # test_dataA = tabA %>% filter(date>=make_date(2022,3,22))
          # test_dataM = tabM %>% filter(date>=make_date(2022,3,22))

MolecularAndHT = full_join(tabA,tabM)

# Add dates to table 
MolecularAndHT = as.data.table(MolecularAndHT) %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(weekday = weekdays(date)) %>% 
  dcast(date + HTPos + HTTotal + HTAvgDaily + MolecularPos + MolecularTotal + MolecularAvgDaily
        ~ weekday, fun.aggregate = length)


# Removes outlier
MolecularAndHT = MolecularAndHT[!(MolecularAndHT$date=="2022-03-12")]
MolecularAndHT = MolecularAndHT[!(MolecularAndHT$date=="2022-06-22")]

MolecularAndHT = MolecularAndHT[!(MolecularAndHT$date=="2022-03-10")]
MolecularAndHT = MolecularAndHT[!(MolecularAndHT$date=="2022-03-09")]

# Smooth hometest positivity rate
MA = forecast::ma(MolecularAndHT$HTAvgDaily, order = 7)

MolecularAndHT$HTAvgDaily = MA

# Divide data into test and train
split1<- sample(c(rep(0, 0.5 * nrow(MolecularAndHT)), rep(1, 0.5 * nrow(MolecularAndHT))))
split1

          # train_dataA <- tabA[split1 == 0, ]
          # train_dataM <- tabM[split1 == 0, ]
          # 
          # test_dataA <- tabA[split1 == 1, ]
          # test_dataM <- tabM[split1 == 1, ]

train_data <- MolecularAndHT[split1 == 0, ]
 
test_data <- MolecularAndHT[split1 == 1, ]

      # train_reg <- subset(tabA, split == "TRUE")
      # test_reg <- subset(tabM, split == "FALSE")

# Create model using (MPos, MNeg) ~ HTAvg

model <- glm(cbind(MolecularPos, MolecularTotal-MolecularPos) ~ HTAvgDaily 
             , data = train_data, family = binomial, weights = HTTotal)

        # model <- glm(MolecularAvgDaily ~ HTAvgDaily
        #              , data = MolecularAndHT, family = binomial, weights = HTTotal)
        
        # model <- glm(cbind(HTPos, HTTotal - HTPos) ~ MolecularAvgDaily
        #               , data = MolecularAndHT, family = binomial(link = "log"))
        #  

# Check model

model
conf_int <- confint(model)
summary(model)

# Prob with test data
probabilities <- predict(model, test_data, type = "response", se.fit = TRUE)


test_data$fit = probabilities$fit
test_data$se = probabilities$se.fit

probabilities
            # probabilities <- unname(probabilities)
            # probabilities
prob = probabilities$fit

mydata <- test_data %>%
  dplyr::select_if(is.numeric)

mydata <- mydata %>%
  ungroup() %>%
  mutate(logit = log(prob / (1 - prob))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

ggplot(mydata, aes(predictor.value, logit)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess")

            # library(boot)
            # inv.logit(-1.2)
            # 
            # fit = model$fit
            # fit2 = model$family$linkinv(fit)
            # 
            # fit3 = data.frame(fit = fit, date = MolecularAndHT$date)
            # 
            # fit4 = full_join(MolecularAndHT, fit3)
            # 
            # ggplot(fit4) +
            #   geom_point(aes(date, HTAvgDaily)) +
            #   geom_point(aes(date, MolecularAvgDaily), color = "red") +
            #   geom_point(aes(date, fit), color = "blue")

ggplot(MolecularAndHT, aes(x = date, y = MolecularAndHT$MolecularAvgDaily)) +
  geom_point() +
  geom_ribbon(data = test_data, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = test_data, aes(y = fit)) 

ggplot(test_data, aes(date, MolecularAvgDaily)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes( ymin = (fit - 1.96 * se), ymax = (fit + 1.96 * se)))

        # plot(test_data$date, prob, ylim = c(0,0.5)) 
        # lines(tabA$date, tabA$HTAvgDaily)
        # lines(tabM$date, tabM$MolecularAvgDaily, ylim = c(0,0.5))

ll.null <- model$null.deviance/-2 
ll.proposed = model$deviance/-2 

(ll.null - ll.proposed)/ ll.null
1-pchisq(model$null.deviance - model$deviance, df = length(model$coefficients)-1)

anova(model, model2, test = "chisq")

pres2 = residuals(model, type = "pearson")
plot(fitted(model)^(1/2), sqrt(pres2))

model$deviance / model$df.residual
sum(residuals(model, type='pearson')^2)/model$df.residual 




# Logistic regression model weekly

tabWeekM = weekData %>% group_by(week) %>%
  filter(testType=="Molecular") %>% summarize(MolecularPos=pos, MolecularTotal=n)
tabWeekA = weekData %>% group_by(week) %>%  
  filter(testType!="Molecular") %>% summarize(HTPos=pos, HTTotal=n)


tabWeekA = tabWeekA %>% group_by(week) %>% filter(week<make_date(2022,6,26)) %>%
  mutate(HTAvgDaily = HTPos/HTTotal)

tabWeekM = tabWeekM %>% group_by(week) %>% filter(week<make_date(2022,6,26)) %>%
  mutate(MolecularAvgDaily = MolecularPos/MolecularTotal)

# Join data
MolecularAndHT = full_join(tabWeekA,tabWeekM)

# Split into test and train
split1<- sample(c(rep(0, 0.5 * nrow(MolecularAndHT)), rep(1, 0.5 * nrow(MolecularAndHT))))
split1


train_data <- MolecularAndHT[split1 == 0, ]

test_data <- MolecularAndHT[split1 == 1, ]


# Model for weekly average
model2 <- glm(cbind(MolecularPos, MolecularTotal-MolecularPos) ~ HTAvgDaily
             , data = train_data, family = binomial(link = "log"))


model2
conf_int <- confint(model2)
summary(model2)

# Probabilities
probabilities <- predict(model2, test_data, type = "response", se.fit = TRUE)


test_data$fit = probabilities$fit
test_data$se = probabilities$se.fit

probabilities
        # probabilities <- unname(probabilities)
        # probabilities

ggplot(MolecularAndHT, aes(x = week, y = MolecularAndHT$MolecularAvgDaily)) +
  geom_point(position = position_jitter(width = 0.05, height = 0.05)) +
  geom_ribbon(data = test_data, aes(y = fit, ymin = fit - 1.96 * (se), ymax = fit + 1.96 * (se)),
              fill = "blue", alpha = 0.3) +
  geom_line(data = test_data, aes(y = fit)) 

ggplot(test_data, aes(week, MolecularAvgDaily)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes( ymin = fit - 1.96 * se, ymax = fit + 1.96 * se))


    # plot(test_data$date, probabilities, ylim = c(0,0.5)) 
    # plot(tabA$date, tabA$HTAvgDaily)
    # plot(tabM$date, tabM$MolecularAvgDaily, ylim = c(0,0.5))

ll.null <- model2$null.deviance/-2 
ll.proposed = model2$deviance/-2 

(ll.null - ll.proposed)/ ll.null
1-pchisq(model2$null.deviance - model2$deviance, df = length(model2$coefficients)-1)

anova(model, model2, test = "chisq")

pres2 = residuals(model2, type = "pearson")
plot(fitted(model2)^(1/2), sqrt(pres2))

model2$deviance / model2$df.residual
sum(residuals(model2, type='pearson')^2)/model2$df.residual 


# GLM but with rollmean of seven days

rollMeanAvg = tab %>% group_by(testType) %>%
  group_by(testType, date) %>%
  filter(date>=make_date(2021,12,15) & date<make_date(2022,6,23)) %>%
  summarize(pos = sum(result == 'positive'), n = n()) %>%
  mutate(seven_day_avg = rollmean(pos/n, 7, align = 'right', fill = 0))

rollMeanM = rollMeanAvg %>% filter(testType=="Molecular") %>% group_by(date) %>% 
  summarize(MPos = pos, MTotal = n, TRMean = seven_day_avg)
rollMeanA = rollMeanAvg %>% filter(testType!="Molecular") %>% group_by(date) %>% 
  summarize(HTPos = pos, HTTotal = n, HTRMean = seven_day_avg)


rollMean_HT_M = full_join(rollMeanM, rollMeanA)

ggplot(rollMean_HT_M, aes(date, HTRMean)) +
  geom_point()


MA = forecast::ma(rollMean_HT_M$HTRMean, order = 7)

rollMean_HT_M$HTRMean = MA

split1<- sample(c(rep(0, 0.5 * nrow(rollMean_HT_M)), rep(1, 0.5 * nrow(rollMean_HT_M))))
split1



train_data <- rollMean_HT_M[split1 == 0, ]

test_data <- rollMean_HT_M[split1 == 1, ]


model3 <- glm(cbind(MPos, MTotal-MPos) ~ HTRMean
             , data = train_data, family = binomial(link = "log"))

# model3 <- glm(cbind(HTPos, HTTotal - HTPos) ~ MolecularAvgDaily
#              , data = MolecularAndHT, family = binomial(link = "log"))


model3
conf_int <- confint(model3)
summary(model3)


probabilities <- predict(model3, test_data, type = "response", se.fit = TRUE)


test_data$fit = probabilities$fit
test_data$se = probabilities$se.fit

probabilities
          # probabilities <- unname(probabilities)
          # probabilities

ggplot(rollMean_HT_M, aes(x = date, y = rollMean_HT_M$TRMean)) +
  geom_point() +
  geom_ribbon(data = test_data, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = test_data, aes(y = fit))



ll.null <- model3$null.deviance/-2 
ll.proposed = model3$deviance/-2 

(ll.null - ll.proposed)/ ll.null
1-pchisq(model3$null.deviance - model3$deviance, df = length(model3$coefficients)-1)

anova(model3, model, test = "Chisq")

pres2 = residuals(model3, type = "pearson")
plot(fitted(model3)^(1/2), sqrt(pres2))

model3$deviance / model3$df.residual
sum(residuals(model3, type='pearson')^2)/model3$df.residual 



# probabilities <- ifelse(probabilities >0.5, 1, 0)
# 
# table(test_dataA$avgDaily, probabilities)
# 
# fit <- model$fit
# fit2 <- model$family$linkinv(fit)
# 
# critval <- 1.96 ## approx 95% CI
# uprHT = mean(MolecularAndHT$HTAvgDaily) + critval * (sd(MolecularAndHT$HTAvgDaily)/sqrt(MolecularAndHT$HTAvgDaily))
# uprM = mean(MolecularAndHT$MolecularAvgDaily) + critval * (sd(MolecularAndHT$MolecularAvgDaily)/sqrt(MolecularAndHT$MolecularAvgDaily))
# 
# lowHT = mean(MolecularAndHT$HTAvgDaily) - critval * (sd(MolecularAndHT$HTAvgDaily)/sqrt(MolecularAndHT$HTAvgDaily))
# lowM = mean(MolecularAndHT$MolecularAvgDaily) - critval * (sd(MolecularAndHT$MolecularAvgDaily)/sqrt(MolecularAndHT$MolecularAvgDaily))
#   
# 
# PCI <- function(x, z) 
#   mean(x) + z * (sd(x)/sqrt(length(x)))
# 
# NCI <- function(x, z)  mean(x) - z * (sd(x)/sqrt(length(x)))
# 
# 
# NCI(binaryData2, critval)
# 
# upr <- preds$fit + (critval * preds$se.fit)
# lwr <- preds$fit - (critval * preds$se.fit)
# fit <- preds$fit
# 
# 
# #GLM
# 
# model = glm(cbind(pos, n-pos) ~ tab4$pos, data = tab4, family=binomial)
# probabilities <- predict(model, type = "response")
# 
# conf_int <- confint(model)
# 
# plot(probabilities)
# 
# mydata <- tab4 %>%
#   dplyr::select_if(is.numeric)
# 
# 
# predictors <- colnames(mydata)
# 
# mydata <- mydata %>%
#   ungroup()%>%
#   mutate(logit = log(probabilities / (1-probabilities))) %>%
#   gather(key = "predictors", value = "predictor.value", -logit)
# 
# ggplot(mydata , aes(logit, predictor.value))+
#   geom_point(size = 0.3, alpha = 0.5) +
#   geom_smooth(method = "loess")
#   
# 
# tab %>% group_by(date) %>% filter(result!="Molecular") %>%
#   summarize(pos = count(result=="positive"))
# 
# 
# test = tab3
# 
# 
# # Other prediction
# 
# probabilities <- predict(model1, type = "response")
# plot(probabilities)
# 
# mydata <- test %>%
#   dplyr::select_if(is.numeric)
# 
# mydata <- mydata %>%
#   ungroup() %>%
#   mutate(logit = log(probabilities / (1 - probabilities))) %>%
#   gather(key = "predictors", value = "predictor.value", -logit)
# 
# ggplot(mydata, aes(logit,predictor.value)) +
#   geom_point(size = 0.5, alpha = 0.5) +
#   geom_smooth(method = "loess")
# 
# 
# split1<- sample(c(rep(0, 0.7 * nrow(tab6)), rep(1, 0.3 * nrow(tab6))))
# table(split1) 
# 
# train <- tab6[split1, ]
# test <- tab6[-split1, ]
# 
# #total molecular test - molecular positive
# 
# model1 <- glm(cbind(pos, n-pos) ~ pos, family = "binomial", data = test)
# 
# preddata <- with(tab6, data.frame(x = seq(min(pos), max(pos), length = 50)))
# preds <- predict(model1, type = "link", se.fit = TRUE)
# 
# critval <- 1.96 ## approx 95% CI
# upr <- preds$fit + (critval * preds$se.fit)
# lwr <- preds$fit - (critval * preds$se.fit)
# fit <- preds$fit
# 
# fit2 <- model1$family$linkinv(fit)
# upr2 <- model1$family$linkinv(upr)
# lwr2 <- model1$family$linkinv(lwr)
# 
# preddata$lwr <- lwr2 
# preddata$upr <- upr2 
# 
# tab6 %>% 
#   ggplot(aes(x=date,y=pos)) +
#   geom_ribbon(aes(ymin=pos-.2, ymax=pos+.2), size=.1, fill='gray') +
#   geom_point(size = 0.5, alpha = 0.5) +
#   geom_smooth(method = "loess")
# 




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

# Confidence interval
# 
# CI = function(x,z) {
#   mean(x) + z * (var(x)/sqrt(x))
# }
# 
# CIN = function(x,z) {
#   mean(x) - z * (var(x)/sqrt(x))
# }
# 
# var(tabA$HTAvgDaily)
# sd(tabA$HTAvgDaily)
# 
# pd = data.frame(date = tabA$date,
#   HTAvgDaily = tabA$HTAvgDaily, Pos = CI(tabA$HTAvgDaily, 1.96), Neg = CIN(tabA$HTAvgDaily, 1.96))
# 
# binaryData <- data.frame(binary = as.integer(as.data.frame(tab)$result == "positive"))
# 
# pd = data.frame(date = tab$date, binary = binaryData$binary, testType = tab$testType)
# 
# 
# 
# pd2 = pd %>% group_by(date, testType) %>% filter(date>=make_date(2021,12,15)) %>%
#   summarize(pos = sum(binary == 1), n=n(), se = se(binary)) %>%
#   mutate(avg = pos/n)
# 
# p1 = pd2 %>% filter(testType == "Molecular")
# p2 = pd2 %>% filter(testType != "Molecular")
# 
# 
# ggplot(p2, aes(x = date, y = avg)) +
#   #geom_point(position = position_jitter(width = 0.05, height = 0.05)) +
#   geom_ribbon(data = p2, aes(y = avg, ymin = avg - 1.96 * se, ymax = avg + 1.96 * se)
#           , fill = "red", alpha = 0.4) +
#   geom_line() +
#   geom_ribbon(data = p1, aes(y = p1$avg, ymin = p1$avg - 1.96 * p1$se, ymax = p1$avg + 1.96 * p1$se)
#               , fill = "blue", alpha = 0.4)
# 
# 
# ggplot(pd2, aes(x = date, y = avg)) +
#   #geom_point(position = position_jitter(width = 0.05, height = 0.05)) +
#   geom_ribbon(data = pd2, aes(y = avg, ymin = avg - 1.96 * se, ymax = avg + 1.96 * se)
#               , fill = "blue", alpha = 0.3) +
#   geom_line() +
#   facet_wrap(~ testType, nrow = 2)
# 
# ggplot() +
#   geom_point(data = pd2, aes(x = date, y = avg), position = position_jitter(width = 0.05, height = 0.05)) +
#   geom_ribbon(data = pd2, aes(y = avg, ymin = avg - 1.96 * se, ymax = avg + 1.96 * se), fill = "blue", alpha = 0.3) +
#   geom_line()
# 
# 
# C = tab %>% group_by(ageRange, testType) %>%
#   summarize(pos = sum(result=="positive"), n=n()) %>%
#   mutate(neg = n - pos)
# 
# len = length(test_data$HTAvgDaily)
# 
# sample(test_data$HTAvgDaily, len, replace = TRUE)
# 
# 
# unique(tab2$result)
# 
# 
# 
# S = duplicated(tab$patientId)
# 
# Ta = tab[S==TRUE, ]
# 
# Ta = Ta %>% filter(result=="positive")
# 

# de = deaths %>% group_by(date) %>% filter(date>=make_date(2021,12,1)) %>%
#   summarize(total = sum(CO_CLASIFICACION=="CONFIRMADO"))
# 
# ggplot(de) +
#   geom_point(aes(date, total))
# 
# ggplot(MolecularAndHT) +
#   geom_point(aes(date, MolecularAvgDaily))


# ## Binomial molecular
# 
# posM = tab %>% filter(testType=="Molecular" & date>=make_date(2021,12,15)) 
# binaryData1 = as.integer(as.data.frame(posM)$result == "positive")
# 
# positive_cases1 <- nrow(posM[posM$result=='positive',])
# n1 <- nrow(posM)
# 
# probability1 <- positive_cases1 / n
# 
# y <- rbinom(positive_cases1, n1, probability1) / n1
# p1 = hist(y)
# plot(p1)
# 
# 
# 
# #Binomial hometest
# 
# posA = tab %>% filter(testType!="Molecular")
# binaryData2 = as.integer(as.data.frame(posA)$result == "positive")
# 
# 
# positive_cases2 <- nrow(posA[posA$result=='positive',])
# n2 <- nrow(posA)
# 
# probability2 <- positive_cases2 / n2
# 
# y <- rbinom(positive_cases2, n2, probability2) / n2
# 
# p2 = hist(y)
# 
# plot(p2)
# 
# 
# # plot( p1, col=rgb(0,0,1,1/4))  # first histogram
# # plot( p2, col=rgb(1,0,0,1/4), add=T) 


# 
# MolecularAndHT %>% group_by(date) %>%
#   mutate(diff = HTAvgDaily - MolecularAvgDaily) %>%
#   ggplot(aes(date, MolecularAvgDaily)) +
#   geom_point() +
#   geom_errorbar(aes( ymin = ifelse(MolecularAvgDaily>HTAvgDaily,HTAvgDaily, MolecularAvgDaily), 
#                      ymax = ifelse(MolecularAvgDaily>HTAvgDaily, MolecularAvgDaily, HTAvgDaily))) +
#   geom_point(aes(date, HTAvgDaily), color = "blue")
# 
# 
# MA = forecast::ma(MolecularAndHT$HTAvgDaily, order = 7)
# 
# MolecularAndHT$HTAvgDaily = MA
# 
# ggplot() +
#   geom_line(aes(1:186,MolecularAndHT$HTAvgDaily)) +
#   geom_line(aes(1:186, MA), color = "blue") +
#   geom_line(aes(1:186, MolecularAndHT$MolecularAvgDaily), color="red")
# 
# lines(MolecularAndHT$HTAvgDaily)
# lines(MA, col="blue", lwd=3)
# 
# library(mgcv)
# 
# regM = gam(MolecularAvgDaily ~ HTAvgDaily, data = MolecularAndHT)
# summary(regM)
# plot(regM)

