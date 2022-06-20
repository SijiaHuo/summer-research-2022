# -- Libraries   
library(tidyverse)
library(lubridate)
library(splines)

if(grepl("fermat|leo", Sys.info()["nodename"])){
  rda_path <- "/Users/michaelterrefortes/Documents/CovidProject"
} else{
  rda_path <- "rdas"
}

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

imputation_delay  <- 2

alpha <- 0.05

## filter by date example: ?createdAtStartDate=2021-09-09T04:00:00Z&createdAtEndDate=2021-09-10T04:00:00Z
caseras_url <- "https://bioportal.salud.pr.gov/api/administration/reports/orders/basic"

#caseras_url_antigens <-  paste0(caseras_url,"?testType=AntigensSelfTest&createdAtStartDate=2022-05-01T04:00:00Z&createdAtEndDate=2022-06-18T04:00:00Z")
#caseras_url_molecular <- paste0(caseras_url,"?testType=MolecularSelfTest&createdAtStartDate=2022-05-01T04:00:00Z&createdAtEndDate=2022-06-18T04:00:00Z")

caseras_url_antigens <-  paste0(caseras_url,"?testType=AntigensSelfTest")
caseras_url_molecular <- paste0(caseras_url,"?testType=MolecularSelfTest")


get_bioportal <- function(url){
  jsonlite::fromJSON(
    rawToChar(
      httr::GET(url, httr::content_type('application/json'),
                httr::add_headers('Accept-Enconding'="br"))$content)
  )
}

test_types <- c("MolecularSelfTest", "AntigensSelfTest")

# Reading and wrangling cases data from database ---------------------------
age_levels <-  paste(seq(0, 125, 5), "to", seq(4, 129, 5))

message("Reading case data.")

caseras_molecular <- get_bioportal(caseras_url_molecular)
caseras_antigens <- get_bioportal(caseras_url_antigens)
caseras <- rbind(caseras_molecular, caseras_antigens)
rm(caseras_molecular, caseras_antigens); gc(); gc()

message("Processing case data.")

caseras <- caseras %>%  
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
  arrange(reportedDate, collectedDate, patientId) %>%
  filter(testType %in% test_types)

## fixing bad dates: if you want to remove bad dates instead, change FALSE TO TRUE
if(FALSE){
  ## remove bad dates
  caseras <- caseras %>% 
    filter(!is.na(collectedDate) & year(collectedDate) %in% the_years & collectedDate <= today()) %>%
    mutate(date = as_date(collectedDate))
} else{
  ## Impute missing dates
 caseras <- caseras %>% 
    mutate(date = if_else(collectedDate > reportedDate, reportedDate, collectedDate)) %>% ## if collectedDate is in the future make it reportedDate
    mutate(date = if_else(is.na(collectedDate), reportedDate - days(imputation_delay),  collectedDate)) %>%
    mutate(date = if_else(!year(date) %in% the_years, reportedDate - days(imputation_delay),  date)) %>%  
    mutate(date = as_date(date)) %>%
    filter(year(date) %in% the_years & date <= today()) %>%
    arrange(date, reportedDate)
}

# Bar plot of home tests reported per day
caseras %>% group_by(date) %>%
  summarize(n=n()) %>%
  ggplot(aes(date, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas por día") +
  xlab("Fecha")

# Bar plot of number of accumulated home tests per day
caseras %>% arrange(date) %>%
  group_by(date) %>% 
  summarize(n=n()) %>%
  mutate(n=cumsum(n)) %>%
  ggplot(aes(date, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas acumuladas") +
  xlab("Fecha")

# home test daily positive rate
caseras %>% group_by(date) %>%
  summarize(n=n(), p = sum(result == "positive")) %>%
  mutate(rate = p/n) %>%
  ggplot(aes(date, rate)) +
  geom_point()  +
  geom_smooth(method="loess", span = 1/10) +
  ylab("Tasa de positividad diaria") +
  xlab("Fecha")

# Bar plot positive cases per day of home tests
caseras %>% group_by(date) %>%
  summarize(p = sum(result == "positive")) %>%
  mutate(rate = p) %>%
  ggplot(aes(date, p)) +
  geom_col()  +
  ylab("Positivos por día") +
  xlab("Fecha")

# bar plot number of reported home tests per week
caseras %>% group_by(week = round_date(date, unit="week")) %>%
  summarize(n=n()) %>%
  ggplot(aes(week, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas por semana") +
  xlab("Fecha")


write.csv(caseras,"/Users/michaelterrefortes/Documents/CovidProject/hometest.csv", row.names = FALSE)

x1 <- caseras %>% group_by(testType) %>%
  summarize(p = sum(result == "positive"))

# Positive results per home test
ggplot(data = x1) +
  geom_col(aes(testType, p))

x2 <- caseras %>% group_by(date) %>%
    summarize(AntigenPos = sum(testType == "AntigensSelfTest" & result == "positive"), 
              MolecularPos = sum(testType == "MolecularSelfTest" & result == "positive"))

# Daily positive cases from both tests
ggplot(data = x2) + 
  geom_line(aes(date, AntigenPos)) +
  geom_line(aes(date, MolecularPos))


#hospital <- hosp %>% group_by(dates) %>%
#  mutate(dates = as.Date(dates)) %>%
#  summarize(hosp = ifelse(is.na(hospitalizaciones), 0, hospitalizaciones)) 

#hospital2 <- filter(hospital, dates > ymd("2021-12-01"))


#ggplot(data = hospital2) + 
#  geom_point(aes(dates, hosp))

#hospital2 %>% group_by(week = round_date(dates, unit="week")) %>%
#  summarize(hosp = sum(hosp)) %>%
#  ggplot(aes(week, hosp)) +
#  geom_smooth()  +
#  ylab("Numero de hospitalizaciones reportadas por semana") +
#  xlab("Fecha")

#hospital2 <- hospital2 %>% 
#  group_by(week = round_date(dates, unit="week")) %>%
#  summarize(hosp = sum(hosp))


#positive <- caseras %>% 
#  group_by(week = round_date(date, unit="week")) %>%
#  summarize(n=n()) 

#cases_hosp <- full_join(hospital2, positive)


#ggplot(data = cases_hosp) + 
#  geom_point(aes(week, hosp), color = "blue") +
#  geom_point(aes(week, n), color = "red")

#reinfection <- data.frame(table(caseras$patientId))
#reinfection <- reinfection[reinfection$Freq > 1,]

# home tests positive cases per age
caseras %>% group_by(ageRange) %>%
  summarize(p = sum(result == "positive")) %>%
  ggplot(aes(ageRange, p)) +
  geom_col()

# Home tests positive per region
caseras %>% group_by(region) %>%
  summarize(p = sum(result == "positive")) %>%
  ggplot(aes(region, p)) +
  geom_col()

# Bar plot of daily home tests 
caseras %>% group_by(date) %>%
  filter (year(date)>=2022) %>%
  summarize(n=n()) %>%
  ggplot(aes(date,n)) +
  geom_col()

# Home tests positive rate per day
caseras %>% group_by(date) %>% filter(result != "other") %>%
  filter (date>make_date(2022,5,1)) %>%
  summarize(n=n(), pos = sum(result=="positive")) %>%
  mutate(wday = wday(date)) %>%
  ggplot(aes(date, pos/n, color = as.factor(wday))) +
  geom_point() 

allTests <- rbind(all_tests_with_id, caseras)

# Positive results per test type
allTests %>% group_by(testType) %>% filter(result != "other") %>%
  filter (date>make_date(2022,5,1)) %>%
  summarize(pos = sum(result=="positive")) %>%
  ggplot(aes(testType, log2(pos))) +
  geom_col()


data <- allTests %>% group_by(date) %>%
  summarize(n=n()) 

allTests %>% group_by(date) %>% filter(result != "other") %>%
  filter (date>make_date(2022,5,1)) %>% filter(testType == "Antigens") %>%
  summarize(pos = sum(result=="positive"))


newData <- allTests %>% group_by(testType, date) %>%
  summarize(pos=sum(result=="positive"), totalType=n()) 

newData <- newData %>% full_join(data, by="date")

caseras %>% group_by(date) %>% filter(testType=="MolecularSelfTest") %>%
  summarize(n=n())

# Positive rate of all test per day
newData %>% group_by(date) %>%
  filter (date>make_date(2022,5,1)) %>%
  mutate(wday = wday(date)) %>%
  ggplot(aes(date, pos/n, color = as.factor(wday), shape=testType)) +
  geom_point()


caseras %>% group_by(date) %>% filter(result != "other") %>%
  filter (date>make_date(2022,5,1)) %>%
  summarize(n=n(), pos = sum(result=="positive")) %>%
  mutate(wday = wday(date)) %>%
  ggplot(aes(date, pos/n, color = as.factor(wday))) +
  geom_point()


ggplot(data = newData) +
  geom_point(mapping = aes(x = date, y = pos/totalType)) +
  facet_wrap(~ testType, nrow = 2)

# Comparison graph of positive rates per type test
newData %>% group_by(date) %>%
  filter (date>make_date(2022,5,1)) %>%
  mutate(wday = wday(date)) %>%
  ggplot(mapping = aes(x = date, y = pos/totalType, color = as.factor(wday))) +
  geom_point() +
  facet_wrap(~ testType, nrow = 2) 

allTests %>% group_by(date) %>%
  filter (date>make_date(2022,5,1)) %>%
  summarize(pos=sum(result=="positive"), n = n()) %>%
  mutate(wday = wday(date)) %>%
  ggplot(aes(date, pos/n,  color = as.factor(wday))) +
  geom_point()
  

deaths %>% group_by(date) %>%
  summarize(n=n()) %>%
  ggplot(aes(date,n)) +
  geom_point()






