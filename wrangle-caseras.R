# -- Libraries   
library(tidyverse)
library(lubridate)
library(splines)

if(grepl("fermat|leo", Sys.info()["nodename"])){
  rda_path <- "/homes10/rafa/dashboard/covidpr/rdas"
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

caseras_url_molecular <-  paste0(caseras_url,"?testType=AntigensSelfTest")
caseras_url_antigens <- paste0(caseras_url,"?testType=MolecularSelfTest")

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


caseras %>% group_by(date) %>%
  summarize(n=n()) %>%
  ggplot(aes(date, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas por día") +
  xlab("Fecha")

caseras %>% arrange(date) %>%
  group_by(date) %>% 
  summarize(n=n()) %>%
  mutate(n=cumsum(n)) %>%
  ggplot(aes(date, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas acumuladas") +
  xlab("Fecha")


caseras %>% group_by(date) %>%
  summarize(n=n(), p = sum(result == "positive")) %>%
  mutate(rate = p/n) %>%
  ggplot(aes(date, rate)) +
  geom_point()  +
  geom_smooth(method="loess", span = 1/10) +
  ylab("Tasa de positividad diaria") +
  xlab("Fecha")

caseras %>% group_by(date) %>%
  summarize(p = sum(result == "positive")) %>%
  mutate(rate = p) %>%
  ggplot(aes(date, p)) +
  geom_col()  +
  ylab("Positivos por día") +
  xlab("Fecha")


caseras %>% group_by(week = round_date(date, unit="week")) %>%
  summarize(n=n()) %>%
  ggplot(aes(week, n)) +
  geom_col()  +
  ylab("Numero de pruebas caseras reportadas por semana") +
  xlab("Fecha")

