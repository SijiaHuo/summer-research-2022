#install.packages('zoo')

library(tidyverse)
library(lubridate)
library(splines)
library(zoo)

theme_set(theme_bw(16))
# initializations

cases_url <- "https://bioportal.salud.pr.gov/api/administration/reports/orders/basic"
age_levels <-  paste(seq(0, 125, 5), "to", seq(4, 129, 5))
the_years <- seq(2020, year(today()))
imputation_delay  <- 2
first_day <- make_date(2020, 3, 12)

get_bioportal <- function(url){
  jsonlite::fromJSON(
    rawToChar(
      httr::GET(url, httr::content_type('application/json'),
                httr::add_headers('Accept-Enconding'="br"))$content)
  )
}

cases_url_molecular <-  paste0(cases_url,"?testType=Molecular&createdAtStartDate=2022-05-01T04:00:00Z&createdAtEndDate=2022-06-17T04:00:00Z")
cases_url_antigens <- paste0(cases_url,"?testType=Antigens&createdAtStartDate=2022-05-01T04:00:00Z&createdAtEndDate=2022-06-17T04:00:00Z")

molecular_tab <- get_bioportal(cases_url_molecular)
antigen_tab <- get_bioportal(cases_url_antigens)
allTests <- rbind(molecular_tab, antigen_tab)

allTests <- allTests %>%  
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
allTests <- allTests %>% 
    mutate(date = if_else(collectedDate > reportedDate, reportedDate, collectedDate)) %>% ## if collectedDate is in the future make it reportedDate
    mutate(date = if_else(is.na(collectedDate), reportedDate - days(imputation_delay),  collectedDate)) %>%
    mutate(date = if_else(!year(date) %in% the_years, reportedDate - days(imputation_delay),  date)) %>%  
    mutate(date = as_date(date)) %>%
    filter(year(date) %in% the_years & date <= today()) %>%
    arrange(date, reportedDate)


allTests %>% group_by(date) %>% filter(result != "other") %>%
  filter(testType!="Molecular") %>%
  filter(date>=make_date(2022,5,1)) %>%
  summarize(n = n(), positive = sum(result=="positive")) %>%
  mutate(wday = wday(date)) %>%
  ggplot(aes(date, positive/n, color = as.factor(wday))) +
  geom_point() + xlab('Date') + ylab('Positivity Rate') + labs(colour="Weekday")

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

t1 <- allTests %>% group_by(date) %>%
  filter(date>=make_date(2022,5,1)) %>%
  summarize(pos = sum(result == 'positive', n = n()))

t1 %>%
  mutate(seven_day_avg = rollmean(t1$pos, 7, align = 'left', fill = 0)) %>%
  ggplot(aes(date, seven_day_avg)) + 
  geom_col(fill="blue") +
  geom_line(aes(y = seven_day_avg), color = "blue", size = .7)
