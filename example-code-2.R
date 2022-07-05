library(tidyverse)
library(lubridate)

# initializations

cases_url <- "https://bioportal.salud.pr.gov/api/administration/reports/orders/basic"

# moving average ----------------------------------------------------------

ma7 <- function(d, y, k = 7) 
  tibble(date = d, moving_avg = as.numeric(stats::filter(y, rep(1/k, k), side = 1)))

sum7 <- function(d, y, k = 7) 
  tibble(date = d, moving_sum = as.numeric(stats::filter(y, rep(1, k), side = 1)))

first_day <- make_date(2021,12,15)

last_complete_day <- today() - 1

the_years <- seq(2020, year(today()))

age_starts <- c(0, 10, 15, 20, 30, 40, 65, 75)
age_ends <- c(9, 14, 19, 29, 39, 64, 74, Inf)

age_levels <- paste(age_starts, age_ends, sep = " a ")
age_levels[length(age_levels)] <- paste0(age_starts[length(age_levels)],"+")

imputation_delay  <- 2


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


tab <- bind_rows(tab1, tab2) %>%  
  as_tibble() %>%
  mutate(collectedDate = ymd_hms(collectedDate, tz = "America/Puerto_Rico"),
         reportedDate = ymd_hms(reportedDate, tz = "America/Puerto_Rico"),
         orderCreatedAt = ymd_hms(orderCreatedAt, tz = "America/Puerto_Rico"),
         resultCreatedAt = ymd_hms(resultCreatedAt, tz = "America/Puerto_Rico"),
         ageRange       = na_if(ageRange, "N/A"),
         age_start = as.numeric(str_extract(as.character(ageRange), "^\\d+")),
         age_end = as.numeric(str_extract(as.character(ageRange), "\\d+$"))) %>%
  mutate(ageRange = age_levels[as.numeric(cut(age_start, c(age_starts, Inf), right = FALSE))]) %>%
  mutate(ageRange = factor(ageRange)) %>%
  mutate(ageRange = forcats::fct_explicit_na(ageRange, "No reportada"),
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


tab <- distinct(tab) ## distinct remove duplicat rows


# wrangling ends,  analysis starts ----------------------------------------

xx <- tab %>% group_by(testType, patientId) %>%
  slice(1)

## data for model
dat <- xx %>% #tab %>% 
  filter(date >= result %in% c("positive", "negative") & date>=first_day) %>%
  group_by(testType, date, ageRange) %>%
  summarize(k = sum(result=="positive"),
            n = n(), .groups = "drop") %>% 
  select(ageRange, testType, date, n, k) %>%
  pivot_wider(names_from=testType, values_from = c(n, k)) %>% 
  group_by(ageRange) %>%
  mutate(k7_AntigensSelfTest = zoo::rollsum(k_AntigensSelfTest, 7, align = "right", na.pad = TRUE),
         n7_AntigensSelfTest = zoo::rollsum(n_AntigensSelfTest, 7, align = "right", na.pad = TRUE)) %>%
  mutate(k7_Molecular = zoo::rollsum(k_Molecular, 7, align = "right", na.pad = TRUE),
         n7_Molecular = zoo::rollsum(n_Molecular, 7, align = "right", na.pad = TRUE)) %>%
  ungroup() %>%
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
  facet_wrap(~ageRange)
  

dat %>% 
  mutate(lower = qbinom(0.025, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest,
         upper = qbinom(0.975, n7_AntigensSelfTest, p7_AntigensSelfTest)/n7_AntigensSelfTest) %>%
  ggplot(aes(date, p7_Molecular)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5) +
  geom_point(aes(y=p7_AntigensSelfTest))  +
  facet_wrap(~ageRange)


## Also note that when lack of fit starts, jump in tests:

dat %>% filter(!ageRange %in% c("0 a 9", "10 a 14", "15 a 19", "No reportada")) %>%
  ggplot(aes(date, n_AntigensSelfTest)) + 
  geom_col() +
  geom_line(aes(y=n7_AntigensSelfTest/7))  +
  facet_wrap(~ageRange)


## excluding minors

dat %>% filter(!ageRange %in% c("0 a 9", "10 a 14", "15 a 19", "No reportada")) %>%
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
fit <- glm(cbind(k7_Molecular, n7_Molecular) ~ p7_AntigensSelfTest, family = "binomial", data = dat)

mutate(dat, p_hat = predict(fit, newdata = dat, type = "response")) %>%
  ggplot(aes(date, p_Molecular)) +
  geom_point() +
  geom_line(aes(y=p_hat))
  


