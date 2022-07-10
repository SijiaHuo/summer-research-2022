library(lubridate)
library(tidyverse)
  
#filter for patient_ids with at least 1 unique test type in its history
pg_tab <- tab %>% group_by(patientId) %>%
  filter(as.integer(n_distinct(testType)) != 1,)

#split df by patient_id and get first 100 df to sample
pg_tab <- split(pg_tab, pg_tab$patientId)[1:100]

for (df in pg_tab) {
  #get differences between dates
  df <- mutate(df, LagDate = lag(df$date), Diff = as.numeric(df$date - LagDate)) 
  
  #date difference greater than 4 days, then remove df
  if (df$Diff[2] > 4) {
    rm(df)
  }
}
