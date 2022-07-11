library(lubridate)
library(tidyverse)

pdat <- tab %>% #tab %>% 
  filter(result %in% c("positive", "negative") & date>=first_day) %>%
  group_by(testType, date) %>%
  summarize(k = sum(result=="positive"),
            n = n(), patientId = patientId, result = result, .groups = "drop") %>%
  ungroup() %>%
  group_by(patientId) %>%
  filter(length(unique(testType))!=1) #gets rid of patients who have data on the same test type (no comparisons possible)

filtered_pdat <- pdat %>% 
  arrange(date) %>%
  mutate(lagDate = lag(date), days_elapsed = as.numeric(date - lagDate)) %>%
  group_by(patientId) %>%
  mutate(y=ifelse(testType=='Molecular', 1, 1.25))

filtered_pdat <- as.data.frame(filtered_pdat)
filtered_pdat$index <- seq.int(nrow(filtered_pdat))

split_data <- split(filtered_pdat, filtered_pdat$patientId)

pairs = data.frame()

for (df in split_data) {
  for (numrow in 1:nrow(df)) {
    rowTestType <- df[numrow,"testType"]
    rowResult  <- df[numrow, "result"]
    rowIndex  <- df[numrow, "index"]
    rowDate  <- df[numrow, "date"]
    
    patientDataExcludingCurrent <- filter(df, index!=rowIndex)
    
    for (numrow_excluded in 1:nrow(patientDataExcludingCurrent)) {
      excludedRowTestType <- patientDataExcludingCurrent[numrow_excluded,"testType"]
      excludedRowResult  <- patientDataExcludingCurrent[numrow_excluded, "result"]
      excludedRowIndex  <- patientDataExcludingCurrent[numrow_excluded, "index"]
      excludedRowDate <- patientDataExcludingCurrent[numrow_excluded, "date"]
      
      if ((excludedRowTestType != rowTestType & abs(difftime(rowDate, excludedRowDate, units='days')) <= 14)) {
        pairs <- rbind(pairs, as.data.frame(rbind(patientDataExcludingCurrent[numrow_excluded,], df[numrow,])))
      }
    }
  }
}

antigenY = 1
molecularY = 1.25

sampleSize = 50

filtered_pdat %>%
  filter(index <= sampleSize & date >= make_date(2021, 12, 22)) %>%
  ggplot(aes(date, y+index, shape=testType, color=result)) +
  geom_text(aes(label=substr(testType, 1, 1)), size=3.5) +
  scale_fill_gradient() +
  theme_bw()

