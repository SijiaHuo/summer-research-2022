library(lubridate)
library(tidyverse)

start_date <- make_date(2021, 12, 15)

pdat <- tab %>% #tab %>% 
  filter(result %in% c("positive", "negative")) %>%
  group_by(testType, date) %>%
  summarize(k = sum(result=="positive"),
            n = n(), patientId = patientId, result = result, .groups = "drop") %>%
  ungroup() %>%
  group_by(patientId) %>%
  filter(length(unique(testType))!=1) #gets rid of patients who have data on the same test type (no comparisons possible)

filtered_pdat <- pdat %>% 
  arrange(date) %>%
  group_by(patientId) %>%
  mutate(y=ifelse(testType=='Molecular', 1, 1.25)) 

filtered_pdat <- as.data.frame(filtered_pdat)
filtered_pdat$index <- seq.int(nrow(filtered_pdat))
filtered_pdat <- transform(filtered_pdat,                                 # Create ID by group
          patientIndex = as.numeric(factor(patientId))) # assign patient index by ID and x 5 to further separate each patient's data on plot

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
      
      if ((excludedRowTestType != rowTestType & abs(difftime(rowDate, excludedRowDate, units='days')) <= 7)) {
        merged_df <- as.data.frame(rbind(patientDataExcludingCurrent[numrow_excluded,], df[numrow,]))
        
        pairs <- rbind(pairs, merged_df)
      }
    }
  }
}

pairs <- distinct(pairs)

pairs <- transform(pairs, patientIndex = as.numeric(factor(patientId)))

pairs <- pairs %>%
  group_by(patientId) %>%
  mutate(y_val=y+patientIndex) %>%
  mutate(bottom=y_val-1, top=y_val+1)

sampleNumber <- 60
textSize <- 3

pairs %>%
  group_by(patientId) %>%
  filter(patientIndex<=sampleNumber) %>%
  filter(date >= start_date) %>%
  ggplot(aes(date, y_val, shape=testType, color=result)) +
  geom_hline(yintercept = seq.int(0.5, sampleNumber+1, 1), alpha=0.3, linetype='dotted') +
  geom_text(aes(label=substr(testType, 1, 1)), size=textSize, position = position_dodge(width=.2), fontface='bold') +
  theme_classic()
  
  
  
  
