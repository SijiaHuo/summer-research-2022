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
  mutate(y=ifelse(testType=='Molecular', 1, 350)) # spacing between A and M for same patient based on test type

filtered_pdat <- as.data.frame(filtered_pdat)
filtered_pdat$index <- seq.int(nrow(filtered_pdat))
filtered_pdat <- transform(filtered_pdat,                                 # Create ID by group
          patientIndex = as.numeric(factor(patientId))) # assign patient index by ID and x 5 to further separate each patient's data on plot

split_data <- split(filtered_pdat, filtered_pdat$patientId)

pairs = data.frame()

p_index_counter = 1

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
        df$patientIndex <- p_index_counter
        patientDataExcludingCurrent$patientIndex <- p_index_counter
        
        merged_df <- as.data.frame(rbind(patientDataExcludingCurrent[numrow_excluded,], df[numrow,]))
        
        pairs <- rbind(pairs, merged_df)
      }
    }
  }
  p_index_counter = p_index_counter+150 # used to help space out each patient's data
}

pairs <- distinct(pairs)

pairs %>%
  group_by(patientId) %>%
  filter(patientIndex<=14000) %>%
  filter(date >= start_date) %>%
  ggplot(aes(date, y+(patientIndex), shape=testType, color=result, group=patientId)) +
  geom_text(aes(label=substr(testType, 1, 1)), size=5, position = position_dodge(width=.2), fontface='bold', family = "Arial", alpha=0.85) +
  theme_bw()

