**Can Antigen Home Tests be Used for Population-Level COVID-19
Surveillance?**

**Members:** <br>
Michael H. Terrefortes Rosado, University of Puerto Rico ‘23 <br>
Javlon Nizomov, University of Florida ‘24
 
**Faculty Mentor:** Rafael Irizarry, Professor of Biostatistics, DFCI, HSPH <br>
**Graduate Student Mentor:** Sijia Huo, HSPH



### Abstract

COVID-19 is an infectious respiratory virus that has had a large impact on global health during the last 2 years. For this research project, we investigated if COVID-19 antigen home tests, which are much more economical than molecular tests, can be used for population-level COVID-19 surveillance. Molecular test results are used to inform the positivity rate in Puerto Rico. However, since December 2021, the Department of Health in Puerto Rico has allowed citizens to report their COVID-19 home test results. We studied both tests to find similarities in the daily and weekly trend of the positivity rate. In addition, we looked at age populations and certain time periods in order to find similar patterns between tests. Although not yet as reliable as molecular tests, we find that self-reported antigen home test results can potentially be used for COVID-19 surveillance.


### Introduction

Can home tests be used as surveillance for COVID-19? Home tests have been reported since December 2021 in Puerto Rico. Home test reporting is important because it provides a way for researchers to measure the positivity rate and provide information to health officials who can take measures. Also, it informs citizens about what measures to take and gives hospitals time to prepare before surges in positive cases. Moreover, it can allow people to find other reliable ways to measure COVID-19 surges that are cheaper and more accessible; for example, using antigen self-tests.

COVID-19 is easily transmissible via air particles and can cause long-term health problems in people who get the virus. Older people are at higher risk for the more severe symptoms of COVID-19, which can result in hospitalization and even death. Long COVID-19 is also becoming problematic, where people experience long-term effects of their infection. While there are different types of COVID-19 tests, for this project, we focused on molecular tests and antigen home tests. 

Molecular tests can be used to detect RNA or DNA. These types of tests are useful for detecting emerging diseases because they are quick to develop and use. Through this, they can detect the presence of COVID-19, even during the early stages of infection. Molecular tests are important because they detect the viral RNA that can accumulate to high levels within a cell before virus particles are created. Some of the limitations of these tests include being too expensive. Also, it can have a positive result long after you are infected because these molecular tests detect the DNA of the virus, which can stay in your body for weeks.

Antigen home tests can detect proteins unique to the COVID-19 virus. Some of the benefits of using these tests include being easy to perform, cheap, and delivering results in a few minutes. However, the limitations are that it can miss newly-infected patients that have virus replicating in their cells that have not yet created measurable amounts of virus particles. Also, it detects a protein that is present with COVID-19 infection. Once the COVID-19 virus is no longer active, this protein will no longer get detected.

### Methods

The methods used to perform the project are focused on the positivity rate of each test type. The positivity rate is the percentage of COVID-19 tests that are positive. This is not an estimate of the population that is infected. It is more likely for an infected person to get a test, so the population that gets tested is not representative. However, the positivity rate appears to follow the same trend as incidence which makes it a useful metric for surveillance. The positivity rate can be simply calculated by dividing the total number of positive cases by the total cases.

We used data from the BioPortal database and to perform data analysis, we stratified the data into separate periods and looked at how well the data correlated. In addition, we stratified the data by region and age group to see how the home tests behave compared to the molecular tests. Pearson correlation was used to measure the correlation between the two tests for different periods in time. 

For our last analysis, we decided to measure how well the tests agree for specific patients. By gathering all the tests done for each patient ID and creating pairs for tests done within a week, we were able to create a plot that allowed us to make conclusions about how well the antigen and molecular test results agreed with each other if conducted within a seven-day window.


### Results

After plotting the daily positivity rate it was found that the data of both tests follow a similar pattern. The home test data is more scattered in comparison to the molecular tests. Also, the antigen home test had some outliers because there is little testing being reported by citizens. For example, on a particular day, there were only two home tests performed and both were positive. Therefore, we decided to smooth the positivity rate by performing a rollmean average of seven days. Figure 1 shows the scatter plot of the daily positivity rate. 

To smoothen out the data, we calculated a rolling average using a seven-day window to better see the correlation between the two positivity rate trends. The data from December to February had a better correlation than any other time period shown. 

In order to investigate potential covariates, we divided the data by age groups to perform a similar analysis. It was found that the correlation between the age ranges 23 to 28 and 30 to 39 was strong from December to February but then weakened afterward. This is likely due to a decrease in the number of home test results reported. 

The correlation between the two tests increased for those 40 to 64 years old, but then weakened significantly for those 65 to 74 years old. This is because those in the older age group reported their test results less frequently.

Then, we stratified the data by region to investigate how different regions compare. There were no observable differences between the positivity rate trends of different regions in Puerto Rico. Figure 5 shows Metro and Ponce.

The same analysis was done for Caguas and Mayagüez in Figure 6 and there were no observable differences found.

The first correlation was done between the months of December and February. This correlation was very good for the surveillance of COVID-19. The correlation had an R of 0.96, this is very good because it shows us how both test positivity rate is very similar in this time period. This can be seen in Figure 7. However, after this time period, the other two correlations did not have a good relationship. In Figure 8, the months of March and April had a correlation of R = 0.78. This means that the home tests should not be used for surveillance. Finally in Figure 9, from the months of May and June, the correlation performed worse than in the previous months and the R-value was 0.62. During this time period, it was seen a big increase in reported cases by citizens. 

To see how well the two tests compared, we plotted paired test data for 60 patients, where each value on the y-axis shows a different patient’s testing history. The blue results indicate a positive result, while the red results indicate a negative result. Based on these results, we concluded that, in general, the two tests had similar results if they were conducted within a week of each other.

### Discussion

We found that for certain time periods, and ages both tests agreed on the positivity rate in Puerto Rico. 
One limitation of our study included being limited to data from the Department of Health in Puerto Rico. This means that our findings may not be applicable to other regions of the world. Another limitation is that there were few home test results reported and this may have affected the trends in the graphs. This also led us to omit molecular home test results since there was not enough data from them. Furthermore, the positivity rate is not entirely representative of those who are infected because it is more likely for an infected person to get a test.


### Conclusion

Based on the plots, we found that the positivity rate appears to follow the same trend as incidence which makes it a useful metric for surveillance. This shows that home tests do have the potential to be used for COVID-19 surveillance. We encourage citizens to report their home test results so that the database can expand and allow for a more in-depth analysis in the future. However, since there was a limitation with having few home test results reported, further research is warranted.


### References


El Nuevo Día. (2022, January 19). ¿Cómo registrar los resultados de las pruebas Caseras del Covid-19 en el Bioportal de Salud? El Nuevo Día. Retrieved from https://www.elnuevodia.com/noticias/locales/notas/como-registrar-los-resultados-de-las-pruebas-caseras-del-covid-19-en-el-bioportal-de-salud/ 

Ledur, J. (2020, September 22). Test positivity: So valuable, so easy to misinterpret. The COVID Tracking Project. Retrieved from https://covidtracking.com/analysis-updates/test-positivity

Molecular tests: Types of COVID-19 tests. COVID-19 Testing Toolkit. (2022, March 1). Retrieved from https://www.centerforhealthsecurity.org/covid-19TestingToolkit/testing-basics/types-of-COVID-19-tests/diagnostic-tests/molecular-tests.html

Winny, A. (2020, November 2). What are all the different kinds of COVID-19 tests? Johns Hopkins Bloomberg School of Public Health. Retrieved from https://publichealth.jhu.edu/2020/what-are-all-the-different-kinds-of-covid-19-tests 
