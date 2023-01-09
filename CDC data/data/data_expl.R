library(tidyverse)
library(readxl)
library(readr)
library(lubridate)

setwd("/Users/madeleinekline/Dropbox (Harvard University)/G1/GradLab/multidemic/CDC data")
#will need to change on your device

flu <- read_csv("AgeViewByWeek.csv")
#flu surv net; > 70 counties and 13 states from EIP and  
#Influenza Hospitalization Surveillance Program
#has flu data by strain
unique(flu$`Age Group`) #only age breakdowns are 0-4, 5-24, 25-64, 65+
#2022 season complete (end of 2021, full 2022)
flu <- flu |> 
  mutate(all_strains = apply(flu[,4:12],1, sum, na.rm = TRUE))
#don't use this, not sure what exactly it is. Better to use EIP and IHSP

flu <- read_csv("FluSurveillance_Custom_Download_Data_header_removed.csv")
#weekly rate and cumulative rate of hospitalization per 100,000
#some race and sex data
flu <- flu[1:26634,]
unique(flu$CATCHMENT) #16 areas, some are EIP some IHSP, so maybe just stick to EIP?
unique(flu$`AGE CATEGORY`) #ages; overall, 0-4, 5-17, 18-49, 50-64, 65+, 65-74, 75-84,
#85+, 18-29, 30-39, 40-49, 5-11, 12-17, <18, >18 so some overlap here
flu_eip <- flu |> filter(NETWORK == "EIP")
unique(flu_eip$`AGE CATEGORY`)

#RSV weekly hospitalization rates
rsv <- read_csv("Weekly_Rates_of_Laboratory-Confirmed_RSV_Hospitalizations_from_the_RSV-NET_Surveillance_System (1).csv")
unique(rsv$State) #not all of them overlap, fewer in RSV than in flu


covid <- read_csv('United_States_COVID-19_Community_Levels_by_County.csv')
#all counties in US, has hospitalization rate need to figure out which counties map onto EIP and filter
#for now can filter to all counties in the states in EIP

cali_flu <- flu_eip |> 
  filter(CATCHMENT == "California") |>
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr"))

# cali_flu |> mutate(new_num = ifelse(`MMWR-WEEK` %in% 1:39, `MMWR-WEEK`,
#                                     ifelse(`MMWR-WEEK` == 52, -1,
#                                            ifelse(`MMWR-WEEK` == 51, -2,
#                                                   ifelse(`MMWR-WEEK` == 50, -3,
#                                                          ifelse(`MMWR-WEEK` == 49, -4,
#                                                                 ifelse(`MMWR-WEEK` == 48, -5,
#                                                                        ifelse(`MMWR-WEEK` == 47, -6,
#                                                                               ifelse(`MMWR-WEEK` == 46, -7,
#                                                                                      ifelse(`MMWR-WEEK` == 45, -8,
#                                                                                             ifelse(`MMWR-WEEK` == 44, -9,
#                                                                                                    ifelse(`MMWR-WEEK` == 43, -10,
#                                                                                                           ifelse(`MMWR-WEEK` == 42, -11,
#                                                                                                                  ifelse(`MMWR-WEEK` == 41, -12, 
#                                                                                                                         ifelse(`MMWR-WEEK` == 40, -13, NA))))))))))))))))

cali_flu |> 
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  mutate("Weekly_rate" = as.numeric(`WEEKLY RATE`)) |>
  ggplot(aes(`MMWR-WEEK`, Weekly_rate, group = `AGE CATEGORY`)) + 
  geom_line(aes(color = `AGE CATEGORY`)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  scale_y_continuous("Weekly Hospitalizations per 100,000 People") +
  theme_bw() + 
  labs(color = "Age Category") + 
  ggtitle("California EIP Site Flu Hospitalization Rates 2022 Season")
#this works for california, do it for all sites with facet

flu_by_age_eip <- flu_eip |> 
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr")) |>
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  mutate("Weekly_rate" = as.numeric(`WEEKLY RATE`)) |>
  mutate(Age = factor(`AGE CATEGORY`, levels = c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr"))) |>
  ggplot(aes(`MMWR-WEEK`, Weekly_rate, group = Age)) + 
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly Hospitalizations per 100,000 People") +
  ggtitle("EIP Site Flu Hospitalization Rates 2022 Season") +
  facet_wrap(~CATCHMENT) 
  

ggsave("Flu_EIP_age_all_sites.pdf", flu_by_age_eip, device = 'pdf', width =15, height = 10)


#make the same plot with finer age stratification
flu_by_age_eip_finer <- flu_eip |> 
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-29 yr", "30-39 yr", 
                               "40-49 yr", "50-64 yr","65+ yr")) |>
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  mutate("Weekly_rate" = as.numeric(`WEEKLY RATE`)) |>
  mutate(Age = factor(`AGE CATEGORY`, levels = c("0-4 yr", "5-17 yr", "18-29 yr", "30-39 yr", 
                                                 "40-49 yr", "50-64 yr","65+ yr"))) |>
  ggplot(aes(`MMWR-WEEK`, Weekly_rate, group = Age)) + 
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly Hospitalizations per 100,000 People") +
  ggtitle("EIP Site Flu Hospitalization Rates 2022 Season") +
  facet_wrap(~CATCHMENT) 
ggsave("Flu_EIP_age_finer_all_sites.pdf", flu_by_age_eip_finer, device = 'pdf', width =15, height = 10)

flu_by_age_eip_finest <- flu_eip |> 
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-29 yr", "30-39 yr", 
                               "40-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+")) |>
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  mutate("Weekly_rate" = as.numeric(`WEEKLY RATE`)) |>
  mutate(Age = factor(`AGE CATEGORY`, levels = c("0-4 yr", "5-17 yr", "18-29 yr", "30-39 yr", 
                                                 "40-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+"))) |>
  ggplot(aes(`MMWR-WEEK`, Weekly_rate, group = Age)) + 
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly Hospitalizations per 100,000 People") +
  ggtitle("EIP Site Flu Hospitalization Rates 2022 Season") +
  facet_wrap(~CATCHMENT) 

ggsave("Flu_EIP_age_finest_all_sites.pdf", flu_by_age_eip_finest, device = 'pdf', width =15, height = 10)

flu_by_age_eip_midcompress <- flu_eip |> 
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+")) |>
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  mutate("Weekly_rate" = as.numeric(`WEEKLY RATE`)) |>
  mutate(Age = factor(`AGE CATEGORY`, levels = c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+"))) |>
  ggplot(aes(`MMWR-WEEK`, Weekly_rate, group = Age)) + 
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly Hospitalizations per 100,000 People") +
  ggtitle("EIP Site Flu Hospitalization Rates 2022 Season") +
  facet_wrap(~CATCHMENT) 
ggsave("Flu_EIP_age_midcompress_all_sites.pdf", flu_by_age_eip_midcompress, device = 'pdf', width =15, height = 10)


#do the same thing for RSV
RSV_age_strat_1 <- rsv |>
  filter(Season == "2022-2023") |>
  filter(Sex == "Overall") |>
  filter(Race == "Overall") |>
  filter(`Age Category` %in% c("0-4 years", "5-17 years", "18-49 years", "50-64 years","65+ years")) |>
  mutate(Age = factor(`Age Category`, levels = c("0-4 years", "5-17 years", "18-49 years", "50-64 years","65+ years"))) |>
  filter(`MMWR Week` %in% c(40:52)) |>
  mutate(Weekly_rate = as.numeric(Rate),
         MMWR_Week = as.numeric(`MMWR Week`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Age)) +
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly RSV Hospitalizations per 100,000 People") +
  ggtitle("EIP Site RSV Hospitalization Rates 2022 Season") +
  facet_wrap(~State) 
ggsave("RSV_EIP_age_strat_1_all_sites.pdf", RSV_age_strat_1, device = 'pdf', width= 15, height = 10)  


rsv_age_strat2 <- rsv |>
  filter(Season == "2022-2023") |>
  filter(Sex == "Overall") |>
  filter(Race == "Overall") |>
  filter(`Age Category` %in% c("----0-<6 months", "----6-<12 months", "----1-<2 years",
                               "----2-4 years","----5-11 years", "----12-17 years",
                               "----18-29 years", "----30-39 years", "----40-49 years", "----50-64 years",
                               "----65-74 years",  "----75-84 years", "----85+ years" )) |>
  mutate(Age =factor(`Age Category`, levels = c("----0-<6 months", "----6-<12 months", "----1-<2 years",
                                                "----2-4 years","----5-11 years", "----12-17 years",
                                                "----18-29 years", "----30-39 years", "----40-49 years", "----50-64 years",
                                                "----65-74 years",  "----75-84 years", "----85+ years" ))) |>
  filter(`MMWR Week` %in% c(40:52)) |>
  mutate(Weekly_rate = as.numeric(Rate),
         MMWR_Week = as.numeric(`MMWR Week`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Age)) +
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly RSV Hospitalizations per 100,000 People") +
  ggtitle("EIP Site RSV Hospitalization Rates 2022 Season") +
  facet_wrap(~State) 

ggsave("RSV_EIP_age_strat_2_all_sites.pdf", rsv_age_strat2, device = 'pdf', width= 15, height = 10) 

rsv_age_strat3 <- rsv |>
  filter(Season == "2022-2023") |>
  filter(Sex == "Overall") |>
  filter(Race == "Overall") |>
  filter(`Age Category` %in% c("0-4 years", "5-17 years", "18-49 years", "50-64 years",  "----75-84 years", "----85+ years" )) |>
  mutate(Age =factor(`Age Category`, levels = c("0-4 years", "5-17 years", "18-49 years", "50-64 years",  "----75-84 years", "----85+ years" ))) |>
  filter(`MMWR Week` %in% c(40:52)) |>
  mutate(Weekly_rate = as.numeric(Rate),
         MMWR_Week = as.numeric(`MMWR Week`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Age)) +
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly RSV Hospitalizations per 100,000 People") +
  ggtitle("EIP Site RSV Hospitalization Rates 2022 Season") +
  facet_wrap(~State) 

ggsave("RSV_EIP_age_strat_3_all_sites.pdf", rsv_age_strat3, device = 'pdf', width= 15, height = 10) 


#now try with COVID; will just need to average over all the states with surveillance
EIP_states <- unique(flu_eip$CATCHMENT)
covid_filt <- covid |> 
  filter(state %in% EIP_states)
covid_EIP_states <- covid_filt |> 
  group_by(state, date_updated) |>
  summarize(average_hosp_rate = mean(covid_hospital_admissions_per_100k, na.rm = TRUE))

#trying to align with RSV data, seems like week 40 would be 22-10-06, and from there they would go to week 52
#use lubridate to make dates

# covid_EIP_states <- covid_EIP_states |>
#   mutate(week_dates = date(date_updated)) |>
#   filter(week_dates >= "2022-10-06") |>
#   select(-date_updated)
# unique(covid_EIP_states$week_dates)
# #number the weeks 
# week_num <- data.frame(week_dates = unique(covid_EIP_states$week_dates), weeks = 40:53)
# covid_EIP_states <- left_join(covid_EIP_states, week_num)
#this data doesn't have age group

#try from covid-net data
covid_net <- read_csv("COVID-19Surveillance_All_Data copy.csv")
#covid-19 associated hospitalization rates per 100,000 population
#do same as with other 2 viruses
COVID_age_strat_1 <- covid_net |>
  filter(`MMWR-YEAR` == 2022) |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  filter(SEX == "Overall") |>
  filter(RACE == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65-74 yr",
                               "75-84 yr", "85+")) |>
  mutate(Age = factor(`AGE CATEGORY`, levels = c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65-74 yr",
                                                 "75-84 yr", "85+"))) |>
  mutate(Weekly_rate = as.numeric(`WEEKLY RATE`),
         MMWR_Week = as.numeric(`MMWR-WEEK`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Age)) +
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly COVID-19 Hospitalizations per 100,000 People") +
  ggtitle("EIP Site COVID-19 Hospitalization Rates 2022 Season") +
  facet_wrap(~CATCHMENT) 
ggsave("COVID_age_strat_1.pdf", COVID_age_strat_1, device = "pdf", width = 15, height = 10)

COVID_age_strat_2 <- covid_net |>
  filter(`MMWR-YEAR` == 2022) |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  filter(SEX == "Overall") |>
  filter(RACE == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-29 yr", "30-39 yr", 
                               "40-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+")) |>
  mutate(Age = factor(`AGE CATEGORY`, levels = c("0-4 yr", "5-17 yr", "18-29 yr", "30-39 yr", 
                                                 "40-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+"))) |>
  mutate(Weekly_rate = as.numeric(`WEEKLY RATE`),
         MMWR_Week = as.numeric(`MMWR-WEEK`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Age)) +
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly COVID-19 Hospitalizations per 100,000 People") +
  ggtitle("EIP Site COVID-19 Hospitalization Rates 2022 Season") +
  facet_wrap(~CATCHMENT)
ggsave("COVID_age_strat_2.pdf", COVID_age_strat_2, device = "pdf", width = 15, height = 10)



COVID_age_strat_3 <- covid_net |>
  filter(`MMWR-YEAR` == 2022) |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  filter(SEX == "Overall") |>
  filter(RACE == "Overall") |>
  filter(`AGE CATEGORY` %in%  c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+")) |>
  mutate(Age = factor(`AGE CATEGORY`, levels = c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65-74 yr", "75-84 yr", "85+"))) |>
  mutate(Weekly_rate = as.numeric(`WEEKLY RATE`),
         MMWR_Week = as.numeric(`MMWR-WEEK`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Age)) +
  geom_line(aes(color = Age)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Age Category") + 
  ylab("Weekly COVID-19 Hospitalizations per 100,000 People") +
  ggtitle("EIP Site COVID-19 Hospitalization Rates 2022 Season") +
  facet_wrap(~CATCHMENT)
ggsave("COVID_age_strat_3.pdf", COVID_age_strat_3, device = "pdf", width = 15, height = 10)


#now plot all 3 viruses together for each age group. Start with 0-4
flu_04 <- flu_eip |> 
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(`AGE CATEGORY` == "0-4 yr") |>
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) 
rsv_04 <- rsv |>
  filter(Season == "2022-2023") |>
  filter(Sex == "Overall") |>
  filter(Race == "Overall") |>
  filter(`Age Category`== "0-4 years") |>
  filter(`MMWR Week` %in% c(40:52)) 

covid_04 <- covid_net |>
  filter(`MMWR-YEAR` == 2022) |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  filter(SEX == "Overall") |>
  filter(RACE == "Overall") |>
  filter(`AGE CATEGORY` =="0-4 yr") |>
  select(CATCHMENT, `MMWR-YEAR`, `MMWR-WEEK`, `AGE CATEGORY`, `WEEKLY RATE`)

#I may want to merge them first then  select age groups
flu_to_merge <- flu_eip |> 
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  select(CATCHMENT, `MMWR-WEEK`, `AGE CATEGORY`, `WEEKLY RATE`) |>
  mutate("Infection" = "Influenza")
covid_to_merge <- covid_net |>
  filter(`MMWR-YEAR` == 2022) |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  filter(SEX == "Overall") |>
  filter(RACE == "Overall") |>
  select(CATCHMENT, `MMWR-WEEK`, `AGE CATEGORY`, `WEEKLY RATE`) |>
  mutate("Infection" = "COVID-19")
rsv_to_merge <- rsv |>
  filter(Season == "2022-2023") |>
  filter(Sex == "Overall") |>
  filter(Race == "Overall") |>
  filter(`MMWR Week` %in% c(40:52)) |>
  select(State, `MMWR Week`, `Age Category`, Rate) |>
  mutate(CATCHMENT = State,
         "MMWR-WEEK" = `MMWR Week`,
         "AGE CATEGORY" = `Age Category`,
         "WEEKLY RATE" = Rate) |>
  select(CATCHMENT,`MMWR-WEEK`, `AGE CATEGORY`,`WEEKLY RATE`) |>
  mutate("Infection"= "RSV")
threevir_allages <- rbind(flu_to_merge, covid_to_merge, rsv_to_merge)

#try plotting across all areas all 3 viruses for 0-4 age group
threevir_0four<- threevir_allages |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "0-4 years")) |>
  filter(CATCHMENT %in% c("California", "Colorado", "Connecticut", "Georgia",
                          "Maryland", "Minnesota", "New Mexico", "Oregon", "Tennessee")) |>
  mutate(Weekly_rate = as.numeric(`WEEKLY RATE`),
         MMWR_Week = as.numeric(`MMWR-WEEK`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Infection)) +
  geom_line(aes(color = Infection)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Infection") + 
  ylab("Weekly Hospitalizations per 100,000 People") +
  ggtitle("EIP Site Hospitalization Rates 2022 Season \n 0-4 Age Group") +
  facet_wrap(~CATCHMENT)

#ggsave("three_viruses_0-4_EIP_sites.pdf", threevir_0four, device = "pdf", width = 15, height = 10)

threevir_fiveele<- threevir_allages |>
  filter(`AGE CATEGORY` %in% c("5-11  yr", "----5-11 years")) |>
  filter(CATCHMENT %in% c("California", "Colorado", "Connecticut", "Georgia",
                          "Maryland", "Minnesota", "New Mexico", "Oregon", "Tennessee")) |>
  mutate(Weekly_rate = as.numeric(`WEEKLY RATE`),
         MMWR_Week = as.numeric(`MMWR-WEEK`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Infection)) +
  geom_line(aes(color = Infection)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Infection") + 
  ylab("Weekly Hospitalizations per 100,000 People") +
  ggtitle("EIP Site Hospitalization Rates 2022 Season \n 5-11 Age Group") +
  facet_wrap(~CATCHMENT)
ggsave("three_viruses_5-11_EIP_sites.pdf", threevir_fiveele, device = "pdf", width = 15, height = 10)
#this one is the most interesting so far!
threevir_twelsev<- threevir_allages |>
  filter(`AGE CATEGORY` %in% c("12-17 yr" ,"----12-17 years")) |>
  filter(CATCHMENT %in% c("California", "Colorado", "Connecticut", "Georgia",
                          "Maryland", "Minnesota", "New Mexico", "Oregon", "Tennessee")) |>
  mutate(Weekly_rate = as.numeric(`WEEKLY RATE`),
         MMWR_Week = as.numeric(`MMWR-WEEK`)) |>
  ggplot(aes(MMWR_Week, Weekly_rate, group = Infection)) +
  geom_line(aes(color = Infection)) +
  scale_x_continuous("2022 MMWR Week", breaks = c(40:52), limits = c(40,52)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() + 
  labs(color = "Infection") + 
  ylab("Weekly Hospitalizations per 100,000 People") +
  ggtitle("EIP Site Hospitalization Rates 2022 Season \n 12-17 Age Group") +
  facet_wrap(~CATCHMENT)

ggsave("three_viruses_12-17_EIP_sites.pdf", threevir_twelsev, device = "pdf", width = 15, height = 10)



