library(tidyverse)
library(energy)
library(viridis)

setwd("~/Desktop/research/multidemic/CDC data/")

flu <- read_csv("data/FluSurveillance_Custom_Download_Data_header_removed.csv")
rsv <- read_csv("data/Weekly_Rates_of_Laboratory-Confirmed_RSV_Hospitalizations_from_the_RSV-NET_Surveillance_System (1).csv")
cov <- read_csv("data/COVID-19Surveillance_All_Data copy.csv")

# flu
flu <- flu |> 
  drop_na() |> 
  
  filter(NETWORK == "EIP") |> 
  
  filter(`SEX CATEGORY` == "Overall") |>
  filter(`RACE CATEGORY` == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr")) |>
  filter(YEAR == "2022-23") |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  mutate(`WEEKLY RATE` = as.numeric(`WEEKLY RATE`)) |>
  mutate(`AGE CATEGORY` = factor(`AGE CATEGORY`, 
                      levels = c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr"))) |> 
  select(CATCHMENT, `MMWR-WEEK`, `AGE CATEGORY`, `WEEKLY RATE`) |>
  mutate(Infection = "Influenza")

names(flu)

# rsv
rsv <- rsv |>
  filter(Sex == "Overall") |>
  filter(Race == "Overall") |>
  filter(`Age Category` %in% c("0-4 years", "5-17 years", "18-49 years", "50-64 years","65+ years")) |>
  filter(Season == "2022-2023") |>
  mutate(`AGE CATEGORY` = factor(`Age Category`, levels = c("0-4 years", "5-17 years", "18-49 years", "50-64 years","65+ years"))) |>
  filter(`MMWR Week` %in% c(40:52)) |>
  mutate(`WEEKLY RATE` = as.numeric(Rate),
         `MMWR-WEEK` = as.numeric(`MMWR Week`),
         CATCHMENT = State) |> 
  select(CATCHMENT,`MMWR-WEEK`,`AGE CATEGORY`,`WEEKLY RATE`) |> 
  mutate(Infection = "RSV")
names(rsv)

cov <- cov |>
  filter(`MMWR-YEAR` == 2022) |>
  filter(`MMWR-WEEK` %in% c(40:52)) |>
  filter(SEX == "Overall") |>
  filter(RACE == "Overall") |>
  filter(`AGE CATEGORY` %in% c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr")) |>
  mutate(`AGE CATEGORY` = factor(`AGE CATEGORY`,
                                 levels = c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr"))) |>
  select(CATCHMENT, `MMWR-WEEK`, `AGE CATEGORY`, `WEEKLY RATE`) |>
  mutate(Infection = "COVID-19")

# common catchments
comm_catch <- Reduce(intersect,list(unique(flu$CATCHMENT),unique(rsv$CATCHMENT),unique(cov$CATCHMENT)))

flu <- flu |> filter(CATCHMENT %in% comm_catch)
rsv <- rsv |> filter(CATCHMENT %in% comm_catch)
cov <- cov |> filter(CATCHMENT %in% comm_catch)

# minimum data for all 3 viruses
cov <- cov |> filter(`MMWR-WEEK` %in% c(40:51))
flu <- flu |> filter(`MMWR-WEEK` %in% c(40:51))

# match naming
rsv <- rsv |> mutate(`AGE CATEGORY` = case_when(
  `AGE CATEGORY` == "0-4 years" ~ "0-4 yr",
  `AGE CATEGORY` == "5-17 years" ~ "5-17 yr",
  `AGE CATEGORY` == "18-49 years" ~ "18-49 yr",
  `AGE CATEGORY` == "50-64 years" ~ "50-64 yr",
  `AGE CATEGORY` == "65+ years" ~ "65+ yr")) |> 
  mutate(`AGE CATEGORY` = factor(`AGE CATEGORY`,
                                 levels = c("0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr","65+ yr")))

# check
unique(flu$`AGE CATEGORY`)
unique(rsv$`AGE CATEGORY`)
unique(cov$`AGE CATEGORY`)

flu_weekly <- flu |> 
  group_by(`MMWR-WEEK`) |> 
  summarise(comb_rate = sum(`WEEKLY RATE`))

rsv_weekly <- rsv |> 
  group_by(`MMWR-WEEK`) |> 
  summarise(comb_rate = sum(`WEEKLY RATE`))

cov_weekly <- cov |> 
  group_by(`MMWR-WEEK`) |> 
  summarise(comb_rate = sum(`WEEKLY RATE`))

cors_by_age <- function(v1,v2) {
  d1 <- v1
  d2 <- v2
  
  i1 <- d1$Infection[1]
  i2 <- d2$Infection[1]
  
  res <- c()
  
  for (i in comm_catch) {
    
    for (j in levels(d1$`AGE CATEGORY`)) {
      d1_strat <- d1 |> filter(`AGE CATEGORY` == j, CATCHMENT == i)
      d2_strat <- d2 |> filter(`AGE CATEGORY` == j, CATCHMENT == i)
      
      dcor_coef <- dcor(d1_strat$`WEEKLY RATE`,d2_strat$`WEEKLY RATE`)
      dcor_T_p <- unname(dcorT.test(d1_strat$`WEEKLY RATE`,d2_strat$`WEEKLY RATE`)$p.value)
      
      res <- rbind(res,c(i,j,i1,i2,dcor_coef,dcor_T_p))
      
    }
    
  }
  
  colnames(res) <- c("Catchment", "Age", "Vir1", "Vir2", "DcorCoef", "Dcor_p")
  res <- as_tibble(res)
  res
  
}

flu_rsv_res <- cors_by_age(flu,rsv)
flu_cov_res <- cors_by_age(flu,cov)
cov_rsv_res <- cors_by_age(cov,rsv)

heatmap_plot <- bind_rows(flu_rsv_res,flu_cov_res,cov_rsv_res) |> 
  mutate(VirPair = paste0(Vir1,"-",Vir2),
         DcorCoef = as.numeric(DcorCoef),
         Dcor_p = as.numeric(Dcor_p)) |> 
  select(VirPair,Catchment,Age,DcorCoef,Dcor_p) |> 
  ggplot(aes(x=Age,y=Catchment)) +
  geom_tile(aes(fill=DcorCoef)) +
  facet_grid(~ VirPair) +
  scale_fill_viridis(option = "inferno")

heatmap_plot





