# Summary stats

#ODE Solver 
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(viridis)
library(ggpubr)

source("Seeding.R")

# Model
Two_Infection <- function(t, state, parameters){
  
  with(as.list(c(state,parameters)),{
    #rate of change
    dS  = -beta_1*S*(I_1 + I_21) - beta_2*S*(I_2 + I_12)
    dE_1 = beta_1*S*(I_1 + I_21) - kappa_1*E_1 
    dE_2 = beta_2*S*(I_2 + I_12) - kappa_1*E_2
    dI_1 = kappa_1*E_1 - delta_1*I_1
    dI_2 = kappa_2*E_2 - delta_2*I_2
    dR_1 = delta_1*I_1 - eta_1*R_1 
    dR_2 = delta_2*I_2 - eta_2*R_2
    dS_12 = -beta_2*S_12*(I_2 + I_12) + eta_1*R_1
    dS_21 = -beta_1*S_21*(I_1 + I_21) + eta_2*R_2
    dE_12 = beta_2*S_12*(I_2 + I_12) - kappa_2*E_12
    dE_21 = beta_1*S_21*(I_1 + I_21) - kappa_1*E_21
    dI_12 = kappa_2*E_12 - delta_2*I_12
    dI_21 = kappa_1*E_21 - delta_1*I_21
    dR = delta_1*I_21 + delta_2*I_12
    
    return(list(c(dS, dE_1, dE_2, dI_1, dI_2, dR_1, dR_2,
                  dS_12, dS_21, dE_12, dE_21, dI_12, dI_21, dR)))
    
  })
}

times <- seq(0, 500, by = 0.01)

plot_I <- function(data){
  ggplot() +
    geom_line(data=data, aes(x=time, y=I_1, col="RSV (naive)")) +
    geom_line(data=data, aes(x=time, y=I_21, col="RSV (recovered from COVID)")) + 
    geom_line(data=data, aes(x=time, y=I_2, col="COVID (naive)")) +
    geom_line(data=data, aes(x=time, y=I_12, col="COVID (recovered from RSV)")) +
    scale_color_manual(name = "Infected", 
                       values = c("RSV (naive)" = "darkblue", 
                                  "RSV (recovered from COVID)" = "lightblue",
                                  "COVID (naive)" = "darkred", 
                                  "COVID (recovered from RSV)" = "pink")) +
    theme_minimal() +
    labs(title = "Dynamics of the number infected over time",
         y = "Number infected",
         x = "Time(Day)")
}

plot_Itot <- function(data){
  ggplot() +
    geom_line(data=data, aes(x=time, y=total_I_1, col="RSV (total)")) +
    geom_line(data=data, aes(x=time, y=total_I_2, col="COVID (total)")) +
    scale_color_manual(name = "Infected", 
                       values = c("RSV (total)" = "darkblue", 
                                  "COVID (total)" = "darkred")) +
    theme_minimal() +
    labs(title = "Dynamics of the number infected over time",
         y = "Number infected",
         x = "Time(Day)")
}

# Scenario A (null) -------------------------------------------------------

state <- c(S=9900, E_1=0, E_2=0, I_1=50, I_2=50, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state)
parameters <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
                   delta_1=0.07142857, delta_2=0.07142857, eta_1=1000, eta_2=1000)

# Model outputs
output_A<- ode(y = state, times = times, func = Two_Infection, parms = parameters)
dfA <- as.data.frame(output_A)
dfA <- dfA %>%
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

pA <- plot_I(data = dfA) + lims(x=c(0,300)) + labs(subtitle = "Null")
pAtot <- plot_Itot(data = dfA) + lims(x=c(0,300)) + labs(subtitle = "Null")

# Scenario B - 5 day refractory, no offset --------------------------------

state <- c(S=9900, E_1=0, E_2=0, I_1=50, I_2=50, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state)
parameters <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
                   delta_1=0.07142857, delta_2=0.07142857, eta_1=0.2, eta_2=0.2)

# Model outputs
output_B<- ode(y = state, times = times, func = Two_Infection, parms = parameters)
dfB <- as.data.frame(output_B)
dfB <- dfB %>%
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

pB <- plot_I(data = dfB) + lims(x=c(0,300)) + labs(subtitle = "5d refractory")
pBtot <- plot_Itot(data = dfB) + lims(x=c(0,300)) + labs(subtitle = "5d refractory")

# Scenario C - 5 day refractory, COV 14 day offset ------------------------

state <- c(S=10000, E_1=0, E_2=0, I_1=0, I_2=0, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state)
parameters <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
                   delta_1=0.07142857, delta_2=0.07142857, eta_1=0.2, eta_2=0.2)

# Set up
duration <- 500
delta <- 0.01

# An empty out matrix to hold the outputs
out <- matrix(NA, nrow = duration + 1, ncol = length(state) + 1)

# Change the seeding time (in days)
tau_2 <- 14

# Initialize the model
state[c("S","I_1")] <- state[c("S","I_1")] + c(-50, 50)

# Model outputs
out[1, 2:ncol(out)] <- state

for(i in 1:duration){
  
  ## If pathogen 2 starts at time i, add 50 I_2 
  if(tau_2 == i){
    
    state[c("S","I_2")] <- state[c("S","I_2")] + c(-50, 50) 
    
  } 
  
  times <- seq(i - 1, i, by = delta)
  
  output <- ode(y = state, 
                times = times, 
                func = Two_Infection, 
                parms = parameters,
                method = "adams")
  
  out[i+1,] <- output[nrow(output),]
  
  state <- output[nrow(output), 2:ncol(output)]
  
}

dfC <- as.data.frame(out)
colnames(dfC) <- c("time","S", "E_1", "E_2", "I_1", "I_2", "R_1", "R_2", 
                   "S_12", "S_21", "E_12", "E_21", "I_12", "I_21", "R") ###### was missing "R"
dfC[1,1] <- 0

dfC <- dfC |>
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

pC <- plot_I(data = dfC) + lims(x=c(0,300)) + labs(subtitle = "5d refractory, COV 14d offset")
pCtot <- plot_Itot(data = dfC) + lims(x=c(0,300)) + labs(subtitle = "5d refractory, COV 14d offset")

# Scenario D - 5 day refractory, RSV 14 day offset ------------------------

state <- c(S=10000, E_1=0, E_2=0, I_1=0, I_2=0, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state)
parameters <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
                   delta_1=0.07142857, delta_2=0.07142857, eta_1=0.2, eta_2=0.2)

# Set up
duration <- 500
delta <- 0.01

# An empty out matrix to hold the outputs
out <- matrix(NA, nrow = duration + 1, ncol = length(state) + 1)

# Change the seeding time (in days)
tau_2 <- 14

# Initialize the model
state[c("S","I_2")] <- state[c("S","I_2")] + c(-50, 50)

# Model outputs
out[1, 2:ncol(out)] <- state

for(i in 1:duration){
  
  ## If pathogen 2 starts at time i, add 50 I_2 
  if(tau_2 == i){
    
    state[c("S","I_1")] <- state[c("S","I_1")] + c(-50, 50) 
    
  } 
  
  times <- seq(i - 1, i, by = delta)
  
  output <- ode(y = state, 
                times = times, 
                func = Two_Infection, 
                parms = parameters,
                method = "adams")
  
  out[i+1,] <- output[nrow(output),]
  
  state <- output[nrow(output), 2:ncol(output)]
  
}

dfD <- as.data.frame(out)
colnames(dfD) <- c("time","S", "E_1", "E_2", "I_1", "I_2", "R_1", "R_2", 
                   "S_12", "S_21", "E_12", "E_21", "I_12", "I_21", "R") ###### was missing "R"
dfD[1,1] <- 0

dfD <- dfD |>
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

pD <- plot_I(data = dfD) + lims(x=c(0,300)) + labs(subtitle = "5d refractory, RSV 14d offset")
pDtot <- plot_Itot(data = dfD) + lims(x=c(0,300)) + labs(subtitle = "5d refractory, RSV 14d offset")


# Plots -------------------------------------------------------------------

StratPlots <- ggarrange(pA,pB,pC,pD,ncol = 2, nrow =2)
TotPlots <- ggarrange(pAtot,pBtot,pCtot,pDtot,ncol = 2, nrow =2)

ggsave("../Figures/StratPlots.jpg",plot = StratPlots,width=11, height=7, dpi=300,bg = "white")
ggsave("../Figures/TotPlots.jpg",plot = TotPlots,width=11, height=7, dpi=300,bg = "white")

# Summary Stats -----------------------------------------------------------

epi_summary_stats <- function(x) {
  x |> 
    summarise(
      
      # seed times
      rsv_seed_t = min(time[which(I_1>0)]),
      cov_seed_t = min(time[which(I_2>0)]),
      
      # RSV naive
      rsv_naive_t_to_peak = time[which.max(I_1)] - min(time[which(I_1>0)]),
      rsv_naive_peak_inf = max(I_1),
      
      # COV naive
      cov_naive_t_to_peak = time[which.max(I_2)] - min(time[which(I_2>0)]),
      cov_naive_peak_inf = max(I_2),
      
      # RSV recovered
      rsv_rec_t_to_peak = time[which.max(I_12)] - min(time[which(I_1>0)]),
      rsv_rec_peak_inf = max(I_12),
      
      # COV recovered
      cov_rec_t_to_peak = time[which.max(I_21)] - min(time[which(I_2>0)]),
      cov_rec_peak_inf = max(I_21),
      
      # RSV total
      rsv_tot_t_to_peak = time[which.max(total_I_1)] - min(time[which(I_1>0)]),
      rsv_tot_peak_inf = max(total_I_1),
      
      # COV total
      cov_tot_t_to_peak = time[which.max(total_I_2)] - min(time[which(I_2>0)]),
      cov_tot_peak_inf = max(total_I_2)
    )
}

dfA_summary <- epi_summary_stats(dfA) |> mutate(scenario = "R0_O0")
dfB_summary <- epi_summary_stats(dfB) |> mutate(scenario = "R5_O0")
dfC_COVoff_summary <- epi_summary_stats(dfC) |> mutate(scenario = "R5_cov14")
dfD_RSVoff_summary <- epi_summary_stats(dfD) |> mutate(scenario = "R5_rsv14")

all_summary <- rbind(dfA_summary,dfB_summary,dfC_COVoff_summary,dfD_RSVoff_summary) |> 
  as_tibble()

names(all_summary)

SummaryPlots <- all_summary |> 
  select(scenario, contains("tot_peak")) |> 
  rename("COVID (total)" = cov_tot_peak_inf,
         "RSV (total)" = rsv_tot_peak_inf) |> 
  pivot_longer(!scenario, names_to = "Virus", values_to = "infections") |> 
  ggplot(aes(x=scenario,y=infections,fill=Virus)) +
  geom_col() +
  scale_fill_manual(values = c("darkred","darkblue")) +
  facet_grid(~ Virus) +
  theme_classic() +
  ggtitle("Summary Statistics by Scenario")

ggsave("../Figures/SummaryPlots.png",plot = SummaryPlots,width=7, height=7, dpi=300,bg = "white")







