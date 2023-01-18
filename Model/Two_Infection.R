#ODE Solver 
library(deSolve)
library(ggplot2)
library(patchwork)

source("Model/Seeding.R")
 
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

 # Scenario 1: with viral interference (Recovery period = 10 days; eta_1 = eta_2 = 0.1)
state <- c(S=10000, E_1=0, E_2=0, I_1=0, I_2=0, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state)
parameters <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
                   delta_1=0.07142857, delta_2=0.07142857, eta_1=0.5, eta_2=0.5)

  # Set up
duration <- 1000
delta <- 0.1

  # An empty out matrix to hold the outputs
out <- matrix(NA, nrow = duration + 1, ncol = length(state) + 1)

  # Change the seeding time (in days)
tau_2 <- 14

  # Initialize the model
state[c("S","I_1")] <- state[c("S","I_1")] + c(-1, 1)

  # Model outputs
out[1, 2:ncol(out)] <- state

for(i in 1:duration){
   
   ## If pathogen 2 starts at time i, add one I_2 
  if(tau_2 == i){
    
    state[c("S","I_2")] <- state[c("S","I_2")] + c(-1, 1) 
    
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

df1 <- as.data.frame(out)
colnames(df1) <- c("time","S", "E", "E_2", "I_1", "I_2", "R_1", "R_2", 
                   "S_12", "S_21", "E_12", "E_21", "I_12", "I_21")

 # Plot
  ## Plot function
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

  ## Plot scenario 1
plot_I(data = df1) + 
  lims(x=c(0,250)) +
  labs(subtitle = "Assuming 2 days of recovery period (with protection from viral interference)")
