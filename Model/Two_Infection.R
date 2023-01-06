#ODE Solver 

library(deSolve)

parameters <- list(beta_1=0.001, beta_2=0.001, kappa_1=0.3, kappa_2=0.3, 
                   delta_1=0.25, delta_2=0.25, eta_1=0.5, eta_2=0.5)

state <- c(S=10e8, E_1=0, E_2=0, I_1=1, I_2=1, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

Two_Infection <- function(t, state, parameters){
  
  with(as.list(c(state,parameters)),{
    #rate of change
    dS = -beta_1*S*(I_1 + I_21) -beta_2*S*(I_2 + I_12)
    dE_1 = beta_1*S*(I_1 + I_21) - kappa_1*E_1 
    dE_2 = beta_2*S*(I_2 + I_12) - kappa_1*E_2
    dI_1 = kappa_1*E_1 - delta_1*I_1
    dI_2 = kappa_2*E_2 - delta_2*I_2
    dR_1 = delta_1*I_1 -eta_1*R_1 
    dR_2 = delta_2*I_2 -eta_2*R_2
    dS_12 = -beta_2*(I_2+ I_12) + eta_1*R_1
    dS_21 = -beta_1*(I_1 + I_21) + eta_2*R_2
    dE_12 = beta_2*(I_2+ I_12) - kappa_2*E_12
    dE_21 = beta_1*(I_1 + I_21) - kappa_1*E_21
    dI_12 = kappa_2*E_12 - delta_2*I_12
    dI_21 = kappa_1*E_21 - delta_1*I_21
    dR = delta_1*I_21 + delta_2*I_12
    
    return(list(c(dS, dE_1, dE_2, dI_1, dI_2, dR_1, dR_2,
                  dS_12, dS_21, dE_12, dE_21, dI_12, dI_21, dR)))
    
  })
}

times <- seq(0, 1000, by = 0.01)
output<- ode(y = state, times = times, func = Two_Infection, parms = parameters)


