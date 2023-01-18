#ODE Solver 
library(deSolve)
library(ggplot2)
library(patchwork)

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

times <- seq(0, 1000, by = 0.01)

# Scenario 1: with viral interference, both start at same time (Recovery period = 2 days; eta_1 = eta_2 = 0.5)
# Parameters
# Initial values
state <- c(S=10000, E_1=0, E_2=0, I_1=50, I_2=50, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state)
parameters <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
                   delta_1=0.07142857, delta_2=0.07142857, eta_1=0.5, eta_2=0.5)


# Model outputs
output<- ode(y = state, times = times, func = Two_Infection, parms = parameters)
df1 <- as.data.frame(output)
df1 <- df1 %>%
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

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
p1 <- plot_I(data = df1) + 
  labs(subtitle = "Assuming 2 days of recovery period assuming no offset")

# Scenario 2: no viral interference (Recovery period = 0 days; eta_1 = eta_2 = 1000)

state_2 <- c(S=10000, E_1=0, E_2=0, I_1=50, I_2=50, R_1=0, R_2=0, 
           S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state_2)
par2 <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
                   delta_1=0.07142857, delta_2=0.07142857, eta_1 = 1000, eta_2=1000)

output_2 <- ode(y = state, times = times, func = Two_Infection, parms = par2)
df2 <- as.data.frame(output_2)
df2 <- df2 %>%
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

## Plot scenario 2
p2 <- plot_I(data = df2) + 
  labs(subtitle = "Assuming 0 days of recovery period and no offset")

plots_to_save <- p1 + p2
ggsave("Figures/Dynamics_hypothetical.png", plot = plots_to_save, width=11, height=4, dpi=300)#ODE Solver 

# Scenario 3: No offset Extended viral interference (Recovery period = 5 days; eta_1 = eta_2 = .2)

state_3 <- c(S=10000, E_1=0, E_2=0, I_1=50, I_2=50, R_1=0, R_2=0, 
             S_12=0, S_21=0, E_12=0, E_21=0, I_12=0, I_21=0, R=0)

N =  sum(state_3)
par3 <- list(beta_1=3.0*0.07142857/N, beta_2=2.7*0.07142857/N, kappa_1=0.2, kappa_2=0.1481481, 
             delta_1=0.07142857, delta_2=0.07142857, eta_1 = .2, eta_2=.2)

output_3 <- ode(y = state, times = times, func = Two_Infection, parms = par3)
df3 <- as.data.frame(output_2)
df3 <- df3 %>%
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

## Plot scenario 3
p3 <- plot_I(data = df3) + 
  labs(subtitle = "Assuming 5 days of recovery period and no offset")

plots_to_save <- p1 + p2 + p3
ggsave("Figures/Dynamics_hypothetical.png", plot = plots_to_save, width=11, height=4, dpi=300)#ODE Solver 


