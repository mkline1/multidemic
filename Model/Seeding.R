#----------------------------------- 
duration <- 1000
delta <- 0.1
out <- matrix(NA, nrow = duration + 1, ncol = length(state) + 1)

tau_1 <- 0
tau_2 <- 14

update_conditions <- function(t){
  
  if(!is.na(tau_1) & tau_1 == t){
    
    state[c("S","I_1")] <- state[c("S","I_1")] + c(-1, 1)
    
  }
  
  if(!is.na(tau_2) & tau_2 == t){
    
    state[c("S","I_2")] <- state[c("S","I_2")] + c(-1, 1) 
    
  } 
  
}

#---------------------------------

update_conditions(t = 0); out[1, 2:ncol(out)] <- state

for(i in 1:duration){
  
  t <- seq(i - 1, i, by = delta)
  
  update_conditions(t = i - 1)
  
  output <- ode(y = state, 
                times = times, 
                func = Two_Infection, 
                parms = parameters, 
                method = "adams"); rm(t)
  
  out[i+1,] <- output[nrow(output),]
  
  state <- output[nrow(output), 2:ncol(output)]; rm(output)
  
}

df1 <- as.data.frame(out)

colnames(df1) <- c("time","S", "E", "E_2", "I_1", "I_2", "R_1", "R_2", 
                   "S_12", "S_21", "E_12", "E_21", "I_12", "I_21")

df1 <- df1 %>%
  mutate(total_I_1 = I_1 + I_21,
         total_I_2 = I_2 + I_12)

df1 <- cbind(time=0:1000, df1)
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
    lims(x=c(0,100),
         y=c(0,5e8)) +
    theme_minimal() +
    labs(title = "Dynamics of the number infected over time",
         y = "Number infected",
         x = "Time(Day)")
}

## Plot scenario 1
p1 <- plot_I(data = df1) + 
  labs(subtitle = "Assuming 2 days of recovery period (with protection from viral interference)")

