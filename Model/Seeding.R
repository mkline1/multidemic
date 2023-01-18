#----------------------------------- 
update_conditions <- function(t){
  
  if(tau_1 == t){
    
    state[c("S","I_1")] <- state[c("S","I_1")] + c(-1, 1)
    
  }
  
  if(tau_2 == t){
    
    state[c("S","I_2")] <- state[c("S","I_2")] + c(-1, 1) 
    
  } 
  
}
