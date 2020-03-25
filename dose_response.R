from_free_dimer <- function(free_dimer, ddF, ddD, ddB, ddT ){
  
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  Ka= 5 *10^7 # ka for dimmer formation 
  DNA= 10^(-9)
  
  deltadimer= -R*Temp*log(Ka) 
  Ka = exp((-(deltadimer+ ddD))/(R*Temp))  
  
  wt_deltaG= -1 # -2.908485 #, set as in the original analysis 
  deltaG_folding1= wt_deltaG+ ddF
  
  #configarations=list(c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1), c(1,1,0), c(1,0,1), c(0,1,1), c(1,1,1)) # 8 configurations OR1, OR2, OR3. 
  wt_BdeltaGs<- c(-11.7 + ddB, -10.1+ ddB, -10.1+ ddB) 
  
  # both on proteins
  
  coop= -2 + ddT  # OR1 and OR2 ; or OR2 and OR3
  
  relations=list(deltaGs=c(0, wt_BdeltaGs[1], wt_BdeltaGs[2], wt_BdeltaGs[3],  wt_BdeltaGs[1]+ wt_BdeltaGs[2]+ coop,  wt_BdeltaGs[1]+ wt_BdeltaGs[3],  wt_BdeltaGs[3]+ wt_BdeltaGs[2]+ coop, sum(wt_BdeltaGs)+ coop),
                 dimer_number=c(0,1,1,1,2,2,2,3)) # thermo_parameters
  
  
  configaration<-c() 
  
  Total_protein <- c()
  Propor_sup<- c()
  sum_config<- c()
  r_formd<- c()
  
  for (i in 1:8) {
    configaration[i]<- exp(-relations[[1]][i]/(R*Temp))*free_dimer^relations[[2]][i]
    sum_config<- sum(configaration)
    config_prob<- configaration/sum_config 
    r_formd<- sum(config_prob[i]*relations[[2]][i])*2*DNA
  }
  Total_protein= (free_dimer/Ka)^0.5+ 2*free_dimer+ r_formd+ (free_dimer/Ka)^0.5*exp(deltaG_folding1/(R*Temp))# monomer + free dimer + bound dimer + unfolded
  
  Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
  Propor_activation=  config_prob[3]+ config_prob[5] # OR2 bound or OR1, OR2 bound
  return(
    list(Total_protein , 1-Propor_sup, Propor_activation)) 
} 
