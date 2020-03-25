# 
######### Forward: Mut -> Phenotype
Forward <- function(total_protein, ddF, ddD, ddB, ddT ){
  
  my_protein = total_protein
  
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  Ka= 5 *10^7 # ka for dimmer formation 
  DNA= 10^(-9)
  
  deltadimer= -R*Temp*log(Ka) 
  Ka = exp((-(deltadimer+ ddD))/(R*Temp))  
  
  wt_deltaG=-1  #-2.908485 #, set as in the original analysis 
  deltaG_folding1= wt_deltaG+ ddF
  
  #configarations=list(c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1), c(1,1,0), c(1,0,1), c(0,1,1), c(1,1,1)) # 8 configurations OR1, OR2, OR3. 
  wt_BdeltaGs<- c(-11.7 + ddB, -10.1+ ddB, -10.1+ ddB) 

  # both on proteins
  
  coop= -2 + ddT  # OR1 and OR2 ; or OR2 and OR3
  
  relations=list(deltaGs=c(0, wt_BdeltaGs[1], wt_BdeltaGs[2], wt_BdeltaGs[3],  wt_BdeltaGs[1]+ wt_BdeltaGs[2]+ coop,  wt_BdeltaGs[1]+ wt_BdeltaGs[3],  wt_BdeltaGs[3]+ wt_BdeltaGs[2]+ coop, sum(wt_BdeltaGs)+ coop),
                 dimer_number=c(0,1,1,1,2,2,2,3)) # thermo_parameters
  
  
  Mytotal_protein <- function(free_dimer) {
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
    abs(log10(my_protein)- log10((free_dimer/Ka)^0.5+ 2*free_dimer+ r_formd+ (free_dimer/Ka)^0.5*exp(deltaG_folding1/(R*Temp))))# monomer + free dimer + bound dimer + unfolded
    # in log scale so that the differerneces are expanded. Ohterwise the differences are too small to keep searching
  }
  
free_dimer=  unlist(optimize(Mytotal_protein, c(10^-40, 10^-6), tol= 10^-23)[1]) # from 10^-21 changed to -40

  configaration<-c() 
  sum_config<- c()
  r_formd<- c() 
  for (i in 1:8) {
    configaration[i]<- exp(-relations[[1]][i]/(R*Temp))*free_dimer^relations[[2]][i]
    sum_config<- sum(configaration)
    config_prob<- configaration/sum_config 
    r_formd<- sum(config_prob[i]*relations[[2]][i])*2*DNA
  }
  
  Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
  Propor_activation=  config_prob[3]+ config_prob[5] # OR2 bound or OR1, OR2 bound
  return(
    list( pr = 1-Propor_sup, prm= Propor_activation)) 
} 

##### Forward 2, which 2. 
###### 
pr_prm_Intramol_protein_withDNA_whichparam <- function(total_protein,whichparam1, param1,whichparam2, param2) { 
  
  my_protein = total_protein
  
  # mydeltas<- c(total_protein1, deltadeltaF1_2, deltadeltadimer3, deltadeltaB1_4, deltadeltacoop5, ddg_or1_6, ddg_or2_7, ddg_or3_8) 
  mydeltas<- c(total_protein, 0, 0, 0, 0, 0, 0, 0)
  if (whichparam1== whichparam2) {
    if (whichparam1>1) {
      mydeltas[whichparam1] = param1+param2 
    } else {
      mydeltas[whichparam1] = param1* param2/ mywtamount
    }
  } else {
    mydeltas[whichparam1] = param1  
    mydeltas[whichparam2] = param2
  }
  
  #1. define parameter values 
  
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  Ka= 5 *10^7 # ka for dimmer formation 
  DNA= 10^(-9)
  
  deltadimer= -R*Temp*log(Ka) 
  Ka = exp((-(deltadimer+ mydeltas[3]))/(R*Temp))  
  
  wt_deltaG=  -1 
  deltaG_folding1= wt_deltaG+ mydeltas[2] 
  
  wt_BdeltaGs<- c(-11.7+mydeltas[4] +mydeltas[6], -10.1+mydeltas[4] +mydeltas[7], -10.1+mydeltas[4] +mydeltas[8] ) 
  
  coop= -2 + mydeltas[5]  # OR1 and OR2 ; or OR2 and OR3
  
  
  #####
  #2. CI-OR thermodynamics model with the updated parameters, empirically search free-dimer vs total concentraiton relationship.  
  
  relations=list(deltaGs=c(0, wt_BdeltaGs[1], wt_BdeltaGs[2], wt_BdeltaGs[3],  wt_BdeltaGs[1]+ wt_BdeltaGs[2]+ coop,  wt_BdeltaGs[1]+ wt_BdeltaGs[3],  wt_BdeltaGs[3]+ wt_BdeltaGs[2]+ coop, sum(wt_BdeltaGs)+ coop),
                 dimer_number=c(0,1,1,1,2,2,2,3)) # thermo_parameters
  
  Mytotal_protein <- function(free_dimer) {
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
    abs(log10(my_protein)- log10((free_dimer/Ka)^0.5+ 2*free_dimer+ r_formd+ (free_dimer/Ka)^0.5*exp(deltaG_folding1/(R*Temp))))# monomer + free dimer + bound dimer + unfolded
    # in log scale so that the differerneces are expanded. Ohterwise the differences are too small to keep searching
  }
  
  free_dimer=  unlist(optimize(Mytotal_protein, c(10^-40, 10^-6), tol= 10^-23)[1]) # from 10^-21 changed to -40
  
  ##
  ### 3.calculate the PR and PRM activity and return 
  configaration<-c() 
  sum_config<- c()
  r_formd<- c() 
  for (i in 1:8) {
    configaration[i]<- exp(-relations[[1]][i]/(R*Temp))*free_dimer^relations[[2]][i]
    sum_config<- sum(configaration)
    config_prob<- configaration/sum_config 
    r_formd<- sum(config_prob[i]*relations[[2]][i])*2*DNA
  }
  
  Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
  Propor_activation=  config_prob[3]+ config_prob[5] # OR2 bound or OR1, OR2 bound
  return(
      list( pr = 1-Propor_sup, prm= Propor_activation)) 
}  


