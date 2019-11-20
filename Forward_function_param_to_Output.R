#####
##### Here are two functions: 
#  " pr_prm_Intramol_protein_withDNA_four_param " & "pr_prm_Intramol_protein_withDNA_whichparam"   
# The two functions are basically the same, but in two different formats for inputting values. 
# First one, with defined orders for 4 parameter inputs (easy for single mutants' effects, with any combiantion of parameter changes) 
# Second one, by inputting which paraemter is changed with what value, is easy for generating lists with defined combination of parameter changes. 

 pr_prm_Intramol_protein_withDNA_four_param <- function(mywtamount, param1,param2, param3, param4) { 
  # in the order of ddG_F, ddG_D, ddG_B , ddG_T
  mydeltas<- c(mywtamount, param1, param2, param3, param4, 0, 0, 0) 
  
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
  
  configaration<-c() 
  free_dimer<- c() # free dimer ones 
  for (i in 2:100) { 
    free_dimer[1]<- 10^(-40) 
    free_dimer[i]<- free_dimer[i-1]*2.4
  }  
  Total_protein <- c()
  Propor_sup<- c()
  sum_config<- c()
  r_formd<- c()
  
  library(rootSolve)
  
  for (j in 1:100) { 
    for (i in 1:8) {
      configaration[i]<- exp(-relations[[1]][i]/(R*Temp))*free_dimer[j]^relations[[2]][i]
      sum_config<- sum(configaration)
      config_prob<- configaration/sum_config 
      r_formd[j]<- sum(config_prob[i]*relations[[2]][i])*2*DNA
    }
    Total_protein[j]= (free_dimer[j]/Ka)^0.5+ 2*free_dimer[j]+ r_formd[j]
  }
  wt_trial<- data.frame(Total_protein, free_dimer)
  loess_linear<- loess(wt_trial$free_dimer~ wt_trial$Total_protein, span=0.3) 
  
  
  fraction_folded1<- exp(-deltaG_folding1/(R*Temp))/(1+exp(-deltaG_folding1/(R*Temp)))
  
  function_protein<- mydeltas[1]*fraction_folded1
  
  free_dimer= predict(loess_linear, function_protein) 
  
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
    list( pr = log2(2^11.761-(2^11.761-23.24125)*Propor_sup)- log2(23.24125), prm= log2(32+ (1024-32)*Propor_activation ) - 5)) 
} 

##########
 ###### 
 pr_prm_Intramol_protein_withDNA_whichparam <- function(mywtamount,whichparam1, param1,whichparam2, param2) { 
   
   # mydeltas<- c(total_protein1, deltadeltaF1_2, deltadeltadimer3, deltadeltaB1_4, deltadeltacoop5, ddg_or1_6, ddg_or2_7, ddg_or3_8) 
   mydeltas<- c(mywtamount, 0, 0, 0, 0, 0, 0, 0)
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
   
   configaration<-c() 
   free_dimer<- c() # free dimer ones 
   for (i in 2:100) { 
     free_dimer[1]<- 10^(-40) 
     free_dimer[i]<- free_dimer[i-1]*2.4
   }  
   Total_protein <- c()
   Propor_sup<- c()
   sum_config<- c()
   r_formd<- c()
   
   library(rootSolve)
   
   for (j in 1:100) { 
     for (i in 1:8) {
       configaration[i]<- exp(-relations[[1]][i]/(R*Temp))*free_dimer[j]^relations[[2]][i]
       sum_config<- sum(configaration)
       config_prob<- configaration/sum_config 
       r_formd[j]<- sum(config_prob[i]*relations[[2]][i])*2*DNA
     }
     Total_protein[j]= (free_dimer[j]/Ka)^0.5+ 2*free_dimer[j]+ r_formd[j]
   }
   wt_trial<- data.frame(Total_protein, free_dimer)
   loess_linear<- loess(wt_trial$free_dimer~ wt_trial$Total_protein, span=0.3) 
   
   
   fraction_folded1<- exp(-deltaG_folding1/(R*Temp))/(1+exp(-deltaG_folding1/(R*Temp)))
   
   function_protein<- mydeltas[1]*fraction_folded1
   
   free_dimer= predict(loess_linear, function_protein) 
   
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
     list( pr = log2(2^11.761-(2^11.761-23.24125)*Propor_sup)- log2(23.24125), prm= log2(32+ (1024-32)*Propor_activation ) - 5)) 
 }  

