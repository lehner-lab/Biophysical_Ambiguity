##### Pleiotropic for three parameters - Folding, Dimerization, and Binding to DNA 
#### Inputting folding & dimerization to look for binding energies that produce the given phenotype 
Pheno_to_param_pleiotropy_FDB<- function(total_protein, ddG1,ddG2, W, whichpheno){ # for folding
  # whichpheno= 1 pr ; whichpheno =2, prm 
  my_protein = total_protein
  # protein amount 8.352848e-07
  inverse = function (f, lower = -2, upper = 10) { # minimum interval and muaximum interval
    function (y) uniroot.all ((function (x) f(x) - y), lower = lower, upper = upper)
  }
  
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  Ka= 5 *10^7 # ka for dimmer formation 
  DNA= 10^(-9)
  
  library(rootSolve)
  
  myphe= whichpheno
  pheno= W
  forward<- function( ddG3) { #  1 folding, 2 dimerization, 3 binding 
    
    wt_deltaG=  -1 # set as in the original analysis 
    deltaG_folding1= wt_deltaG+ ddG1
    
    #configarations=list(c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1), c(1,1,0), c(1,0,1), c(0,1,1), c(1,1,1)) # 8 configurations OR1, OR2, OR3. 
    wt_BdeltaGs<- c(-11.7, -10.1,-10.1 ) 
    wt_BdeltaGs<- wt_BdeltaGs + ddG3
    # both on proteins
    
    coop= -2+ ddG2  # OR1 and OR2 ; or OR2 and OR3
    
    relations=list(deltaGs=c(0, wt_BdeltaGs[1], wt_BdeltaGs[2], wt_BdeltaGs[3],  wt_BdeltaGs[1]+ wt_BdeltaGs[2]+ coop,  wt_BdeltaGs[1]+ wt_BdeltaGs[3],  wt_BdeltaGs[3]+ wt_BdeltaGs[2]+ coop, sum(wt_BdeltaGs)+ coop),
                   dimer_number=c(0,1,1,1,2,2,2,3)) # thermo_parameters
    
    ##### empirical search
    
    ##### empirical search
    
    Mytotal_protein<- function(free_dimer) {
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
    
    free_dimer=   unlist(optimize(Mytotal_protein, c(10^-40, 10^-6), tol= 10^-23)[1]) # from 10^-21 changed to -40
    
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
    
    if (myphe== 1) {
      return(pr = 1- Propor_sup)
    } else {
      return( prm= Propor_activation) 
    }
  }
  
  forward=Vectorize(forward) # added in case uniroot.all fails because the function is not vectorized
  ddG_cal= inverse(forward, -2, 10 )
  ddG3= unlist(ddG_cal(pheno))
  
  return(ddG3)
}