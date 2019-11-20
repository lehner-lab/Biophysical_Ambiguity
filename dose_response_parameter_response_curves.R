##### Data generation for model ##### 
setwd('~/dir') # set to the right director where the codes are. 
source('Forward_function_param_to_Output.R')
source('Reverse_function_pheno_to_biochem.R')

### 1. Generating data for dose-response curves for PR and PRM 

myamount<- 10^(seq(-9, -5, by=0.1))    

doses<- data.frame(amount= myamount, PR=rep(NA, length(myamount)), PRM= rep(NA, length(myamount)))
for (i in 1:length(myamount) ) { 
  
  x<- pr_prm_Intramol_protein_withDNA_four_param(myamount[i],0, 0,0, 0) 
  doses[i,'PR'] = x[[1]] 
  doses[i,'PRM'] =x[[2]] 
  
}

### 2. Parameter-response curves for PR and PRM 
# wild type set as 8.352848e-07
Gs_forward<- seq(-2, 7, by= 0.1) # 19 vals 
phes<- list(pr= list(), prm= list())

for (i in 1:4) {   
  for (j in 1:2) {
    phes[[j]][[i]]= sapply(Gs_forward,function(v){pr_prm_Intramol_protein_withDNA_whichparam(8.352848e-07, i+1, v,5, 0)}[[j]])
  }
}

forwards<- list(pr= list(), prm=list())
param_vals= c('1.f', '2.d','3.b','4.t') 
for (j in 1:2){
  for (i in 1:4){
    forwards[[j]][[i]] = data.frame(param= Gs_forward, pheno= phes[[j]][[i]]) 
    forwards[[j]][[i]]$whichW= names(forwards)[j] 
    forwards[[j]][[i]]$whichParam= param_vals[i] 
    
  }
}

for (j in 1:2) {
  names(forwards[[j]]) = c('1.f', '2.d','3.b','4.t') 
}

save(forwards, file= 'Forwards_4_param_list.RData') 
# saved to compare with the Reverse-function generated data later. 

####
#### to plot 

library(ggplot2)
library(ggpubr)
##### 1. Dose-response curve for wild type 

pr= ggplot(data= doses) + geom_line(aes(x= amount, y=PR)) +
  theme_classic() + scale_x_log10() + 
  geom_vline(xintercept = 8.352848e-07, lty=2, col='gray')  + ylim(0, 7.5)
prm= ggplot(data= doses) + geom_line(aes(x= amount, y=PRM)) +
  theme_classic() + scale_x_log10() + 
  geom_vline(xintercept = 8.352848e-07, lty=2, col='gray') + ylim(0, 7.5) # 8.352848e-07 is the value selected for the wild type expression where mtuations are modeled 
ggarrange(pr, prm, nrow=1, ncol=2 ) 

#### 2. Parameter-response curve
mymaxs<- c(7.5, 5)
mygs<- list(pr= list(), prm=list())
for(i in 1:2) {
  for (j in 1:4) {
    mygs[[i]][[j]] = ggplot(data=forwards[[i]][[j]]) + 
      labs(x='∆∆G', y=paste('Phenotype', names(mygs)[i]), title= param_vals[j]) + 
      geom_vline(xintercept = 0, lty=2) + theme_classic() + xlim(-2, 7) + ylim(-0.5, mymax[i]) + 
      geom_line(aes(x=param, y= pheno), size=1) 
    
  }
}

ggarrange(mygs[[1]][[1]], mygs[[1]][[2]],mygs[[1]][[3]],mygs[[1]][[4]],
          mygs[[2]][[1]], mygs[[2]][[2]],mygs[[2]][[3]],mygs[[2]][[4]], 
          nrow=2, ncol=4)
