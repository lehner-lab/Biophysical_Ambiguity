setwd('~/dir')
### Dose-response curve to resolve biochemical ambiguities?
### For Fig. S6
source('Forward_function_param_to_Output.R')

load('Even_space_4_param_singles.RData')
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # first gray, second near yellow, third near blue


totry<- list(
  pr= sings[[1]][sings[[1]]$phen==1.5, ], # 4 rows
  prm=sings[[2]][sings[[2]]$phen==2.7, ]  # 7 rows
)
mydeltas<- list(
  f= c(totry[[1]][totry[[1]]$whichparam=='f', 'param'], totry[[2]][totry[[2]]$whichparam=='f', 'param']), 
  d=  c(totry[[1]][totry[[1]]$whichparam=='d', 'param'], totry[[2]][totry[[2]]$whichparam=='d', 'param']), 
  b= c(totry[[1]][totry[[1]]$whichparam=='b', 'param'], totry[[2]][totry[[2]]$whichparam=='b', 'param']), 
  t=  c(totry[[1]][totry[[1]]$whichparam=='t', 'param'], totry[[2]][totry[[2]]$whichparam=='t', 'param'], 0) # to includ WT  
)

pars<- unlist(mydeltas)
phes<- c(rep(c('pr_1.5', 'prm_2.7', 'prm_2.7'),3), 'pr_1.5', 'prm_2.7', 'wt')

# with above parameters, check with changing concentraiton of the protein, how dose-response changes 
myamount<- 10^(seq(-9, -5, by=0.1)) # 41 vals 
doses<- expand.grid(param= mydeltas[[1]], amount=myamount)
doses$whichParam= 1 # folding 
for (i in 2:4){
  x<- expand.grid(param= mydeltas[[i]], amount=myamount)
  x$whichParam = i
  doses<- rbind(doses, x)
} #252
for (i in 1:length(rownames(doses))) {
  doses[i,'group'] = phes[which(pars==doses[i, 'param'])]
}


rm(pars)
rm(phes)

for (i in 1:length(rownames(doses))) {
  x<- pr_prm_Intramol_protein_withDNA_whichparam(doses[i,'amount'], doses[i,'whichParam']+1, doses[i,'param'],5, 0) 
  doses[i,'PR'] = x[[1]] 
  doses[i, 'PRM'] =x[[2]]
}

save(doses, file='dose_response_mutation_example.RData')
# plot 

ggplot(data= doses[doses$group!='prm_2.7', ]) + geom_line(aes(x= amount, y=PR, col=as.factor(whichParam), lty= group)) + 
  theme_classic() + scale_x_log10() +  scale_color_manual(values= cbPalette[c(2, 6:8)]) + 
  geom_vline(xintercept = 8.352848e-07, lty=2, col='gray') + geom_hline(yintercept = 1.5, lty=2, col='gray')

ggplot(data= doses[doses$group!='pr_1.5', ]) + geom_line(aes(x= amount, y=PRM, col=as.factor(whichParam), lty= group, shape=as.factor(param))) + 
  theme_classic() + scale_x_log10() +  scale_color_manual(values= cbPalette[c(2, 6:8)]) +
  geom_vline(xintercept = 8.352848e-07, lty=2, col='gray') + geom_hline(yintercept = 2.7, lty=2, col='gray')

