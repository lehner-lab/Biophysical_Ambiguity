##### Data generation for model ##### 
setwd('~/dir') # set to the right director where the codes are. 
source('Forward_function_param_to_Output.R')
source('Reverse_function_pheno_to_biochem.R')

library(rootSolve)

######## 1. Single mutants with evenly distributed ddGs 
myws<- list(
  pr= seq(0.05, 7.25, by=0.05), 
  prm= seq(0.05, 5.00, by= 0.05)
)
#### to store parameters temporarly 
param_vals<- list(
  pr= list(f= list(), d=list(), b=list(), t=list()), 
  prm=list(f= list(), d=list(), b=list(), t=list())
)

param_vals_unlis<- list(
  pr= list(f= c(), d=c(), b=c(), t=c()), 
  prm=list(f= c(), d=c(), b=c(), t=c())
)

#### 
for(i in 1:2){
  for(j in 1:4){
    param_vals[[i]][[j]]<- lapply(myws[[i]], function(v){Pheno_to_param_list[[j]](8.352848e-07,v, i)}) 
    names(param_vals[[i]][[j]]) = myws[[i]]
    param_vals_unlis[[i]][[j]] = unlist(param_vals[[i]][[j]])
  }
}

even_space_sing<- list(
  pr=list(), 
  prm=list()
)

for(i in 1:2) {
  for(j in 1:4){
    even_space_sing[[i]][[j]]<-data.frame(param_vals_unlis[[i]][[j]] ) 
    names(even_space_sing[[i]][[j]])= 'param' 
    even_space_sing[[i]][[j]]$phen_arb= as.numeric(rownames(even_space_sing[[i]][[j]])) 
    even_space_sing[[i]][[j]]$whichparam= names(param_vals_unlis[[i]])[j]
  }
} # names arbitarily with some errors, when two different values are asigned for each name. 
### It happens due to the ambiguity in prm predction. We will correct them by correcting it through numeric maniputlation.  
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

for (j in 1:4) {
  even_space_sing[[2]][[j]][ even_space_sing[[2]][[j]]$phen_arb>10, 'phen_arb'] = as.numeric(specify_decimal(even_space_sing[[2]][[j]][even_space_sing[[2]][[j]]$phen_arb>10, 'phen_arb']*0.1, 0))
}

for(j in 1:4) {
  for (k in 1:length(rownames(even_space_sing[[2]][[j]]))) {
    if (as.numeric(specify_decimal(even_space_sing[[2]][[j]][k, 'phen_arb']*100, 0)) %% 5 ==0) {
      even_space_sing[[2]][[j]][k, 'phen'] = as.numeric( specify_decimal(even_space_sing[[2]][[j]][k, 'phen_arb'], 2))
    } else {
      even_space_sing[[2]][[j]][k, 'phen'] = as.numeric( specify_decimal(even_space_sing[[2]][[j]][k, 'phen_arb'], 1))
    }
  }
} 

for (j in 1:4) {
  even_space_sing[[1]][[j]]$phen = even_space_sing[[1]][[j]]$phen_arb
}

sings<- list(pr= data.frame(), prm=data.frame())
for(i in 1:2){
  sings[[i]]= even_space_sing[[i]][[1]][, c('param',  'whichparam', 'phen')]
  for (j in 2:4){
    sings[[i]]= rbind(sings[[i]], even_space_sing[[i]][[j]][, c('param',  'whichparam', 'phen')])
  } 
  cat(i)
  print(table(sings[[i]]$whichparam))
  print(summary(sings[[i]]))
}
save(sings, file= 'Even_space_4_param_singles.RData')

### 2. Generate the double mutants : 
######## 2.1 Double mutants with the known biochemical types, with its own 
library(reshape2) 

doubs<- list(pr= data.frame(), prm=data.frame())
for (i in 1:2) {
  xx<- expand.grid(pheno1= sings[[i]]$phen,pheno2= sings[[i]]$phen )
  x1<- expand.grid(type1= sings[[i]]$whichparam,type2= sings[[i]]$whichparam )
  x2<- expand.grid(param1= sings[[i]]$param,param2= sings[[i]]$param )
  
  doubs[[i]] = cbind(xx, x1, x2) 
  cat(i) 
  print(head(doubs[[i]]))
}

for(i in 1:2) {
  doubs[[i]]$paramtyp= NA 
  doubs[[i]][doubs[[i]]$type1== doubs[[i]]$type2, 'paramtyp'] = 'same'
  doubs[[i]][doubs[[i]]$type1!= doubs[[i]]$type2, 'paramtyp'] = 'diff'
  
  cat(i) 
  print(table(doubs[[i]]$paramtyp))
}

rm(xx)
rm(x1)
rm(x2)

#######
######  calculate double 

t1<- unique(as.character(doubs[[2]]$type1)) # 4 
t1_key<- c(2,3,4,5) # to input to 'whichparameter' of the Forward Function  


mywts<-8.352848e-0 #  high, low as in the earlier. 
mywts_out<- list()

for(j in 1:2){
  x<- pr_prm_Intramol_protein_withDNA_whichparam(mywts,0,0,0,0)
  mywts_out[[j]]= x[[j]]
  
  print(mywts_out[[j]])
} #  pr= 0.03, prm= 2.75 


pr_prm_Intramol_protein_withDNA_whichparam(mywts, t1_key[which(t1== doubs[[1]][1,'type1'])], param1 = doubs[[1]][1,'param1'], 
                                           t1_key[which(t1== doubs[[1]][1,'type2'])], param2 = doubs[[1]][1,'param2'])[[1]] 

### function works 
for (i in 1:2) {
  doubs[[i]]$pheno = apply(doubs[[i]], 1, function (v) {
    pr_prm_Intramol_protein_withDNA_whichparam (mywts, t1_key[which(t1== v[3]) ], as.numeric( v[5]), 
                                                t1_key[which(t1== v[4]) ],  as.numeric(v[6])) [[i]] 
  }
  )
}
for (i in 1:2) {
  doubs[[i]]$whichparam= paste(doubs[[i]]$type1, doubs[[i]]$type2, sep=',')
}

save(doubs, file= 'Even_space_4_param_doubs.RData')
#######
###### 3. make the ambiguity data frame

mycombos_ci<- list(pr=list(ones= list(),twos=list(), threes=list(), fours=list()), 
                   prm=list(ones= list(), twos=list(), threes=list(), fours=list()))
### make the ambiguity dataframe 

mytypCombs<- list()
for(i in c(1:4)){
  mytypCombs[[i]]= combn(t1, i, simplify = F)
}

mycombo<-list()

for (i in 1:2){  
  mycombo[[i]]<- expand.grid(p1= as.numeric(unique(doubs[[i]]$pheno1)), p2=as.numeric(unique(doubs[[i]]$pheno1))) }

for (i in 1:2){  
  for (j in 1:4){  
    for(k in 1:length(mytypCombs[[j]]))
      mycombos_ci[[i]][[j]][[k]]= mycombo[[i]]
  }
}

for (i in 1:2){ 
  
  for(l in 1:length(rownames(mycombo[[i]]))) {
    x1 = doubs[[i]][ doubs[[i]]$pheno1== mycombo[[i]][l,'p1'] & doubs[[i]]$pheno2== mycombo[[i]][l,'p2'], c('type1','type2','pheno') ]  
    
    
    for (j in 1:4){
      for(k in 1:length(mytypCombs[[j]])) {
        xx= x1[x1$type1 %in% mytypCombs[[j]][[k]] & x1$type2 %in% mytypCombs[[j]][[k]], ]
        
        if (length(rownames(xx))>=1) {
          mycombos_ci[[i]][[j]][[k]][l, 'ambiNum'] = length(unique(xx$pheno)) 
          mycombos_ci[[i]][[j]][[k]][l, 'ambiDiff']= max(xx$pheno,na.rm=T) - min(xx$pheno,na.rm=T) 
          
        }  else {
          mycombos_ci[[i]][[j]][[k]][l, 'ambiNum'] = NA
          mycombos_ci[[i]][[j]][[k]][l, 'ambiDiff']= NA
        }
      }
    }
  }
}


for (i in 1:2){ 
  for (j in 1:4){  
    for(k in 1:length(mytypCombs[[j]])) {   
      mycombos_ci[[i]][[j]][[k]]<- mycombos_ci[[i]][[j]][[k]][ !is.na(mycombos_ci[[i]][[j]][[k]]$ambiNum), ]
    }
  }
}


library(stringr)
for (i in 1:2){ 
  for (j in 1:4){  
    for(k in 1:length(mytypCombs[[j]])) { 
      mycombos_ci[[i]][[j]][[k]]$allowed_biochem = j 
      
      mycombos_ci[[i]][[j]][[k]]$which_biochem = paste(paste(mytypCombs[[j]][[k]][1],mytypCombs[[j]][[k]][2],sep=''), paste(mytypCombs[[j]][[k]][3], mytypCombs[[j]][[k]][4], sep=""), sep='')
    }
  }
}

mycombo_ci_inrows<- list(pr=data.frame(), prm=data.frame())

for(i in 1:2) { 
  mycombo_ci_inrows[[i]] = mycombos_ci[[i]][[1]][[1]]
  for (j in 2:4) {
    mycombo_ci_inrows[[i]]= rbind(mycombo_ci_inrows[[i]],mycombos_ci[[i]][[j]][[1]])
  }
  
  for (j in 1:3) {
    for (k in 2:length(mytypCombs[[j]])) {
      mycombo_ci_inrows[[i]] = rbind(mycombo_ci_inrows[[i]], mycombos_ci[[i]][[j]][[k]] )
    }
  }
}

for(i in 1:2) {
  mycombo_ci_inrows[[i]]$involvF= 0 
  mycombo_ci_inrows[[i]][str_sub(mycombo_ci_inrows[[i]]$which_biochem, 1,1)=='f', 'involvF'] =1
}

mycombo_ci_inrows[[1]]$whichpheno= 'PR'
mycombo_ci_inrows[[2]]$whichpheno= 'PRM'


bytypes<- list(pr=data.frame(), prm=data.frame())
mytps<- unique(mycombo_ci_inrows[[1]]$which_biochem) 
bioche_no<- c(1,2,3,4,1,1,1,2,2,2,2,2,3,3,3)
invol_f<- c(1,1,1,1,0,0,0,1,1,0,0,0,1,1,0)

for (i in 1:2) {
  bytypes[[i]]= data.frame(mytps, bioche_no, invol_f)
  for (l in 1:length(mytps)) {
    xx= mycombo_ci_inrows[[i]][mycombo_ci_inrows[[i]]$which_biochem==mytps[l],  ] 
    
    bytypes[[i]][l, 'Max_ambiNum'] = max(xx$ambiNum)
    bytypes[[i]][l, 'Med_ambiNum'] = median(xx$ambiNum)
    bytypes[[i]][l, 'Min_ambiNum'] = min (xx$ambiNum)
    
    bytypes[[i]][l, 'Max_ambiDiff'] = max(xx$ambiDiff)
    bytypes[[i]][l, 'Med_ambiDiff'] = median(xx$ambiDiff)
    bytypes[[i]][l, 'Min_ambiDiff'] = min (xx$ambiDiff)
    
    bytypes[[i]][l, 'Prop_ambig'] = length(rownames(xx[xx$ambiNum>1,])) / length(rownames(xx))
    
  }
  
}

biochem_ambi<- list(keys= mytypCombs, val_lists= mycombos_ci, inrows_with1_4= mycombo_ci_inrows, bytypes_stat= bytypes)
save(biochem_ambi, file= 'biochem_ambiguity_4_param_doub_lists.RData') 

## ambiguity of PRM with known two different paraemter-combinations, for Fig.S3 
prm_2_fix<- doubs[[2]][doubs[[2]]$paramtyp=='diff',  ]
mycombo<- expand.grid(p1= as.numeric(unique(prm_2_fix$pheno1)), p2=as.numeric(unique(prm_2_fix$pheno1))) 
mycombos<- list()
library(stringr)
for(j in 1:length(unique(as.character(prm_2_fix$types)))) {
  mycombos[[j]]= mycombo
  mycombos[[j]]$types= unique(as.character(prm_2_fix$types))[j] 
  mycombos[[j]]$type1= str_sub(mycombos[[j]]$types, 1,1)
  mycombos[[j]]$type2= str_sub(mycombos[[j]]$types, -1,-1)
}

for(l in 1:length(rownames(mycombo))) {
  x1 = prm_2_fix[prm_2_fix$pheno1== mycombo[l,'p1'] & prm_2_fix$pheno2== mycombo[l,'p2'],]  
  for (j in 1:length(unique(as.character(prm_2_fix$types)))) {
    
    xx= x1[x1$types == unique(as.character(prm_2_fix$types))[j], ] 
    
    if (length(rownames(xx))>=1) {
      mycombos[[j]][l, 'ambiNum'] = length(unique(xx$pheno)) 
      mycombos[[j]][l, 'ambiDiff']= max(xx$pheno,na.rm=T) - min(xx$pheno,na.rm=T) 
      
    }  else {
      mycombos[[j]][l, 'ambiNum'] = NA
      mycombos[[j]][l, 'ambiDiff']= NA
    }
  }
  
}

for (i in 1: length(mycombos)) {
  mycombos[[i]]= mycombos[[i]][!is.na(mycombos[[i]]$ambiNum), ]
} 
mycombs_inrow= mycombos[[1]]
for (i in 2:12){
  mycombs_inrow= rbind(mycombs_inrow, mycombos[[i]])
}
mycombs_inrow$type1= factor(mycombs_inrow$type1, levels=c('f', 'd', 'b', 't'))
mycombs_inrow$type2= factor(mycombs_inrow$type2, levels=c('f', 'd', 'b', 't'))
mycombs_inrow$types=  factor(mycombs_inrow$types, levels=levels(prm_2_fix$types))

save(mycombs_inrow, file='biochem_ambiguity_doubs_twodiff_params_lists.RData')

### 4. Pleiotropic mutants
######## 4.1 Phenotypes to Parameter combinations 
source('Reverse_function_for_pleiotropic_muts.R')

#4.1.1 when mtuations are pleiotropic for two paraemters - folding and binding, how do they interact with others? 
# Fig. 4 
myws<- list(c(0.1, 1,3,5,7), c(0.1, 1, 2.7, 4)) # pr and prm: to try a range of the phenotypes 
## But in the figure only with the phenotype =1 is plotted
mydgs<- seq(-1, 4.5, by=0.1) # This is the folding energy

pleio2<- list(
  pr=list(p0.1=c(),p1=c(), p3=c(), p5=c(),p7=c() ), 
  prm=list(p0.1=c(),p1=c(), p2.7=c(), p4=c()))

for (i in 1:2){
  for (j in 1:length(myws[[i]])) {
    parms<- list()
    for (k in 1:length(mydgs)) {
      parms[[k]]<- Pheno_to_param_pleiotropy_FDB(8.352848e-07,mydgs[k], 0,myws[[i]][j], i) 
    }
    names(parms)<- mydgs
    pleio2[[i]][[j]] = unlist(parms)
    
  }
  
}
pleio<- list()
for (i in 1:2){
  pleio[[i]] = data.frame() 
  xx<- list()
  for (k in 1:length(names(pleio2[[i]]))) {
    xx[[k]]<- data.frame(ddG2= pleio2[[i]][[k]] , ddG1= as.numeric(names(pleio2[[i]][[k]] ) )) 
    xx[[k]]$pheno = as.numeric(rep(myws[[i]][k], length(rownames(xx[[k]]))))
  }
  pleio[[i]]= xx[[1]] 
  for (k in 2: length(names(pleio2[[i]]))) {
    pleio[[i]]= rbind(pleio[[i]], xx[[k]])
  }
}
pleio[[1]]$whichPhen= 'PR'
pleio[[2]]$whichPhen= 'PRM'

m(xx)
##### names are given as the phenotypes, but when the same phenotype is assigned to two or more parameters, they were wrong. 
# Correct them. 
pleio[[2]][pleio[[2]]$ddG1< -1, 'ddG1'] = pleio[[2]][pleio[[2]]$ddG1< -1, 'ddG1']*0.1
pleio[[2]][pleio[[2]]$ddG1< -1, 'ddG1'] = as.numeric(specify_decimal(pleio[[2]][pleio[[2]]$ddG1< -1, 'ddG1'], 0))

pleio[[2]][pleio[[2]]$ddG1> 10, 'ddG1'] = as.numeric(specify_decimal(pleio[[2]][pleio[[2]]$ddG1> 10, 'ddG1']*0.1, 0))

pleio[[2]]$ddG1 = as.numeric(specify_decimal(pleio[[2]]$ddG1, 1))
save(pleio, file='pleiotropic_FB_list.RData')

library(reshape2)
whichParam= c('f', 'b', 'd', 't')
doubs_with_pleio<- list(pr=list(), prm=list())
for (i in 1:2) {
  for(j in 1:4){
    p2<- sings[[i]][sings[[i]]$whichparam==whichParam[j], ]
    pp<- expand.grid(P1=pleio[[i]]$pheno, P2= p2$phen)
    params1<- expand.grid(P1_ddG_f=pleio[[i]]$ddG1, P2_ddG= p2$param)
    params2<- expand.grid(P1_ddG_b=pleio[[i]]$ddG2, P2_ddG= p2$param)
    doubs_with_pleio[[i]][[j]]= cbind(pp, params1, params2) 
    doubs_with_pleio[[i]][[j]]$P2_which= whichParam[j] 
    doubs_with_pleio[[i]][[j]]$promoter = names(doubs_with_pleio)[i]
  }
}
for (i in 1:2) {
  for(j in 1:4){
    
    doubs_with_pleio[[i]][[j]]= doubs_with_pleio[[i]][[j]][, c(1:3,5:8) ]
  }
}

for (i in 1:2) {
  doubs_with_pleio[[i]][[1]]$ddG_f= doubs_with_pleio[[i]][[1]]$P1_ddG_f+ doubs_with_pleio[[i]][[1]]$P2_ddG
  doubs_with_pleio[[i]][[2]]$ddG_b= doubs_with_pleio[[i]][[2]]$P1_ddG_b+ doubs_with_pleio[[i]][[2]]$P2_ddG } 

for (i in 1:2) {
  doubs_with_pleio[[i]][[1]]$doubP= apply(  doubs_with_pleio[[i]][[1]], 1, function(v) {pr_prm_Intramol_protein_withDNA_four_param(8.352848e-07,as.numeric(v[8]), 0,  as.numeric(v[4]),0)[[i]]})
  doubs_with_pleio[[i]][[2]]$doubP= apply(  doubs_with_pleio[[i]][[2]], 1, function(v) {pr_prm_Intramol_protein_withDNA_four_param(8.352848e-07,as.numeric(v[3]),0,as.numeric(v[8]),0)[[i]]})
  
  doubs_with_pleio[[i]][[3]]$doubP= apply(  doubs_with_pleio[[i]][[3]], 1, function(v) {pr_prm_Intramol_protein_withDNA_four_param(8.352848e-07,as.numeric(v[3]), as.numeric(v[5]), as.numeric(v[4]), 0)[[i]]})
  doubs_with_pleio[[i]][[4]]$doubP= apply(  doubs_with_pleio[[i]][[4]], 1, function(v) {pr_prm_Intramol_protein_withDNA_four_param(8.352848e-07,as.numeric(v[3]), 0, as.numeric(v[4]), as.numeric(v[5]))[[i]]})
  
}
save(doubs_with_pleio, file='pleiotropic_doub_with_one_not_list.RData')

###### two pleiotropic mutants' combinations 
##### pleio + pleio pheno combo. 
doubs_with_pleio_pleio<- list()
for (i in 1:2) {
  pp<- expand.grid(P1=pleio[[i]]$pheno, P2= pleio[[i]]$pheno)
  params1<- expand.grid(P1_ddG_f=pleio[[i]]$ddG1, P2_ddG_f= pleio[[i]]$ddG1)
  params2<- expand.grid(P1_ddG_b=pleio[[i]]$ddG2, P2_ddG_b= pleio[[i]]$ddG2)
  doubs_with_pleio_pleio[[i]]= cbind(pp, params1, params2) 
  doubs_with_pleio_pleio[[i]]$promoter = names(doubs_with_pleio)[i] 
  
  doubs_with_pleio_pleio[[i]]$ddG_f= doubs_with_pleio_pleio[[i]]$P1_ddG_f +doubs_with_pleio_pleio[[i]]$P2_ddG_f 
  doubs_with_pleio_pleio[[i]]$ddG_b= doubs_with_pleio_pleio[[i]]$P1_ddG_b +doubs_with_pleio_pleio[[i]]$P2_ddG_b 
}

for (i in 1:2) {
  doubs_with_pleio_pleio[[i]]$pheno = apply(doubs_with_pleio_pleio[[i]], 1, function(v) {pr_prm_Intramol_protein_withDNA_four_param(8.352848e-07,as.numeric(v[8]),0,as.numeric(v[9]),0)[[i]]})
} 

save(doubs_with_pleio_pleio, file='pleiotropic_doub_both.RData')


######## 4.2 Pleiotropic mutants affecting two different parameters, to examine the phenotypic landscape
dds<- seq(-1,5, by=0.1) 
dat<- expand.grid(ddx= dds, ddy= dds)
myPleioSing<- list()

whichtwos<- list(
  FD= c(1,2), 
  FB=c(1,3), 
  FT= c(1,4), 
  DB= c(2,3), 
  DT= c(2,4), 
  BT= c(3,4)
)
for(i in 1:6){
  myPleioSing[[i]]= dat 
  for (j in 1:length(rownames(myPleioSing[[i]]))) {
    x1<- pr_prm_Intramol_protein_withDNA_whichparam(8.352848e-07,whichtwos[[i]][1]+1, dat[j,1] ,whichtwos[[i]][2]+1,dat[j,2])
    myPleioSing[[i]][j,'PR'] = x1[[1]]
    myPleioSing[[i]][j,'PRM'] = x1[[2]]
  }
}
for(i in 1:6){
  #cat(i)
  #print(summary(myPleioSing[[i]]))
  myPleioSing[[i]]$id= paste(as.character(myPleioSing[[i]]$ddx), as.character(myPleioSing[[i]]$ddy), sep=',')
}
### use singles to make the landscape 
### use the doubles to plot the area
names(myPleioSing)<- names(whichtwos)

source('Function_decimal.R')
for(i in 1:6){
  myPleioSing[[i]]$PR_round= as.numeric(specify_decimal(myPleioSing[[i]]$PR, 2))
  myPleioSing[[i]]$PRM_round= as.numeric(specify_decimal(myPleioSing[[i]]$PRM, 2)) 
  myPleioSing[[i]]$whichTwo = names(whichtwos)[i]
}
save(myPleioSing, file='Pleio_sings_forward_all_for_later.RData' )

######## 4.3 Pleiotropic mutants affecting three different parameters for the summary plot 
## combinations of ddG folding and ddG dimerizations to search for ddG binding energy that produces the given phenotypes: 
# PR= 0.1 (for mutation A ) and PR=3 (for mutation B)
mydgs<- seq(-1, 5, by=0.05) 
mydDs<- seq(-1, 5, by=0.05)  

pleio3<- list(
  p0.1=data.frame(),
  p3= data.frame()
)

for (i in 1:2){
  pleio3[[i]]= data.frame(expand.grid(mydgs, mydDs)) 
}
for (i in 1:2){
  names(pleio3[[i]]) = c('ddF', 'ddD')  
  pleio3[[i]]=pleio3[[i]][pleio3[[i]]$ddF+pleio3[[i]]$ddD<5.8,  ] 
  # the condition is set to remove those that can not reach the given phentoype with whichever ddB based on the forward- modeling 
  pleio3[[i]]$phen= myws[[i]]
} # 10986 values 

mypars<- list(p0.1=c(), p3=c())

for (j in 1:length(myws)) {
  parms<- list()
  for (k in 1:length(rownames(pleio3[[j]]))) {
    parms[[k]]<- Pheno_to_param_pleiotropy_FDB(8.352848e-07,as.numeric(pleio3[[j]][k,1]), 
                                               as.numeric(pleio3[[j]][k,2]),myws[j], 1) 
  }
  names(parms)<- seq(1:10986)
  mypars[[j]] = unlist(parms)
  
}

for (i in 1:2){
  pleio3[[i]]$ddB=NA
  for(j in 1:length(mypars[[i]])){ 
    pleio3[[i]][as.numeric(as.character(names(mypars[[i]][j]))), 'ddB'] = mypars[[i]][j]
  }
  pleio3[[i]]=pleio3[[i]][!is.na(pleio3[[i]]$ddB), ]
  
}

save(pleio3, file='Reverse_anyof3_param_0.1_3_as_prPheno.RData')

ABs<- cbind(
  expand.grid(ddF1=pleio3[[1]]$ddF, ddF2=pleio3[[2]]$ddF), 
  expand.grid(ddD1=pleio3[[1]]$ddD, ddD2=pleio3[[2]]$ddD), 
  expand.grid(ddB1=pleio3[[1]]$ddB, ddB2=pleio3[[2]]$ddB) 
) # 75,219,260
ABs$ddF= ABs$ddF1+ ABs$ddF2
ABs$ddD= ABs$ddD1+ ABs$ddD2
ABs$ddB= ABs$ddB1+ ABs$ddB2

for (i in 1:length(rownames(ABs))) {
  x1<- pr_prm_Intramol_protein_withDNA_four_param(8.352848e-07, ABs[i,'ddF1'] + ABs[i,'ddF2'], ABs[i,'ddD1'] + ABs[i,'ddD2'],ABs[i,'ddB1'] + ABs[i,'ddB2'], 0)
  ABs[i, 'PR'] = x1[[1]] 
  ABs[i, 'PRM'] =x1[[2]]
} 

save(ABs, file='Reverse_anyof3_param_0.1_3_as_prPheno_doub_comb_pr_prm75219260.RData')



