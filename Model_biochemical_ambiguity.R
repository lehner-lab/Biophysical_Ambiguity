rm(list=ls())
setwd('~/dir/')
##### 
######### functions 

source('Forward_function_param_to_Output.R')
source('Reverse_function_pheno_to_biochem.R')
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
source('dose_response.R')
source('Function_reverse_for_pleiotropic_20200310.R')



####### For Figure 1 and S1. The relationships between parameters and the outputs

mywts<-8.352848e-07 #  high, low as in the earlier. 

########
##### For Fig1. dose-response for the wt 
Total<- 10^(seq(-10, -4, by=0.1)) # free ones 
pr<-c() 
prm<-c()
for (i in 1:length(Total)) {
  x<- Forward(Total[i], 0, 0, 0, 0) 
  pr[i]=x[[1]] 
  prm[i] = x[[2]]
  }

dose_response<- data.frame(Total_protein= Total, PR= pr, prm=prm)
save(dose_response, file='Dose_response_WT.RData')


rm(dose_response) 

#########
#### ∆∆ G vs phenotypes 
myd<-seq(-2, 10, by=0.2) # 61 vals 
myds<- data.frame(ddg= rep(myd, 4), whichPar= c(rep(2, length(myd)), rep(3, length(myd)),rep(4, length(myd)),rep(5, length(myd))))

for (i in 1:length(rownames(myds))) {
  x<- pr_prm_Intramol_protein_withDNA_whichparam(mywts, myds[i,'whichPar'] ,myds[i, 'ddg'], 2, 0)
  myds[i,'pr'] = x[[1]] 
  myds[i,'prm'] = x[[2]]
}
save(myds, file='ddG_to_output.RData')

######  For Fig. S1, how do different parameters affect dose-response curves?  
########

free_dimer<- c() # free ones 
for (i in 2:50) { 
  free_dimer[1]<- 10^(-20) 
  free_dimer[i]<- free_dimer[i-1]*2
}  
myddgs<- c(-2, 0, 2, 4) 
fre2<- expand.grid(free_dimer= free_dimer, ddG1= myddgs, ddG2= myddgs, ddG3= myddgs, ddG4= myddgs)
rm(i)
for (i in 1:length(rownames(fre2))) { 
  x=from_free_dimer(fre2[i, 1], fre2[i, 2],fre2[i, 3],fre2[i, 4],fre2[i, 5])  
  fre2[i, 'Total_protein'] = x[[1]] 
  fre2[i, 'pr'] = x[[2]] 
  fre2[i, 'prm'] = x[[3]] 
}

save(fre2, file= 'starting_from_free_dimer.RData')
rm(fre2) 


#############
########### generating log-evenly distributed phenotypes to study. 
myw_outs<- Forward(8.352848e-07, 0,0,0,0)
mysampl<- seq(-13.5, 0, by= 0.1)# log even,136 values

myws<- list(
  pr= 2^mysampl,  # seq(0.0001, 1, by=0.01), 
  prm= 2^mysampl # seq(0.0001, 1, by= 0.01)
) # 100 values



param_vals<- list(
  pr= list(f= list(), d=list(), b=list(), t=list()), 
  prm=list(f= list(), d=list(), b=list(), t=list())
)

param_vals_unlis<- list(
  pr= list(f= c(), d=c(), b=c(), t=c()), 
  prm=list(f= c(), d=c(), b=c(), t=c())
)

even_space_sing<- list(
  pr=list(), 
  prm=list()
)

#### 
for(i in 1:2){
  for(j in 1:4){
    param_vals[[i]][[j]]<- lapply(myws[[i]], function(v){Pheno_to_param_list[[j]](8.352848e-07,v, i)}) 
     names(param_vals[[i]][[j]]) = log2(myws[[i]])# to check later
    myPhes<- c() 
    myPhes= c(rep(as.numeric(names(param_vals[[i]][[j]])[1]), length(param_vals[[i]][[j]][[1]])))
    
  for (k in 2:length(myws[[i]])) { 
    myPhes= c(myPhes, rep(as.numeric(names(param_vals[[i]][[j]])[k]), length(param_vals[[i]][[j]][[k]])))
  }
    
    even_space_sing[[i]][[j]]= data.frame(param= unlist(param_vals[[i]][[j]]), phen= myPhes, 
                                          whichparam = names(param_vals_unlis[[i]])[j])
  }
}


sings<- list(pr= data.frame(), prm=data.frame())
for(i in 1:2){
  sings[[i]]= data.frame(even_space_sing[[i]][[1]])
  for (j in 2:4){
    sings[[i]]= rbind(sings[[i]], data.frame(even_space_sing[[i]][[j]]))
  } 
  cat(i) 
  print(table(sings[[i]]$whichparam))
  print(summary(sings[[i]]))
}


length(unique(sings[[1]]$phen))
length(unique(sings[[2]]$phen)) # due to the fact that the same phenotype can be mapped to two params

for (i in 1:2){
  names(sings[[i]])[3] = 'type1'
} # for later

save(sings, file= 'Even_space_4_param_singles_0_1_logeven_20200310.RData') # pheno in log scale
rm(param_vals)
rm(param_vals_unlis)
rm(i); rm(j); rm(mysampl)

###########
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

rm(xx); rm(x1) ; rm(x2)

#######
#######
######  calculate double 

t1<- unique(as.character(doubs[[2]]$type1)) # 4 
t1_key<- c(2,3,4,5) # to input to 'whichparameter' of the Forward Function  



### function works 
for (i in 1:2) {
  doubs[[i]]$pheno = apply(doubs[[i]], 1, function (v) {
    pr_prm_Intramol_protein_withDNA_whichparam (mywts, t1_key[which(t1== v[3]) ], as.numeric( v[5]), 
                                                t1_key[which(t1== v[4]) ],  as.numeric(v[6])) [[i]] 
  }
  )
}
for (i in 1:2) { 
  doubs[[i]]$phen_log2= log2(doubs[[i]]$pheno) 
  doubs[[i]]$whichparam= paste(doubs[[i]]$type1, doubs[[i]]$type2, sep=',')
  doubs[[i]]$whichparam= factor(doubs[[i]]$whichparam, 
                                levels= c('f,f', 'd,d', 'b,b', 't,t', 
                                          'f,d', 'f,b', 'f,t', 
                                          'd,b', 'd,t', 'b,t', 
                                           'd,f', 'b,f', 't,f', 
                                          'b,d', 't,d', 't,b'))
  
}


save(doubs, file= 'Even_space_4_param_doubs0_1_logeven_20200310.RData')



########
###### 3. make the ambiguity data frame


### decimal place 2, in order not to call the biologically irrelevantly small phenotypic differences as 
# different phenotypes 

for (i in 1:2){
  doubs[[i]]$phen_log2_decimal= as.numeric(specify_decimal( doubs[[i]]$phen_log2, 2))
}


### make the ambiguity dataframe 

mytypCombs<- list()
for(i in c(1:4)){
  mytypCombs[[i]]= combn(t1, i, simplify = F)
}

mycombo<-list()

for (i in 1:2){  
  mycombo[[i]]<- data.frame(expand.grid(p1= as.numeric(unique(doubs[[i]]$pheno1)), p2=as.numeric(unique(doubs[[i]]$pheno1)))) }

for (i in 1:2){  
  for (j in 1:4){  
    for(k in 1:length(mytypCombs[[j]]))
      mycombos_ci[[i]][[j]][[k]]= mycombo[[i]]
  }
}

for (i in 1:2){ 
  for(l in 1:length(rownames(mycombo[[i]]))) {
    x1 = doubs[[i]][ doubs[[i]]$pheno1== mycombo[[i]][l,'p1'] & doubs[[i]]$pheno2== mycombo[[i]][l,'p2'], c('type1','type2','phen_log2_decimal') ]  
    
    
    for (j in 1:4){
      for(k in 1:length(mytypCombs[[j]])) {
        xx= x1[x1$type1 %in% mytypCombs[[j]][[k]] & x1$type2 %in% mytypCombs[[j]][[k]], ]
        
        if (length(rownames(xx))>=1) {
          mycombos_ci[[i]][[j]][[k]][l, 'ambiNum'] = length(unique(xx$phen_log2_decimal)) 
          mycombos_ci[[i]][[j]][[k]][l, 'ambiDiff']= max(xx$phen_log2_decimal,na.rm=T) - min(xx$phen_log2_decimal,na.rm=T) 
          
        }  else {
          mycombos_ci[[i]][[j]][[k]][l, 'ambiNum'] = NA
          mycombos_ci[[i]][[j]][[k]][l, 'ambiDiff']= NA
        }
      }
    }
  }
}
for ( i in 1:2){
  for(j in 1:4){
    for (k in 1:length(mytypCombs[[j]])) {
      cat(i); cat(j); cat(k)
      print(summary(mycombos_ci[[i]][[j]][[k]]))
    }
  }
}


for (i in 1:2){ 
  for (j in 1:4){  
    for(k in 1:length(mytypCombs[[j]])) {   
      mycombos_ci[[i]][[j]][[k]]= mycombos_ci[[i]][[j]][[k]][!(is.na(mycombos_ci[[i]][[j]][[k]]$ambiNum)), ]
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
####### 
## ambiguity of PRM with known two different paraemter-combinations, for Fig.S3 
prm_2_fix<- doubs[[2]][doubs[[2]]$paramtyp=='diff',  ]
mycombo<- expand.grid(p1= as.numeric(unique(prm_2_fix$pheno1)), p2=as.numeric(unique(prm_2_fix$pheno1))) 

mycombos<- list()

for(j in 1:length(unique(as.character(prm_2_fix$whichparam)))) {
  mycombos[[j]]= mycombo
  mycombos[[j]]$whichparam= unique(as.character(prm_2_fix$whichparam))[j] 
  mycombos[[j]]$type1= str_sub(mycombos[[j]]$whichparam, 1,1)
  mycombos[[j]]$type2= str_sub(mycombos[[j]]$whichparam, -1,-1)
}

for(l in 1:length(rownames(mycombo))) {
  x1 = prm_2_fix[prm_2_fix$pheno1== mycombo[l,'p1'] & prm_2_fix$pheno2== mycombo[l,'p2'],]  
  for (j in 1:length(unique(as.character(prm_2_fix$whichparam)))) {
    
    xx= x1[x1$whichparam == unique(as.character(prm_2_fix$whichparam))[j], ] 
    
    if (length(rownames(xx))>=1) {
      mycombos[[j]][l, 'ambiNum'] = length(unique(xx$phen_log2_decimal)) 
      mycombos[[j]][l, 'ambiDiff']= max(xx$phen_log2_decimal,na.rm=T) - min(xx$phen_log2_decimal,na.rm=T) 
      
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
mycombs_inrow$whichparam=  factor(mycombs_inrow$whichparam, levels=levels(prm_2_fix$whichparam))

save(mycombs_inrow, file='biochem_ambiguity_doubs_twodiff_params_lists_prm2_with_log2_0.01_precision.RData')



##########
############
bytypes<- list(pr=data.frame(), prm=data.frame())
mytps<- unique(mycombo_ci_inrows[[1]]$which_biochem) 
bioche_no<- c(1,2,3,4,1,1,1,2,2,2,2,2,3,3,3)
invol_f<- c(1,1,1,1,0,0,0,1,1,0,0,0,1,1,0)

for (i in 1:2) {
  bytypes[[i]]= data.frame(mytps, bioche_no, invol_f)
  for (l in 1:length(mytps)) {
    xx= mycombo_ci_inrows[[i]][mycombo_ci_inrows[[i]]$which_biochem==mytps[l] ,  ] 
    
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
for (i in 1:2){
  biochem_ambi[[3]][[i]]$which_biochem= 
    factor(biochem_ambi[[3]][[i]]$which_biochem, 
           levels= c('fNANANA', 'dNANANA','bNANANA','tNANANA', 
                     'fdNANA', 'fbNANA', 'ftNANA', 
                     'dbNANA', 'dtNANA', 'btNANA', 
                     'fdbNA', 'fdtNA', 'fbtNA', 'dbtNA', 'fdbt'))
}


save(biochem_ambi, file= 'biochem_ambiguity_4_param_doub_01_20200324_with_precision0.01_lists.RData') 



######## 
######### with one pleiotropic mutant, how do mutations interact? 
# pick a mutant with PR=2^-10, and its corresponding phenotype in PRM 
myPR= -10
ddG1= Pheno_to_param_list[[1]](mywts, 2^myPR, 1) 
myPRM= log2(Forward(mywts, ddG1,0,0,0)[[2]])
myphes<- list(myPR, myPRM) 


ddF<- seq(-2, 5, by=0.1) # This is the folding energy, 71 different values 
ddBs<- list()
for (i in 1:2) {
  ddBs[[i]]= sapply(ddF, function(v){Pheno_to_param_pleiotropy_FDB(mywts,v,0, 2^myphes[[i]], i)})
  names(ddBs[[i]]) = ddF
  }
pleioSingles<- list(
  pr= data.frame(ddB= unlist(ddBs[[1]])), 
  prm= data.frame(ddB= unlist(ddBs[[2]]))
)

for (i in 1:2){
  pleioSingles[[i]]$ddF= as.numeric(rownames(pleioSingles[[i]]))
  cat(i) 
  print(length(rownames(pleioSingles[[i]])))
  print(unique(pleioSingles[[i]]$ddF))
}

pleioSingles[[2]][pleioSingles[[2]]$ddF< -2, 'ddF'] = as.numeric(specify_decimal(pleioSingles[[2]][pleioSingles[[2]]$ddF< -2, 'ddF']*0.1, 0))
pleioSingles[[2]][pleioSingles[[2]]$ddF> 10, 'ddF'] = as.numeric(specify_decimal(pleioSingles[[2]][pleioSingles[[2]]$ddF> 10, 'ddF']*0.1, 0))
pleioSingles[[2]]$ddF = as.numeric(specify_decimal(pleioSingles[[2]]$ddF, 1))

save(pleioSingles, file='pleiotropic_FB_single_list.RData')
######## how do they combine with other mutants? 
library(reshape2)
whichParam= c('f', 'b', 'd', 't', 'f,b')
doubs_with_pleio<- list(pr=list(), prm=list())
for (i in 1:2) {
  for(j in 1:4){
    p2<- sings[[i]][sings[[i]]$type1==whichParam[j], ]
    params1<- expand.grid(P1_ddG_f=pleioSingles[[i]]$ddF, P2_ddG= p2$param)
    params2<- expand.grid(P1_ddG_b=pleioSingles[[i]]$ddB, P2_pheno= p2$phen)
    doubs_with_pleio[[i]][[j]]= cbind(params1, params2) 
    doubs_with_pleio[[i]][[j]]$P2_which= whichParam[j] 
    doubs_with_pleio[[i]][[j]]$promoter = names(doubs_with_pleio)[i]
  }
  px<- doubs[[i]][doubs[[i]]$whichparam=='f,b',]
  p2= px[sample(nrow(px), 135), ]
  P2_pheno= expand.grid(P1_ddG_b=pleioSingles[[i]]$ddB, P2_pheno = p2$pheno)
  params1<- expand.grid(P1_ddG_f=pleioSingles[[i]]$ddF, P2_ddG_f= p2$param1)
  params2<- expand.grid(P1_ddG_b=pleioSingles[[i]]$ddB, P2_ddG_b= p2$param2)
  
  doubs_with_pleio[[i]][[5]]= cbind( params1, params2, P2_pheno) 
  doubs_with_pleio[[i]][[5]]$promoter = names(doubs_with_pleio)[i]
} 

for (i in 1:2) {
    
    doubs_with_pleio[[i]][[5]]= doubs_with_pleio[[i]][[5]][, c(1:4,6:7) ]

}

for (i in 1:2) {
  doubs_with_pleio[[i]][[1]]$ddG_f= doubs_with_pleio[[i]][[1]]$P1_ddG_f+ doubs_with_pleio[[i]][[1]]$P2_ddG
  doubs_with_pleio[[i]][[2]]$ddG_b= doubs_with_pleio[[i]][[2]]$P1_ddG_b+ doubs_with_pleio[[i]][[2]]$P2_ddG 
  doubs_with_pleio[[i]][[5]]$ddG_f= doubs_with_pleio[[i]][[5]]$P1_ddG_f+ doubs_with_pleio[[i]][[5]]$P2_ddG_f 
  doubs_with_pleio[[i]][[5]]$ddG_b= doubs_with_pleio[[i]][[5]]$P1_ddG_b+ doubs_with_pleio[[i]][[5]]$P2_ddG_b
} 

for (i in 1:2) {
  doubs_with_pleio[[i]][[1]]$doubP= apply(  doubs_with_pleio[[i]][[1]], 1, function(v) {Forward(8.352848e-07,as.numeric(v[7]), 0,  as.numeric(v[3]),0)[[i]]})
  doubs_with_pleio[[i]][[2]]$doubP= apply(  doubs_with_pleio[[i]][[2]], 1, function(v) {Forward(8.352848e-07,as.numeric(v[1]),0,as.numeric(v[7]),0)[[i]]})
  doubs_with_pleio[[i]][[3]]$doubP= apply(  doubs_with_pleio[[i]][[3]], 1, function(v) {Forward(8.352848e-07,as.numeric(v[1]), as.numeric(v[2]), as.numeric(v[3]), 0)[[i]]})
  doubs_with_pleio[[i]][[4]]$doubP= apply(  doubs_with_pleio[[i]][[4]], 1, function(v) {Forward(8.352848e-07,as.numeric(v[1]), 0, as.numeric(v[3]), as.numeric(v[2]))[[i]]})
  doubs_with_pleio[[i]][[5]]$doubP= apply(  doubs_with_pleio[[i]][[5]], 1, function(v) {Forward(8.352848e-07,as.numeric(v[7]), 0, as.numeric(v[8]),0)[[i]]})
  
  }


for (i in 1:2){
  doubs_with_pleio[[i]][[5]]$P2_pheno= log2( doubs_with_pleio[[i]][[5]]$P2_pheno)
   doubs_with_pleio[[i]][[5]]$P2_which= 'fb_pleio'
}

save(doubs_with_pleio, file='pleiotropic_with1_-10_-1.484815_with_one_not_list_log0_1_distribut.RData')


########
#########
############# Fitness landscape of the mutation-combinations 
myd<-seq(-2, 10, by=0.1) # 61 vals 
myds<- data.frame(ddg= rep(myd, 4), whichPar= c(rep(2, length(myd)), rep(3, length(myd)),rep(4, length(myd)),rep(5, length(myd))))

for (i in 1:length(rownames(myds))) {
  x<- pr_prm_Intramol_protein_withDNA_whichparam(mywts, myds[i,'whichPar'] ,myds[i, 'ddg'], 2, 0)
  myds[i,'pr'] = x[[1]] 
  myds[i,'prm'] = x[[2]]
}
save(myds, file='ddG_to_output.RData')

doubs_forward<- data.frame()

xx<- expand.grid(ddG1= myds$ddg,ddG2= myds$ddg )
x1<- expand.grid(type1= myds$whichPar,type2= myds$whichPar )
x2<- expand.grid(pheno1_pr= myds$pr,pheno2_pr= myds$pr)
x3<- expand.grid(pheno1_prm= myds$prm,pheno2_prm= myds$prm)

doubs_forward = cbind(xx, x1, x2,x3) 

print(head(doubs_forward))


doubs_forward$paramtyp= paste(doubs_forward$type1, doubs_forward$type2, sep=',')

rm(xx); rm(x1) ; rm(x2); rm(x3)
for (i in 1:length(rownames(doubs_forward))) {
  x=  pr_prm_Intramol_protein_withDNA_whichparam (
    mywts, as.numeric(doubs_forward[i,3]), as.numeric( doubs_forward[i,1]), 
    as.numeric(doubs_forward[i,4]),  as.numeric(doubs_forward[i,2])) 
  
  doubs_forward[i,'pr'] = x[[1]] 
  doubs_forward[i,'prm'] = x[[2]] 
  
}
x1<- doubs_forward[doubs_forward$paramtyp=='2,2', 'pr']
hist(log2(x1))
hist(log2(myds$pr))
x1<- myds[myds$whichPar=='2', 'pr']

save(doubs_forward, file='ddG_to_output_doubs.RData')



######## 
######## pleio +pleio example 
mydgs<- seq(-2, 5, by=0.2) # 36  
mydDs<- seq(-2, 5, by=0.2)  # 36

pleio3<- list(
  p_10=data.frame(),
  p_1= data.frame()
)

for (i in 1:2){
  pleio3[[i]]= data.frame(expand.grid(mydgs, mydDs)) 
}
myps<- c(2^-10, 2^-1)

for (i in 1:2){
  names(pleio3[[i]]) = c('ddF', 'ddD')  
   pleio3[[i]]=pleio3[[i]][pleio3[[i]]$ddF+pleio3[[i]]$ddD<5.8,  ] 
  # the condition is set to remove those that can not reach the given phentoype with whichever ddB based on the forward- modeling 
  pleio3[[i]]$phen= myps[i]
} 

mypars<- list(p_10=c(), p_1=c())

for (j in 1:length(myps)) {
  parms<- list()
  for (k in 1:length(rownames(pleio3[[j]]))) {
    parms[[k]]<- Pheno_to_param_pleiotropy_FDB(8.352848e-07,as.numeric(pleio3[[j]][k,1]), 
                                               as.numeric(pleio3[[j]][k,2]),myps[j], 1) 
  }
  names(parms)<- seq(1:length(rownames(pleio3[[j]])))
  mypars[[j]] = unlist(parms)
}

for (i in 1:2){
  pleio3[[i]]$ddB=NA
  for(j in 1:length(mypars[[i]])){ 
    pleio3[[i]][as.numeric(as.character(names(mypars[[i]][j]))), 'ddB'] = mypars[[i]][j]
  }
  # pleio3[[i]]=pleio3[[i]][!is.na(pleio3[[i]]$ddB), ]
}

for (i in 1:2){
  pleio3[[i]] = pleio3[[i]][!is.na(pleio3[[i]]$ddB), ]
  }

save(pleio3, file='Reverse_anyof3_param_log_even_minus10_minus1_as_prPheno.RData')
rm(mypars)

ABs<- cbind(
  expand.grid(ddF1=pleio3[[1]]$ddF, ddF2=pleio3[[2]]$ddF), 
  expand.grid(ddD1=pleio3[[1]]$ddD, ddD2=pleio3[[2]]$ddD), 
  expand.grid(ddB1=pleio3[[1]]$ddB, ddB2=pleio3[[2]]$ddB) 
) # 957,151

ABs$ddF= ABs$ddF1+ ABs$ddF2
ABs$ddD= ABs$ddD1+ ABs$ddD2
ABs$ddB= ABs$ddB1+ ABs$ddB2

ABs2<- ABs[sample(rownames(ABs), 90000), ] # 1/100 sampling, otherwise too big to plot
#### to start with the sampled ones
for (i in 1:length(rownames(ABs2))) {
  x1<- Forward(mywts, ABs2[i,'ddF'] , ABs2[i,'ddD'] ,ABs2[i,'ddB'], 0)
  ABs2[i, 'PR'] = x1[[1]] 
  ABs2[i, 'PRM'] =x1[[2]]
}
save(ABs2, file='log_even0.000977_0.5_as_prPheno_doub_sampled_90000.RData')
