
setwd('~/dir/')
###### Later we can update it. 
library(viridis)
library(ggplot2)
library(ggpubr)
library(reshape2)

########################## 1 Protein_with_another_molecule_for generalizability ##############
# P1 amount = Molecule amount 
# P1 stability = marginal - 1kcal/mol
# P1-P2 binding = -2 
### 1. ddG_folding -> PM & ### 2. ddG_binding -> PM 
# make mutations 
myppi<- function(pr2,deltadelta_G1, deltadelta_Gbind){ 
  # try pr2 = less than pr1 (0.5). same as pr1 (1). More than pr1 (2)   
  # deltag1_wt or deltag2_wt in a range of -3~ -1 is reasonable 
  # you may try deltadelta_G1 (binding energy) in a range of [-5, -1] 
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15 # this is 37 degrees!!! 
  deltG_bind_wt = -2 # this to update 
  
  pr1=1  
  deltag1_wt = - 1 # this to update
  deltag2_wt = deltag1_wt # this, may update
  
  deltaG1= deltag1_wt + deltadelta_G1
  deltaG2= deltag2_wt 
  
  deltaGb= deltG_bind_wt + deltadelta_Gbind
  
  k1 = exp(-deltaG1/(R*Temp)) # fraction of the folded protein = k1/(1+k1)
  k2= exp(-deltaG2/(R*Temp))
  k3= exp(-deltaGb/(R*Temp))
  
  a1= k3/((1+k1^-1)* (1+k2^-1))
  b1= -(pr1+pr2)*k3/((1+k1^-1)* (1+k2^-1))-1
  c1= k3*pr1*pr2/((1+k1^-1)* (1+k2^-1))
  
  comp2= (-b1-(b1^2-4*a1*c1)^0.5)/(2*a1)
  comp1= (-b1+(b1^2-4*a1*c1)^0.5)/(2*a1) 
  
  if (comp2>0 && comp2< min(pr1, pr2)) {
    sol=comp2
  } else if (comp1>0 && comp1<min(pr1, pr2)){
    sol=comp1
  } else {
    sol=NA
  }
  return(sol)
}

##### make the dataframe for mutations 
##### first for WT 
pr2<- c(0.25, 1, 4)
wts<-c()
wts<- sapply(pr2, function(v) {myppi(v,0,0)})

dds<- seq(-2,10, by=0.1) 
dat<- expand.grid(ddx= dds, ddy= dds)
myPleioSing<- list()
for(i in 1:3){
  myPleioSing[[i]]= dat 
  myPleioSing[[i]]$Compl= apply(myPleioSing[[i]], 1, function(v){myppi(pr2[i], v[1], v[2])})
  myPleioSing[[i]]$rel_Compl = myPleioSing[[i]]$Compl/ wts[i]
  myPleioSing[[i]]$LigandAmount= pr2[i] 
} # takes seconds to calculate 

myPleios<- rbind(myPleioSing[[1]], myPleioSing[[2]], myPleioSing[[3]])
save(myPleios, file='ppi_temp.RData')

#####
##### plot

ggplot(data=myPleios[myPleios$ddx<=7 & myPleios$ddy<=7 & myPleios$LigandAmount==1,  ],aes(x= ddx, y= ddy, fill= log2(Compl))) + geom_tile()+ 
  theme_classic() + labs(x= '∆∆Gf', y='∆∆Gb') + scale_fill_viridis( name= 'Complex(abs.log2)') +   
  geom_contour(aes(x = ddx, y = ddy, z =log2(Compl),colour = ..level..),binwidth=1 , col='gray')+ 
  geom_point(aes(x=0, y=0), shape= 4,size=2 , col='white') 


fols<- ggplot(data= myPleios[myPleios$ddy==0 & myPleios$LigandAmount==1,], aes(x=ddx, y= log2(Compl))) +  
  geom_line(size=1) + 
  theme_classic() + geom_vline(xintercept = 0, col='gray', lty=2)+ 
  labs(x='∆∆Gf', y='Complex(log2)')
binds<- ggplot(data= myPleios[myPleios$ddx==0& myPleios$LigandAmount==1,], aes(x=ddy, y= log2(Compl))) +
  geom_line(size=1) + 
  theme_classic() + geom_vline(xintercept = 0, col='gray', lty=2)+ 
  labs(x='∆∆Gb', y='Complex(log2)')

ggarrange(fols, binds,nrow=1, ncol=2)

xx<- myPleios[myPleios$ddx<=7 & myPleios$ddy<=7 & myPleios$LigandAmount==1,  ]
yy= xx[log2(xx$Compl)>= -1.01 & log2(xx$Compl)<= -0.99, c('ddx', 'ddy')]
zz= cbind (expand.grid(yy$ddx, yy$ddx), expand.grid(yy$ddy, yy$ddy)) 
names(zz)= c('ddx1', 'ddx2', 'ddy1', 'ddy2')

ggplot(data=xx,aes(x= ddx, y= ddy)) + geom_tile(fill= 'white')+ 
  theme_classic() + labs(x= '∆∆Gf', y='∆∆Gb') + 
  geom_contour(aes(x = ddx, y = ddy, z =log2(Compl),colour = ..level..),binwidth=1 , col='gray')+ 
  geom_point(data= zz, aes(x= ddx1+ddx2, y= ddy1+ddy2), shape= 4,size=2 , col='red',  alpha=0.5)  
  



########################## 2 Protein_protein_dimer_for_generalizability ######### ##############
### 1. ddG_folding -> PP, both mutations are shared , one-copy in each
### 2. ddG_binding -> PP,both mutations are shared , one-copy in each
### 3. ddG_folding + ddG_binding -> fitness, PP 



###### PPI dimerization 
myppi_dim<- function(pr2,deltadelta_G1, deltadelta_Gbind){ 
  # try pr2 = less than pr1 (0.5). same as pr1 (1). More than pr1 (2)   
  # deltag1_wt or deltag2_wt in a range of -3~ -1 is reasonable 
  # you may try deltadelta_G1 (binding energy) in a range of [-5, -1] 
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15 # this is 37 degrees!!! 
  deltG_bind_wt = -2 # this to update 
  
  pr1=0.5* pr2 
  pr2= pr1
  deltag1_wt = - 1 # this to update
  deltag2_wt = deltag1_wt # this, may update
  
  deltaG1= deltag1_wt + deltadelta_G1
  deltaG2= deltag2_wt + deltadelta_G1
  
  deltaGb= deltG_bind_wt + deltadelta_Gbind
  
  k1 = exp(-deltaG1/(R*Temp)) # fraction of the folded protein = k1/(1+k1)
  k2= exp(-deltaG2/(R*Temp))
  k3= exp(-deltaGb/(R*Temp))
  
  a1= k3/((1+k1^-1)* (1+k2^-1))
  b1= -(pr1+pr2)*k3/((1+k1^-1)* (1+k2^-1))-1
  c1= k3*pr1*pr2/((1+k1^-1)* (1+k2^-1))
  
  comp2= (-b1-(b1^2-4*a1*c1)^0.5)/(2*a1)
  comp1= (-b1+(b1^2-4*a1*c1)^0.5)/(2*a1) 
  
  if (comp2>0 && comp2< min(pr1, pr2)) {
    sol=comp2
  } else if (comp1>0 && comp1<min(pr1, pr2)){
    sol=comp1
  } else {
    sol=NA
  }
  return(sol)
}

# make mutations 

##### make the dataframe for mutations 
##### first for WT 
pr<- c(1, 2, 4)
wts<-c()
wts<- sapply(pr2, function(v) {myppi_dim(v,0,0)})

dat2<- expand.grid(ddx= dds, ddy= dds)
myPleioSing2<- list()
for(i in 1:3){
  myPleioSing2[[i]]= dat 
  myPleioSing2[[i]]$Compl= apply(myPleioSing2[[i]], 1, function(v){myppi_dim(pr2[i], v[1], v[2])})
  myPleioSing2[[i]]$rel_Compl = myPleioSing2[[i]]$Compl/ wts[i]
  myPleioSing2[[i]]$Total_amount= pr[i] 
} # takes seconds to calculate 

myPleios2<- rbind(myPleioSing2[[1]], myPleioSing2[[2]], myPleioSing2[[3]])
save(myPleios2, file='ppi_dimerization_temp.RData')
#####
##### plot

ggplot(data=myPleios2[myPleios2$ddx<=5 & myPleios2$ddy<=5 & myPleios2$Total_amount==2,  ],aes(x= ddx, y= ddy, fill= log2(Compl))) + geom_tile()+ 
  theme_classic() + labs(x= '∆∆Gf', y='∆∆Gb') + scale_fill_viridis( name= 'Complex(abs.log2)') +   
  geom_contour(aes(x = ddx, y = ddy, z =log2(Compl),colour = ..level..),binwidth=1 , col='gray')+ 
  geom_point(aes(x=0, y=0), shape= 4,size=2 , col='white') 


fols<- ggplot(data= myPleios2[myPleios2$ddy==0 & myPleios2$Total_amount==2 & myPleios2$ddx<=7,], aes(x=ddx, y= log2(Compl))) +  
  geom_line(size=1) + 
  theme_classic() + geom_vline(xintercept = 0, col='gray', lty=2)+ 
  labs(x='∆∆Gf', y='Complex(log2)')
binds<- ggplot(data= myPleios2[myPleios2$ddx==0& myPleios2$Total_amount==2& myPleios2$ddy<=7,], aes(x=ddy, y= log2(Compl))) +
  geom_line(size=1) + 
  theme_classic() + geom_vline(xintercept = 0, col='gray', lty=2)+ 
  labs(x='∆∆Gb', y='Complex(log2)')

ggarrange(fols, binds,nrow=1, ncol=2)

xx<- myPleios2[myPleios2$ddx<=5 & myPleios2$ddy<=5 & myPleios2$Total_amount==2,  ]
yy= xx[log2(xx$Compl)>= -2.01 & log2(xx$Compl)<= -1.99, c('ddx', 'ddy')]
zz= cbind (expand.grid(yy$ddx, yy$ddx), expand.grid(yy$ddy, yy$ddy)) 
names(zz)= c('ddx1', 'ddx2', 'ddy1', 'ddy2')

ggplot(data=xx,aes(x= ddx, y= ddy)) + geom_tile(fill= 'white')+ 
  theme_classic() + labs(x= '∆∆Gf', y='∆∆Gb') + 
  geom_contour(aes(x = ddx, y = ddy, z =log2(Compl),colour = ..level..),binwidth=1 , col='gray')+ 
  geom_point(data= zz, aes(x= ddx1+ddx2, y= ddy1+ddy2), shape= 4,size=2 , col='red',  alpha=0.5)  



#