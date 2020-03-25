setwd('~/dir/')
###### 
########## for plotting
library(ggplot2) 
library(ggpubr)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # first gray, second near yellow, third near blue
mygreenPalette<- c('#C3F780', '#8BC34A', '#76A63F', '#648C35', '#486627', '#364D1D', '#243313')
mygreyPalette<- c('#EEEEEE', '#DDDDDD', '#CCCCCC', '#BBBBBB', '#AAAAAA', '#999999', '#777777', '#555555', '#333333','#111111')
myredPalette<- c('#FFEBEE', '#FFCDD2', '#EF9A9A', '#E57373', '#EF5350', '#F44336', '#E53935', '#D32F2F', '#C62828', '#B71C1C')
mypinkPalette<- c('#fef2f9', '#fcdef0', '#fac6e5', '#ffb0e1',  '#ff95d6',  '#ff80ce',  '#e76eb1',  '#cd5a91', '#b44772',  '#96304c')

library(viridis)
##### to plot dose responses 

####### For Fig 1 for WT 
load('Dose_response_WT.RData') 


### linear scale

ggplot(data= dose_response) + 
  geom_line(aes(x=Total_protein, y=pr)) + scale_x_log10() + 
  geom_vline(xintercept = mywts, col='red') + 
  theme_classic()  + 
  geom_line(aes(x=Total_protein, y=prm), lty=2) + labs(x= 'CI concentration (M)', y= 'Gene expression') # 4*3
### log scale
ggplot(data= dose_response) + 
  geom_line(aes(x=Total_protein, y=log2(pr))) + scale_x_log10() + 
  geom_vline(xintercept = mywts, col='red') + 
  theme_classic()  + 
  geom_line(aes(x=Total_protein, y=log2(prm)), lty=2) + labs(x= 'CI concentration (M)', y= 'Gene expression (log2)') # 4*3

rm(dose_response) 

####### Fig 1 and Fig S1. The relationships between parameters and the outputs

load ('ddG_to_output.RData')
ggplot(data= myds)+ geom_line(aes(x=  ddg, y= pr), col='red') + 
  geom_line(aes(x=  ddg, y= prm), col='red', lty=2) + 
  theme_classic() +
  geom_vline(xintercept = 0, lty=3) + facet_wrap(whichPar~. ) + 
  ylim(0, 1) + labs(x='??????G (kcal/mol)', y='Gene expression') + xlim(-2, 10)

p1<- ggplot(data= myds)+ geom_line(aes(x=  ddg, y= log2(pr))) + 
  geom_line(aes(x=  ddg, y= prm), col='red', lty=2) + 
  theme_classic() +
  geom_vline(xintercept = 0, lty=2, col='gray') + facet_grid(.~whichPar) + 
  xlim(-2, 10) + ylim(-20, 0)+ labs(x='??????G (kcal/mol)', y='Gene expression (log2)')

p2<- ggplot(data= myds)+   geom_line(aes(x=  ddg, y= log2(prm))) + 
  theme_classic() +
  geom_vline(xintercept = 0, lty=2, col='gray') + facet_grid(.~whichPar) + 
  xlim(-2, 10) + labs(x='??????G (kcal/mol)', y='Gene expression (log2)')

ggarrange(p1, p2, nrow=2, ncol=1) # 
rm(myds); rm(p1); rm(p2); rm(myd)
rm(pr); rm(prm) ; rm(x)

########## Fig. S1 
########## Relationships between WT protein & phenotype relationship in each set of ∆∆G changes. 
load('starting_from_free_dimer.RData') 

out11<- ggplot(data= fre2[fre2$ddG2==0 & fre2$ddG3==0 & fre2$ddG4==0 & fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG1') + 
  geom_line(aes(x=Total_protein, y=pr, col= as.factor(ddG1)), size=0.7) + scale_x_log10() + 
   geom_vline(xintercept = mywts,  col='red') + labs(x='CI concentation (M)', y= 'Gene expression') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=prm, col= as.factor(ddG1)), lty=2, size=0.7) 

out22<- ggplot(data= fre2[fre2$ddG1==0 & fre2$ddG3==0 & fre2$ddG4==0& fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG2') + 
  geom_line(aes(x=Total_protein, y=pr, col= as.factor(ddG2)), size=0.7) + scale_x_log10() + 
   geom_vline(xintercept = mywts, col='red') + labs(x='CI concentation (M)', y= 'Gene expression') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=prm, col= as.factor(ddG2)), lty=2, size=0.7)

out33<-ggplot(data= fre2[fre2$ddG1==0 & fre2$ddG2==0 & fre2$ddG4==0& fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG3') + 
  geom_line(aes(x=Total_protein, y=pr, col= as.factor(ddG3)), size=0.7) + scale_x_log10() + labs(x='CI concentation (M)', y= 'Gene expression') + 
  geom_vline(xintercept = mywts, col='red') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=prm, col= as.factor(ddG3)), lty=2, size=0.7)

out44<-ggplot(data= fre2[fre2$ddG1==0 & fre2$ddG2==0 & fre2$ddG3==0& fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG4') + 
  geom_line(aes(x=Total_protein, y=pr, col= as.factor(ddG4)), size=0.7) + scale_x_log10() + labs(x='CI concentation (M)', y= 'Gene expression') + 
  geom_vline(xintercept = mywts, col='red') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=prm, col= as.factor(ddG4)), lty=2, size=0.7)


ggarrange(
          out11, out22, out33, out44, 
          nrow=2, ncol=2, common.legend = T)

rm(fre2)
###### y axis in log 2? 
out11<- ggplot(data= fre2[fre2$ddG2==0 & fre2$ddG3==0 & fre2$ddG4==0 & fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG1') + 
  geom_line(aes(x=Total_protein, y=log2(pr), col= as.factor(ddG1)), size=0.7) + scale_x_log10() + 
  geom_vline(xintercept = mywts,  col='red') + labs(x='CI concentation (M)', y= 'Gene expression') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=log2(prm), col= as.factor(ddG1)), lty=2, size=0.7) 

out22<- ggplot(data= fre2[fre2$ddG1==0 & fre2$ddG3==0 & fre2$ddG4==0& fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG2') + 
  geom_line(aes(x=Total_protein, y=log2(pr), col= as.factor(ddG2)), size=0.7) + scale_x_log10() + 
  geom_vline(xintercept = mywts, col='red') + labs(x='CI concentation (M)', y= 'Gene expression') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=log2(prm), col= as.factor(ddG2)), lty=2, size=0.7)

out33<-ggplot(data= fre2[fre2$ddG1==0 & fre2$ddG2==0 & fre2$ddG4==0& fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG3') + 
  geom_line(aes(x=Total_protein, y=log2(pr), col= as.factor(ddG3)), size=0.7) + scale_x_log10() + labs(x='CI concentation (M)', y= 'Gene expression') + 
  geom_vline(xintercept = mywts, col='red') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=log2(prm), col= as.factor(ddG3)), lty=2, size=0.7)

out44<-ggplot(data= fre2[fre2$ddG1==0 & fre2$ddG2==0 & fre2$ddG3==0& fre2$Total_protein>=10^-11 & fre2$Total_protein<=0.1, ]) + labs(title='ddG4') + 
  geom_line(aes(x=Total_protein, y=log2(pr), col= as.factor(ddG4)), size=0.7) + scale_x_log10() + labs(x='CI concentation (M)', y= 'Gene expression') + 
  geom_vline(xintercept = mywts, col='red') + 
  scale_color_manual(values= c(cbPalette[2], 'black', cbPalette[3], cbPalette[4])) + theme_classic()  + 
  geom_line(aes(x=Total_protein, y=log2(prm), col= as.factor(ddG4)), lty=2, size=0.7)
ggarrange(
  out11, out22, out33, out44, 
  nrow=2, ncol=2, common.legend = T)
rm(out11); rm(out22);rm(out33); rm(out44);rm(fre2) 


#### Sanity check for the forward & reverse functions. 
load('Even_space_4_param_singles_0_1_logeven_20200310.RData') # generated by the reverse function 
load('Even_space_4_param_doubs0_1_logeven_20200310.RData')  # generated by the forward function 
####### check whether the ∆∆G - phenotype relationship from the two functions give the same relationsihp 
####  linear first

pr= ggplot(data= doubs[[1]][doubs[[1]]$paramtyp=='same', ]) + 
  geom_line( aes(x= param1+param2, y= pheno), size=2) + 
  geom_point(data= sings[[1]], aes(x= param, y= 2^phen), col='red', size=1,  shape=19) + 
  facet_grid(.~type1) + theme_classic() + labs(x='∆∆G (kcal/mol)', y='Gene expression') + 
  geom_vline(xintercept = 0, lty=2) + xlim(-2, 10)

prm= ggplot(data= doubs[[2]][doubs[[2]]$paramtyp=='same', ]) + 
  geom_line( aes(x= param1+param2, y= pheno), size=2) + 
  geom_point(data= sings[[2]], aes(x= param, y= 2^phen), col='red', size=1,  shape=19) + 
  facet_grid(.~type1) + theme_classic() + labs(x='∆∆G (kcal/mol)', y='Gene expression') + 
  geom_vline(xintercept = 0, lty=2) + xlim(-2, 10)
ggarrange(pr, prm, nrow=2, ncol=1)

##### log 2 cale 

pr= ggplot(data= doubs[[1]][doubs[[1]]$paramtyp=='same', ]) + 
  geom_line( aes(x= param1+param2, y= log2(pheno)), size=2) + 
  geom_point(data= sings[[1]], aes(x= param, y= phen), col='red', size=1,  shape=19) + 
  facet_grid(.~type1) + theme_classic() + labs(x='∆∆G (kcal/mol)', y='Gene expression (log2)') + 
  geom_vline(xintercept = 0, lty=2) + xlim(-2, 10)

prm= ggplot(data= doubs[[2]][doubs[[2]]$paramtyp=='same', ]) + 
  geom_line( aes(x= param1+param2, y= log2(pheno)), size=2) + 
  geom_point(data= sings[[2]], aes(x= param, y= phen), col='red', size=1, shape=19) + 
  facet_grid(.~type1) + theme_classic() + labs(x='∆∆G (kcal/mol)', y='Gene expression (log2)') + 
  geom_vline(xintercept = 0, lty=2) + xlim(-2, 10)
ggarrange(pr, prm, nrow=2, ncol=1)
rm(pr); rm(prm)
dev.off()
##############
###############
################

###### with PR and PRM with tetramerization mutants, generate the fitness landscape 
##### 2.1 tile plots of double mutants' phenotypes, when there is no ambiguity 
# Fig.  2 and Fig. S2 
ggplot(data=doubs[[1]][doubs[[1]]$paramtyp=='same', ]) + geom_tile(aes(x= pheno1, y= pheno2, fill= phen_log2))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B', title='PR') + scale_fill_viridis( name= 'Mutation AB') + 
  geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept =   log2(myw_outs[[1]]), lty=2)+ facet_wrap(.~whichparam) +
   xlim(-13.5, 0) + ylim (-13.5, 0) # 5.6*4.7

ggplot(data=doubs[[1]][doubs[[1]]$paramtyp=='same', ]) + geom_tile(aes(x= pheno1, y= pheno2, fill= phen_log2))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B', title='PR') + scale_fill_viridis( name= 'Mutation AB') + 
  geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept =   log2(myw_outs[[1]]), lty=2)+ facet_wrap(.~whichparam)  + 
 geom_contour(aes(x= pheno1, y= pheno2, z =phen_log2,colour = ..level..),colour = "black",alpha=0.5 ,binwidth=1 )+ 
  xlim(-13.5, 0) + ylim (-13.5, 0) #  5.6*4.7

ggplot(data=doubs[[2]][doubs[[2]]$whichparam=='t,t', ]) + geom_tile(aes(x= pheno1, y= pheno2, fill= phen_log2))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B', title='PRM') + scale_fill_viridis( name= 'Mutation AB') + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept =   log2(myw_outs[[2]]), lty=2) +    xlim(-13.5, 0) + ylim (-13.5, 0) # 3.7*2.59

ggplot(data=doubs[[2]][doubs[[2]]$whichparam=='t,t', ]) + geom_tile(aes(x= pheno1, y= pheno2, fill= phen_log2))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B', title='PRM') + scale_fill_viridis( name= 'Mutation AB') + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept =   log2(myw_outs[[2]]), lty=2) +    xlim(-13.5, 0) + ylim (-13.5, 0)+
  geom_contour(aes(x= pheno1, y= pheno2, z =phen_log2,colour = ..level..),colour = "gray",binwidth=1 )

ggplot(data=doubs[[1]][doubs[[1]]$paramtyp!='same', ]) + geom_tile(aes(x= pheno1, y= pheno2, fill= phen_log2))+ 
    theme_classic() + labs(x= 'Mutation A', y='Mutation B', title='PR') + scale_fill_viridis( name= 'Mutation AB') + 
    geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept =   log2(myw_outs[[1]]), lty=2)+ facet_wrap(type1~type2)  
  
ggplot(data=doubs[[1]][doubs[[1]]$paramtyp!='same', ]) + geom_tile(aes(x= pheno1, y= pheno2, fill= phen_log2))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B', title='PR') + scale_fill_viridis( name= 'Mutation AB') + 
  geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept =   log2(myw_outs[[1]]), lty=2)+ facet_wrap(type1~type2) +  
 geom_contour(aes(x= pheno1, y= pheno2, z =phen_log2 ,colour = ..level..),colour = "black",alpha=0.5 ,binwidth=1 )

dev.off()  

#######
##### ambiguity plots Fig.2, and Fig S2,   
############ ambiguity with PRM, with two different known parameters
load('biochem_ambiguity_doubs_twodiff_params_lists_prm2_with_log2_0.01_precision.RData') 

ambiPlot<- list()

ambiPlot[[1]] <- ggplot(data=mycombs_inrow) + geom_tile(aes(x= p1, y= p2, fill=as.factor(ambiNum)))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_manual(values= cbPalette[c(1,2,4)]) + 
  labs(fill='') + facet_grid(type1~type2) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept =  log2(myw_outs[[2]]), lty=2) + 
  geom_abline(intercept = 0, slope = 1, col='gray')


ambiPlot[[2]] <- ggplot(data=mycombs_inrow) + geom_tile(aes(x= p1, y= p2, fill= ambiDiff))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_gradientn(colours = c(mygreyPalette[2],mypinkPalette[c(2,4, 6)] ,myredPalette[c(8,10)]))  + 
  labs(fill='') + facet_grid(type1~type2) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept =  log2(myw_outs[[2]]), lty=2) + 
  geom_abline(intercept = 0, slope = 1, col='gray')


ggarrange(ambiPlot[[1]] ,ambiPlot[[2]] , nrow=1, ncol=2)  # 13*5.5

rm(ambiPlot) ; rm(mycombs_inrow)

###### Fig. 3 and Fig S3
load('biochem_ambiguity_4_param_doub_01_20200324_with_precision0.01_lists.RData') 

#### when allow biochemical ones are two.
### 10*4 
# PR 2s
ambiPlot<- list()
ambiPlot[[1]] <- ggplot(data=biochem_ambi[[3]][[1]][biochem_ambi[[3]][[1]]$allowed_biochem==2, ]) + geom_tile(aes(x= p1, y= p2, fill=as.factor(ambiNum)))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_manual(values= c(cbPalette[1:4], mypinkPalette)) + 
  labs(fill='') + facet_grid(. ~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[1]]), lty=2)

ambiPlot[[2]] <- ggplot(data=biochem_ambi[[3]][[1]][ biochem_ambi[[3]][[1]]$allowed_biochem==2, ]) + geom_tile(aes(x= p1, y= p2, fill= ambiDiff))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_gradientn (colours = c(mygreyPalette[2],mypinkPalette[c(2,4, 6)] ,myredPalette[c(8,10)]))  + 
  labs(fill='') + facet_grid(.~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[1]]), lty=2)
ggarrange(ambiPlot[[1]], ambiPlot[[2]], nrow=2, ncol=1) # 10*4

######### PRM 2s together
ambiPlot[[1]] <- ggplot(data=biochem_ambi[[3]][[2]][biochem_ambi[[3]][[2]]$allowed_biochem==2, ]) + geom_tile(aes(x= p1, y= p2, fill=as.factor(ambiNum)))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_manual(values= c(cbPalette[1:4], mypinkPalette[2:10], myredPalette[c(5:10)])) + 
  labs(fill='') + facet_grid(. ~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[2]]), lty=2)


ambiPlot[[2]] <- ggplot(data=biochem_ambi[[3]][[2]][biochem_ambi[[3]][[2]]$allowed_biochem==2, ]) + geom_tile(aes(x= p1, y= p2, fill= ambiDiff))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_gradientn (colours = c(mygreyPalette[2],mypinkPalette[c(2,4,6)] ,myredPalette[c(8,10)]))  + 
  labs(fill='') + facet_grid(.~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[2]]), lty=2)
ggarrange(ambiPlot[[1]], ambiPlot[[2]], nrow=2, ncol=1) # 10*4

######### PRM 1s together
ambiPlot[[1]] <- ggplot(data=biochem_ambi[[3]][[2]][biochem_ambi[[3]][[2]]$allowed_biochem==1, ]) + geom_tile(aes(x= p1, y= p2, fill=as.factor(ambiNum)))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_manual(values= c(cbPalette[1:4], mypinkPalette, myredPalette[c(5:10)])) + 
  labs(fill='') + facet_grid(. ~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[2]]), lty=2) #8*4.7


ambiPlot[[2]] <- ggplot(data=biochem_ambi[[3]][[2]][ biochem_ambi[[3]][[2]]$allowed_biochem==1, ]) + geom_tile(aes(x= p1, y= p2, fill= ambiDiff))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_gradientn (colours = c(mygreyPalette[2],mypinkPalette[c(2,4,6)] ,myredPalette[c(8,10)]))  + 
  labs(fill='') + facet_grid(.~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[2]]), lty=2)
ggarrange(ambiPlot[[1]], ambiPlot[[2]], nrow=2, ncol=1) # 10*4

###### PR 3 and 4s 
ambiPlot<- list()
ambiPlot[[1]] <- ggplot(data=biochem_ambi[[3]][[1]][biochem_ambi[[3]][[1]]$allowed_biochem %in% c(3,4), ]) + geom_tile(aes(x= p1, y= p2, fill=as.factor(ambiNum)))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_manual(values= c(cbPalette[1:4], mypinkPalette[c(2:10)])) + 
  labs(fill='') + facet_grid(. ~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[1]]), lty=2)

ambiPlot[[2]] <- ggplot(data=biochem_ambi[[3]][[1]][ biochem_ambi[[3]][[1]]$allowed_biochem %in% c(3,4), ]) + geom_tile(aes(x= p1, y= p2, fill= ambiDiff))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_gradientn (colours = c(mygreyPalette[2],mypinkPalette[c(2,4, 6)] ,myredPalette[c(8,10)]))  + 
  labs(fill='') + facet_grid(.~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[1]]), lty=2)
ggarrange(ambiPlot[[1]], ambiPlot[[2]], nrow=2, ncol=1) # 8*4

######### PRM 3, 4s together
ambiPlot[[1]] <- ggplot(data=biochem_ambi[[3]][[2]][biochem_ambi[[3]][[2]]$allowed_biochem %in% c(3,4), ]) + geom_tile(aes(x= p1, y= p2, fill=as.factor(ambiNum)))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_manual(values= c(cbPalette[c(2, 3,4)], mypinkPalette[2:10], myredPalette[c(4:10)])) + 
  labs(fill='') + facet_grid(. ~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[2]]), lty=2)


ambiPlot[[2]] <- ggplot(data=biochem_ambi[[3]][[2]][ biochem_ambi[[3]][[2]]$allowed_biochem %in% c(3,4), ]) + geom_tile(aes(x= p1, y= p2, fill= ambiDiff))+ 
  theme_classic() + labs(x= 'Mutation A', y='Mutation B') + scale_fill_gradientn (colours = c(mypinkPalette[c(2,4,6)] ,myredPalette[c(8,9,10)]))  + 
  labs(fill='') + facet_grid(.~which_biochem) + 
  geom_vline(xintercept = log2(myw_outs[[2]]), lty=2) + geom_hline(yintercept = log2(myw_outs[[2]]), lty=2)
ggarrange(ambiPlot[[1]], ambiPlot[[2]], nrow=2, ncol=1) # 8*4
rm(ambiPlot)

##### Example plots for Fig. 2 and 3 panels 
########
############# ####### examples of how two mutations combine, with ambiguity 
### Fig. 1, 2. Check 3 mutations i.e. 
myphenos<- c(-10, -5, -1)

tmp<- data.frame(doubs[[2]][doubs[[2]]$paramtyp=='same' & doubs[[2]]$pheno1%in% myphenos, ])
for (i in 1:length(rownames(tmp))) { 
  tmp[i, 'ambiNum'] = biochem_ambi[[3]][[2]][str_sub(biochem_ambi[[3]][[2]]$which_biochem,1,1)==tmp[i,'type1'] & 
                                               biochem_ambi[[3]][[2]]$allowed_biochem==1 & 
                                               biochem_ambi[[3]][[2]]$p1== tmp[i,'pheno1'] & 
                                               biochem_ambi[[3]][[2]]$p2== tmp[i,'pheno2'], 'ambiNum']
  
}

ggplot(data= tmp) + 
  geom_point(aes(x= pheno2, y=phen_log2, col= as.factor(ambiNum)), shape=4, size=0.5) + 
  facet_grid(type1~ as.factor(pheno1)) + theme_classic() + scale_colour_manual(values=c(cbPalette[c(1:4)],mypinkPalette[2])) # 4*6

tmp2<- list(pr=list(), 
            prm=list()) 

mysets<- list(
  FD= c('d,f','f,d', 'd,d', 'f,f'), 
  FB=c('b,f','f,b', 'b,b', 'f,f') , 
  FT=c('t,f','f,t', 't,t', 'f,f'), 
  DB=c('b,d','d,b', 'd,d', 'b,b'), 
  DT=c('t,d','d,t', 't,t', 'd,d'), 
  BT=c('t,b','b,t', 't,t', 'b,b')
)

mysecond_pheno<- list(
  unique(doubs[[1]]$pheno2), 
  unique(doubs[[2]]$pheno2)
) 

for (i in 1:2) {

  for (j in 1:6){
    tmp2[[i]][[j]]= data.frame(doubs[[i]][doubs[[i]]$pheno1%in% myphenos &
                                            doubs[[i]]$whichparam %in%  mysets[[j]] , ])
    
    
    for (k in myphenos){
      for (l in mysecond_pheno[[i]]) {
        
        tmp2[[i]][[j]][tmp2[[i]][[j]]$pheno1==k & tmp2[[i]][[j]]$pheno2== l, 'ambiNum']=
          length(unique(tmp2[[i]][[j]][tmp2[[i]][[j]]$pheno1==k & tmp2[[i]][[j]]$pheno2== l ,'pheno'])) 
      }
    }

    tmp2[[i]][[j]]$allowed_2_group = names(mysets)[j]
        }
  
}

for (i in 1:2) {
  for (j in 1:6){
    tmp2[[i]][[j]]$allowed_2_group = names(mysets)[j]
  }
}

ambiguity_examples<- list()
for(i in 1:2){
  ambiguity_examples[[i]]= rbind(tmp2[[i]][[1]], tmp2[[i]][[2]],tmp2[[i]][[3]],tmp2[[i]][[4]],tmp2[[i]][[5]],tmp2[[i]][[6]])
  ambiguity_examples[[i]]$allowed_2_group= factor(ambiguity_examples[[i]]$allowed_2_group, 
                                                  levels= c('FD', 'FB', 'FT', 'DB', 'DT', 'BT'))
}
###
examplPlot<- list()
for (i in 1:2){
    examplPlot[[i]]=ggplot(data= ambiguity_examples[[i]]) + 
      geom_point(aes(x= pheno2, y=phen_log2, col= as.factor(ambiNum)), shape=4, size=0.5) + 
      facet_grid(allowed_2_group~ as.factor(pheno1)) + theme_classic() + 
      scale_colour_manual(values=c(cbPalette[c(1:4)],mypinkPalette[3:10])) + 
      labs(x='', y= '') # 6*8 

} # 
examplPlot[[1]]
examplPlot[[2]]





rm(examplPlot)
dev.off()


##########  
##### 2.4 Ambiugity jitter plots  Fig. 3 and S3

mygs<- list() # ambiguity Number
mygs2<- list() # ambiguity range
dev.off()
for (i in 1:2) {
  mygs[[i]] = ggplot(data=biochem_ambi[[4]][[i]], aes(x=as.factor(bioche_no), y= Max_ambiNum, color=as.factor(invol_f))) + 
    labs(x="# allowed biochemical changes", y= "Max_ambiguities", title= names(biochem_ambi[[3]])[i]) +theme_classic() + 
    stat_summary(fun="mean",geom="crossbar", 
                 mapping=aes(ymin=..y.., ymax=..y..), width=0.5)+ 
    geom_jitter(alpha=0.4 , size=3, position =position_jitterdodge()) + 
    scale_color_manual(values= c('black', 'red'))
}

### values 
for (i in 1:2) {
  mygs2[[i]] = ggplot(data=biochem_ambi[[4]][[i]], aes(x=as.factor(bioche_no), y= Max_ambiDiff, color=as.factor(invol_f))) + 
    labs(x="# allowed biochemical changes", y= "Max_ambiguiDiff", title= names(biochem_ambi[[3]])[i]) +theme_classic() + 
    stat_summary(fun="mean",geom="crossbar", 
                 mapping=aes(ymin=..y.., ymax=..y..), width=0.5)+ 
    geom_jitter(alpha=0.4 , size=3,  width = 0.1) + 
    scale_color_manual(values= c('black', 'red'))
}
ggarrange(mygs[[1]], mygs2[[1]],mygs[[2]], mygs2[[2]], nrow=1, ncol=4, common.legend = T) #12*3
############## what about median? 
for (i in 1:2) {
  mygs[[i]] = ggplot(data=biochem_ambi[[4]][[i]], aes(x=as.factor(bioche_no), y= Med_ambiDiff, color=as.factor(invol_f))) + 
    labs(x="# allowed biochemical changes", y= "Med_ambiguities", title= names(biochem_ambi[[3]])[i]) +theme_classic() + 
    stat_summary(fun="mean",geom="crossbar", 
                 mapping=aes(ymin=..y.., ymax=..y..), width=0.5)+ 
    geom_jitter(alpha=0.4 , size=3, position =position_jitterdodge()) + 
    scale_color_manual(values= c('black', 'red'))
}

### values 
for (i in 1:2) {
  mygs2[[i]] = ggplot(data=biochem_ambi[[4]][[i]], aes(x=as.factor(bioche_no), y= Med_ambiDiff, color=as.factor(invol_f))) + 
    labs(x="# allowed biochemical changes", y= "Med_ambiguiDiff", title= names(biochem_ambi[[3]])[i]) +theme_classic() + 
    stat_summary(fun="mean",geom="crossbar", 
                 mapping=aes(ymin=..y.., ymax=..y..), width=0.5)+ 
    geom_jitter(alpha=0.4 , size=3,  width = 0.1) + 
    scale_color_manual(values= c('black', 'red'))
}
ggarrange(mygs[[1]], mygs2[[1]],mygs[[2]], mygs2[[2]], nrow=1, ncol=4, common.legend = T)

rm(i);rm(j); rm(k); rm(l) ; 

rm(tmp) ; rm(tmp2) ; rm(to_check) ;  
rm(mygs); rm(mysecond_pheno) ; rm(samplphens)
rm(mygs2)

#######
####### Fig. 4
####### to plot how one pleiotropic mutation combine with other types of mutants

load('pleiotropic_with1_-10_-1.484815_with_one_not_list_log0_1_distribut.RData') 
dobs<- list()
dobs[[1]] = rbind(doubs_with_pleio[[1]][[1]][,c('P2_pheno','promoter','P2_which', 'doubP')], 
                  doubs_with_pleio[[1]][[2]][,c('P2_pheno','promoter', 'P2_which', 'doubP')])
dobs[[2]] = rbind(doubs_with_pleio[[2]][[1]][,c('P2_pheno', 'promoter','P2_which', 'doubP')], 
                  doubs_with_pleio[[2]][[2]][,c('P2_pheno', 'promoter','P2_which', 'doubP')])
for (i in 1:2) {
  for (j in 3:4) {
    dobs[[i]]= rbind(dobs[[i]], doubs_with_pleio[[i]][[j]][,c('P2_pheno', 'promoter','P2_which', 'doubP')])
  } 
  dobs[[i]]$P2_which= factor(dobs[[i]]$P2_which, levels= c('f','d','b','t'))
}

subs<- list(doubs[[1]][doubs[[1]]$pheno1==-10 & doubs[[1]]$type1%in% c('f', 'b'), ], 
            doubs[[2]][doubs[[2]]$pheno1==-1.5 & doubs[[2]]$type1%in% c('f', 'b'), ])# 'doubs' is the nonpleiotropic data frame 
for (i in 1:2){ 
  names(subs[[i]])[c(1:4)] = c('P1', 'P2', 'P1_which', 'P2_which')
}


# PR 
ggplot(data= dobs[[1]]) + 
  geom_point(aes(x= P2_pheno, y=log2(doubP)), alpha=0.2, shape=4, size=1, col='gray') + 
  facet_grid(.~ P2_which) + theme_classic() + 
  #scale_color_manual(values = c('red', 'blue', 'magenta', 'orange', 'gray')) + 
  labs(x='Mutation B', y='Mutation AB', title= 'A: Folding+Binding') + 
  geom_point(data= subs[[1]][subs[[1]]$P1_which=='f' , ], aes(x= P2, y=phen_log2), shape=3, alpha=0.5, size=0.5, col='red') + 
  geom_point(data= subs[[1]][subs[[1]]$P1_which=='b' , ], aes(x= P2, y=phen_log2), shape=3, alpha=0.5, size=0.5 ,col='blue')+ 
  xlim(-15, 0)
# -10, , 9*3
# PRM 
ggplot(data= dobs[[2]]) + 
  geom_point(aes(x= P2_pheno, y=log2(doubP)), alpha=0.2, shape=4,size=1, col='gray') + 
  facet_grid(.~ P2_which) + theme_classic() + 
  #scale_color_manual(values = c('red', 'blue', 'magenta', 'orange', 'gray')) + 
  labs(x='Mutation B', y='Mutation AB', title= 'A: Folding+Binding') + 
  geom_point(data= subs[[2]][subs[[2]]$P1_which=='f' , ], aes(x= P2, y=phen_log2), shape=3, alpha=0.5, size=0.5, col='red') + 
  geom_point(data= subs[[2]][subs[[2]]$P1_which=='b' , ], aes(x= P2, y=phen_log2), shape=3, alpha=0.5, size=0.5 ,col='blue')+
  xlim(-15, 0)

rm(subs)

######## Fig. 4. how do two different parameter-changes combine? 
##########  
load (doubs_forward, file='ddG_to_output_doubs.RData')

# for the main figure first: Folding and binding 
toplot_type<- c('2,4', '2,3', '2,5', '3,4', '3,5', '4,5') 
xx<-doubs_forward[doubs_forward$paramtyp %in% toplot_type & 
                    doubs_forward$ddG1>=-2.5 & doubs_forward$ddG2>=-2.5 & 
                    doubs_forward$ddG1<=7 & doubs_forward$ddG2<=7, ]


pr= ggplot(data= xx) + 
    geom_tile(aes(x= ddG1, y= ddG2, fill= log2(pr)))+ 
    theme_classic() + labs(x= '', y='', title='PR' ) + scale_fill_viridis( name= 'Mutation AB') + 
    facet_grid(.~paramtyp) +xlim(-2.5, 7) + ylim(-2.5, 7) +  
    geom_contour(aes(x= ddG1, y= ddG2,  z =log2(pr),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) + 
    geom_point(aes(x= 0, y= 0), shape=4, col='white')
pr    

prm= ggplot(data= xx) + 
  geom_tile(aes(x= ddG1, y= ddG2, fill= log2(prm)))+ 
  theme_classic() + labs(x= '', y='', title='PRM' ) + scale_fill_viridis( name= 'Mutation AB') + 
  facet_grid(.~paramtyp) +xlim(-2.5, 7) + ylim(-2.5, 7) +  
  geom_contour(aes(x= ddG1, y= ddG2,  z =log2(prm),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) + 
  geom_point(aes(x= 0, y= 0), shape=4, col='white')

ggarrange(pr, prm, nrow=2, ncol =1)

##### with one example phenotype, what combinations of parameters do they lead to ?  
######## check for PR = -10 or PRM = -2, what combinations of ddGs lead to the phenotypes. 

xx$pr_for_area= abs(log2(xx$pr)- (-10))
xx$prm_for_area= abs(log2(xx$prm)- (-2))
x2<- list(xx[xx$pr_for_area<0.1, ], 
          xx[xx$prm_for_area<0.1, ])

x3<-list()
for (i in 1:2){ 
  x3[[i]]= data.frame(ddG1=NA, ddG2=NA, paramtyp=NA)
  for (j in 1:length(toplot_type)){
    ddG1<- expand.grid(ddG1_1= x2[[i]][x2[[i]]$paramtyp==toplot_type[j], "ddG1"], ddG1_2= x2[[i]][x2[[i]]$paramtyp==toplot_type[j], "ddG1"]) 
    ddG2<- expand.grid(ddG2_1= x2[[i]][x2[[i]]$paramtyp==toplot_type[j], "ddG2"], ddG2_2= x2[[i]][x2[[i]]$paramtyp==toplot_type[j], "ddG2"])  
    
    myx<- data.frame(ddG1= ddG1$ddG1_1+ ddG1$ddG1_2, 
                     ddG2= ddG2$ddG2_1+ ddG2$ddG2_2, 
                     paramtyp = rep(toplot_type[j], length(rownames(ddG1)))) 
    
    x3[[i]]<- rbind(x3[[i]], myx)
    
  } 
  
  x3[[i]]= x3[[i]][!is.na(x3[[i]]$ddG1), ]
}
for (i in 1:2){
  x3[[i]]= x3[[i]][x3[[i]]$ddG1>= -2.5 & x3[[i]]$ddG1<= 7 & 
                     x3[[i]]$ddG2>= -2.5 & x3[[i]]$ddG2<= 7 ,  ]
}



rm(myx); rm(ddG1); rm(ddG2)

xx$pr_yes= 0
xx$prm_yes=0

for (i in 1:length(rownames(x3[[1]]))) {
  xx[ xx$ddG1== x3[[1]][i, 'ddG1'] & xx$ddG2==x3[[1]][i, 'ddG2'] & xx$paramtyp==x3[[1]][i, 'paramtyp'],'pr_yes' ] = 1 
}
for (i in 1:length(rownames(x3[[2]]))) {
  xx[ xx$ddG1== x3[[2]][i, 'ddG1'] & xx$ddG2==x3[[2]][i, 'ddG2'] & xx$paramtyp==x3[[2]][i, 'paramtyp'],'prm_yes' ] = 1 
}

pr= ggplot(data= xx)+geom_tile(aes(x= ddG1, y= ddG2, fill= as.factor(pr_yes)))+ 
  theme_classic() + labs(x= '', y='', title='PR' ) + scale_fill_manual(values=c('white', 'red')) + 
  facet_grid(.~paramtyp) +xlim(-2.5, 7) + ylim(-2.5, 7) +  
  #geom_contour(aes(x= ddG1, y= ddG2,  z =log2(pr),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) + 
  geom_point(aes(x= 0, y= 0), shape=4, col='black')
pr    

prm= ggplot(data= xx) + 
  geom_tile(aes(x= ddG1, y= ddG2, fill= as.factor(prm_yes)))+ 
  theme_classic() + labs(x= '', y='', title='PRM' ) + scale_fill_manual(values=c('white', 'red')) + 
  facet_grid(.~paramtyp) +xlim(-2.5, 7) + ylim(-2.5, 7) +  
  #geom_contour(aes(x= ddG1, y= ddG2,  z =log2(prm),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) + 
  geom_point(aes(x= 0, y= 0), shape=4, col='black')
ggarrange(pr, prm, nrow=2, ncol=1)

###### Fig.5 

####### pleio + pleio -10 & -1 

p= ggplot() + 
  geom_line(data=dose_response, aes(x=log2(PR), y=log2(prm)), size=0.5) + 
  theme_classic() + geom_vline(xintercept = log2(myw_outs[[1]]), lty=2) + 
  geom_vline(xintercept = c(-10, -0.5), col='gray', lty=2) + labs(x='PR', y='PRM') + xlim(-17, 0) + ylim(-17, 0)
  
p
 
p + geom_point(data=ABs2, aes(x= log2(PR), y= log2(PRM)), col='orange', alpha=0.8,  size=2) + xlim(-20, 0) + ylim(-20, 0)



#
##### the package geometry
#### make alpha shape 
library(geometry)
library(alphashape3d)

##
#####  3D 

load('log_even0.000977_0.5_as_prPheno_doub_sampled_90000.RData')
tmp3<-  as.matrix(ABs2[, c(7:9)])


tmp<- as.matrix(pleio3[[1]][, c(1:2,4)])
ashape3d.obj <- ashape3d(tmp, alpha = 1,pert =TRUE ) # alpha increases, with density of the surface
plot(ashape3d.obj, col = 'red', xlim= c(-5, 12), ylim= c(-5, 11), zlim= c(-5, 6),vertices = FALSE)
bg3d("white")
#plot3d(ashape3d.obj, col = 'red', xlim= c(-5, 12), ylim= c(-5, 11), zlim= c(-5, 6))
writeWebGL(filename = file.path("As_2.html"))

tmp2<- as.matrix(pleio3[[2]][, c(1:2,4)])
ashape3d.obj2 <- ashape3d(tmp2, alpha =1) # more sparse, to smooth further
plot(ashape3d.obj2, col = 'cyan', xlim= c(-5, 12), ylim= c(-5, 11), zlim= c(-5, 6),vertices = FALSE)
bg3d("white")

#plot3d(ashape3d.obj2, col = 'cyan', xlim= c(-5, 12), ylim= c(-5, 11), zlim= c(-5, 6))
writeWebGL(filename = file.path("Bs_2.html"))

ashape3d.obj3 <- ashape3d(tmp3, alpha =1)
plot(ashape3d.obj3, col = 'orange', xlim= c(-5, 12), ylim= c(-5, 11), zlim= c(-5, 6),edges = FALSE)
bg3d("white")
# plot3d(ashape3d.obj3, col = 'orange', xlim= c(-5, 12), ylim= c(-5, 11), zlim= c(-5, 6))
writeWebGL(filename = file.path("ABs_2.html"))

