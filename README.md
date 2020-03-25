# Biochemical_Ambiguity
# The R scripts here are used to generate the modelling data and figures for the manuscript: 
"Biochemical ambiguities prevent accurate genetic prediction"  
RData generated from the scripts are not included.  

# Brief description of the study: 
It is widely assumed that the combined effects of genetic variants can be determined from their individual effects on a trait. Here, we show that this assumption is incorrect and that even with perfect measurements and a mechanistic understanding of a system it is often impossible to predict what happens when two mutations are combined without additional information.  This apparent paradox arises because mutations can have many different biochemical effects to cause the same change in a phenotype. When combining mutations, the outcome can be very different depending upon what these hidden biochemical changes actually are.  Using Lambda repressor (CI) - Operator system, we show that accurate genetic prediction of phenotypes and disease will sometimes not be possible unless these biochemical ambiguities can be resolved by making additional measurements.


## Here is the description for each function: 

1) Forward_function_param_to_Output.R 
This function uses parameter values to calculate downstream phenotypes (expressions from PR and PRM promoters)

2) Reverse_function_pheno_to_biochem.R
This function uses phenotypic values to calculate mutational effects on indicated biochemical changes.  

3) Reverse_function_for_pleiotropic_muts.R
This function uses phenotypic values, protein-folding and dimerization parameters as input to calculate mutational effects on DNA-binding parameter.   

4) dose_response.R
This code generates data to plot dose-response curves for PR and PRM expressions as a function of wild type CI levels. 

6) Model_biochemical_ambiguity.R 
This code generates data for plotting how mutations affecting different biochemical parameters combine. 

7) Plotting_for_Model_biochemical_ambiguity.R
All the codes for plotting are listed here.

