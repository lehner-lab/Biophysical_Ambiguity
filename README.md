# Biochemical_Ambiguity
# The R scripts here are used to generate the modelling data and figures for the manuscript: 
"Biochemical ambiguities prevent accurate genetic prediction"  
RData generated from the scripts are not included.  
Here is the description for each function: 

1) Forward_function_param_to_Output.R 
This function uses parameter values to calculate downstream phenotypes (expressions from PR and PRM promoters)

2) Reverse_function_pheno_to_biochem.R
This function uses phenotypic values to calculate mutational effects on indicated biochemical changes.  

3) Reverse_function_for_pleiotropic_muts.R
This function uses phenotypic values, protein-folding and dimerization parameters as input to calculate mutational effects on DNA-binding parameter.   

4) dose_response_parameter_response_curves.R
This code generates data to plot dose-response curves for PR and PRM expressions as a function of wild type CI levels. 
(Fig. 1A)

5) dose_response_for_mutants_example.R 
This code generates dose-response curves for mutations with the same phenotypes (but with different biochemical effects) at an expression level of CI. 
(Fig. S6)

6) Model_biochemical_ambiguity.R 
This code generates data for plotting how mutations affecting different biochemical parameters combine. 

7) Plotting_for_Model_biochemical_ambiguity.R
All the codes for plotting are listed here, except the dose-response curve codes that are included in the dose-response R curves (No.4 and No.5 aobve)

# To be noted: 
Some functions may give errors but it can be safely ignored. 
