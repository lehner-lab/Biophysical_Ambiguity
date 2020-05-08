# Biochemical_Ambiguity
The R scripts here are used to generate the modelling data and figures for the manuscript: "Biochemical ambiguities prevent accurate genetic prediction".

# Brief description of the study: 
It is widely assumed that the combined effects of genetic variants can be determined from their individual effects on a trait. Here, we show that this assumption is incorrect and that even with perfect measurements and a mechanistic understanding of a system it is often impossible to predict what happens when two mutations are combined without additional information.  This apparent paradox arises because mutations can have many different biochemical effects to cause the same change in a phenotype. When combining mutations, the outcome can be very different depending upon what these hidden biochemical changes actually are.  Using Lambda repressor (CI) - Operator system, we show that accurate genetic prediction of phenotypes and disease will sometimes not be possible unless these biochemical ambiguities can be resolved by making additional measurements.

# 1. System requirements:
* **1.1** R >= v 3.3.3 (with packages: rootSolve, optimize, ggplot2, ggpubr, viridis, reshape2, stringr)
* **1.2** Versions the software has been tested on: R (v 3.3.3) and RStudio (v 1.1.463) is used for this study, to run R scripts.

# 2. Installation guide: 
No installation necessary, functions loaded on-the-fly with "source('FUNCTION')" (see below). 

# 3. Demo: 
The two R scripts listed below ("Model_biochemical_ambiguity.R" and "Plotting_for_Model_biochemical_ambiguity.R") reproduce the simulated datasets and the plots in the manuscript to draw the conclusion. The run-time is not more than a couple of hours in total on a normal desktop computer.

For simple run-test of the custom functions, three demo files are included (see below Demo data list). With the data frame, the functions are applied by rows. Some cases NA will be returned, if the given phenotype is not possible with indicated mutation type. 

# 4. Instructions for use
The listed scripts (#1-6) are specific to the lambda repressor – operator system (OR1, OR2, and OR3), based on the experimental and theoretical knowledge on the gene-regulatory system. 
The script #7 ("Generalizability.R") can be used for any protein-protein interaction pair where measured phenotype is linear to the protein-protein complex concentration. In that case, parameters for the ∆∆G of binding, and ∆∆G of folding energy of both proteins need to be updated to suit the users’ interest.

## Details of each script (including functions for generating datasets used in the study): 

* **1. Forward_function_param_to_Output.R** This function uses parameter values to calculate downstream phenotypes (expressions from PR and PRM promoters)
* **2. Reverse_function_pheno_to_biochem.R** This function uses phenotypic values to calculate mutational effects on indicated biochemical changes.
* **3. Reverse_function_for_pleiotropic_muts.R** This function uses phenotypic values, protein-folding and dimerization parameters as input to calculate mutational effects on DNA-binding parameter.
* **4. dose_response.R** This function generates data to plot dose-response curves for PR and PRM expressions as a function of wild type CI levels.
* **5. Model_biochemical_ambiguity.R** This script generates data for plotting how mutations affecting different biochemical parameters combine, as an example.
* **6. Plotting_for_Model_biochemical_ambiguity.R** Script used for plotting are listed here, as an example.
* **7. Generalizability.R** Script to generate datasets and plot to examine whether biochemical ambiguities generate phenotypic unpredictability for a protein-protein interaction pair. 

## List of demo datasets. 

* **1. Demo1a.RData**  (to run ‘pr_prm_Intramol_protein_withDNA_whichparam’  in the Forward_function_param_to_Output.R. This allows specifying which two parameters are changed by which amount for a given mutant.) 
* **2. Demo1b.RData**  (to run ‘Forward’ function in the Forward_function_param_to_Output.R. This inputs protein total amount expressed, and all the parameter changes for each mutant. ) 
* **3. Demo2.RData**  (to run Reverse_function_pheno_to_biochem.R To be noted: Some phenotypes are not reachable with mutations affecting tetramerization alone, or folding-alone. In that case, the function will return NA ) 
* **4. Demo3.RData**  (to run Reverse_function_for_pleiotropic_muts.R)
* **5. Demo4.RData**  (to run dose_response.R) 
* **6. Demo5.RData**  (to run ‘myppi’ function in the file  Generalizability.R)

