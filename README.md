# FusioMR-analysis

This repository contains simulation and data analysis R code for the FusioMR manuscript. 

Folder `Functions`: data generation and utility functions used in simulation and real data analysis.

You need to change the directory for inputs or outputs to your own.

Codes are largely self-explained. If you have any quesions, please feel free to leave a message here or contact kbw@uchciago.edu. 

## Files description
* `simu_code_full.R`: Main file for the first set of simulation studies with autoregressive LD structure, including scenarios 1 and 2, and the weak-IV scenario.
* `weakly_invalid_IV.R`: The additional weakly invalid IV scenario.
* `Iy_weakly_associated_Ix.R`: The additional scenario where SNPs in $\mathcal{I}_X$ and SNPs in $\mathcal{I}_Y$ are uncorrelated or only weakly correlated.
* `simu2_s1.R`: The second set of simulation studies, scenario 1.
* `simu2_s2.R`: The second set of simulation studies, scenario 2.
* `process_ARIC.R`: Code to preprocess real data in the proteome-wide application.
* `Analysis_redo.R`: Code to perform cisMR-cML on the processed real data.
* `Analysis_GIVW.R`: Code to perform LDA-Egger, GIVW, and GEgger on the processed real data.
* `coloc.R`: Code to run colocalization in real data analysis.
