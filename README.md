# File repository for 13C-metabolic flux analysis of yeasts
This repository provides Supplementary files regarding 13C-metabolic flux analysis (13C-MFA) that is part of the manuscript:<br>
**Proteome capacity constraints favor respiratory ATP generation**<br>
Yihui Shen, Hoang V. Dinh, Edward Cruz, Catherine M. Call, Heide Baron, Rolf-Peter Ryseck, Jimmy Pratas, Arjuna Subramanian, Zia Fatma, Daniel Weilandt, Sudharsan Dwaraknath, Tianxia Xiao, John I. Hendry, Vinh Tran, Lifeng Yang, Yasuo Yoshikuni, Huimin Zhao, Costas D. Maranas, Martin Wuhr, Joshua D. Rabinowitz<br>
bioRxiv 2022.08.10.503479; doi: https://doi.org/10.1101/2022.08.10.503479
<br>
#### File descriptions
**1) mfa_scripts**<br>
Source MATLAB scripts to run 13C-MFA. Available at https://github.com/maranasgroup/SteadyState-MFA.

**2) run_scripts**<br>
Scripts that combine and arrange source MATLAB 13C-MFA scripts.

**3) resources**<br>
Collected resources regarding carbon mappings and other setups to run 13C-MFA scripts

**4) mfa**<br>
13C-MFA run files and results, arranged by directory for each organism.<br>
* ***./mfa/organism/metabolic_rxns_mappings***: Metabolic network and carbon mappings of metabolic reactions only. The input files are generated from this metabolic_rxns_mappings and ./resources/dilutions_network_modules.xlsx
* ***./mfa/organism/run_files***: list of input files for 13C-MFA software run
* ***./mfa/organism/result_files***: processed output files containing simulated 13C-labeling patterns and best-fit metabolic fluxes
    
There are the following organisms in this repository:
* *Saccharomyces cerevisiae*: 14 datasets
* *Issatchenkia orientalis*: 16 datasets
