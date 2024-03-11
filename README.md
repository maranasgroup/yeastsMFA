# File repository for 13C-metabolic flux analysis of yeasts
This repository provides Supplementary files regarding 13C-metabolic flux analysis (13C-MFA) that are part of the manuscripts:<br>
**Mitochondrial ATP generation is more proteome efficient than glycolysis**<br>
Yihui Shen, Hoang V. Dinh, Edward R. Cruz, Zihong Chen, Caroline R. Bartman, Tianxia Xiao, Catherine M. Call, Rolf-Peter Ryseck, Jimmy Pratas, Daniel Weilandt, Heide Baron, Arjuna Subramanian, Zia Fatma, Zong-Yen Wu, Sudharsan Dwaraknath, John I. Hendry, Vinh G. Tran, Lifeng Yang, Yasuo Yoshikuni, Huimin Zhao, Costas D. Maranas, Martin Wühr & Joshua D. Rabinowitz<br>
Nature Chemical Biology. 2024. [doi: https://doi.org/10.1101/2022.08.10.503479](https://doi.org/10.1038/s41589-024-01571-y)
<br>
<br>
**Comparative study of two *Saccharomyces cerevisiae* strains with kinetic model at genome-scale**<br>
Mengqi Hu, Hoang V. Dinh, Yihui Shen, Patrick F. Suthers, Charles Foster, Catherine M. Call, Xuanjia Ye, Zia Fatma, Huimin Zhao, Joshua D. Rabinowitz, and Costas D. Maranas<br>
Metabolic Engineering. 76:1-17. 2023. https://doi.org/10.1016/j.ymben.2023.01.001<br>
<br>
**References notice**: Please cite the paper corresponding to the dataset that you use. The details are provided in the README.md file within the sub-directories. If you use the common resources and scripts, please cite Shen et al., 2022.<br>

#### File descriptions
**1) mfa_scripts**<br>
Source MATLAB scripts to run 13C-MFA. Available at https://github.com/maranasgroup/SteadyState-MFA.

**2) run_scripts**<br>
Scripts that combine and arrange source MATLAB 13C-MFA scripts.

**3) resources**<br>
Collected resources regarding carbon mappings and other setups to run 13C-MFA scripts

**4) mfa**<br>
13C-MFA run files and results, arranged by directory for each organism.<br>
* ***./mfa/\<dataset\>/metabolic_rxns_mappings***: Metabolic network and carbon mappings of metabolic reactions only. The input files are generated from this metabolic_rxns_mappings and ./resources/dilutions_network_modules.xlsx
* ***./mfa/\<dataset\>/run_files***: list of input files for 13C-MFA software run
* ***./mfa/\<dataset\>/result_files***: processed output files containing simulated 13C-labeling patterns and best-fit metabolic fluxes
    
There are the following \<dataset\> in this repository:
* S_cerevisiae: 14 datasets on S. cerevisiae (Shen et al., 2022)
* I_orientalis: 16 datasets on I. orientalis (Shen et al., 2022)
* Hu2022_KFIT_SC: 9 datasets on S. cerevisiae (Hu et al., 2023)
