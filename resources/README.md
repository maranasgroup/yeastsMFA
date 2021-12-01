#### File descriptions
**1) carbon_mappings_and_annotations.xlsx (version: 2021-05-01)**<br>
Contain mappings of carbon atoms between metabolites in a reactions. In addition, annotations of database entries for reactions and the information source of carbon mappings are also provided (this is specifically for yeasts, not for other organisms). Annotations of functional group and adjacent carbon atoms are provided for carbon atoms to facilitate future usage and re-examiniation.

**2) dilutions_network_modules.xlsx**<br>
Dilution network modules correspond to metabolites in labeling dataset need to be added to the 13C-MFA run input files before running. These dilution network modules allow the program to map simulated labeling data to experimental labeling data. This file is specific to yeasts (which account for cytosol and mitochondria compartments for some metabolites).

**3) metabolite_ID_mapping.txt**<br>
Mapping of metabolite ID in the model (first column) to the metabolite ID in the raw dataset (second column). Names, formulae, and mapping to KEGG database entries are also provided.

**4) tracer_info.txt**<br>
Run settings regarding the tracers. This file depends on how experiments are set up.

**5) moiety_mapping.py**<br>
Unique number index assignments for some cofactors to help with mapping readability since these cofactors can contain a large number of carbons.

