# soil_microbe_GEMs
Compilation of published genome-scale metabolic models (GEMs) for constraint-based modeling of the soil microbiome

Each organism has its own directory. GEMs are uploaded in various formats depending on the original publication (usually XML/SBML). An R script was used to enforce a common nomenclature for all exchange metabolites (`01_modify_metabolites.r`). Python was used to verify COBRA compliance and ensure that growth rates match the rates described in publications. (`/02_make_cobra_compliant.ipynb`). GEMs have also been evaluated for growth on various carbon substrates (`evaluate_carbon_usage.ipynb`). 

Simulation-ready GEMs can be integrated into community modeling frameworks of various flavors:<br>
    [![COBRApy](https://img.shields.io/badge/-COBRApy-028?&logo=GitHub)](https://github.com/opencobra/cobrapy)&nbsp; standard Flux-Balance Analysis (FBA), Flux-Variability Analysis (FVA) 
<br>
    [![COMETSpy](https://img.shields.io/badge/-COMETSpy-028?&logo=GitHub)](https://github.com/segrelab/cometspy)&nbsp; dynamic FBA, spatiotemporal simulations
<br>
    [![MICOM](https://img.shields.io/badge/-MICOM-028?&logo=GitHub)](https://github.com/micom-dev)&nbsp; community trade-offs, growth rate estimation

60+ more GEMs to come soon!
