# soil_microbe_GEMs

Compilation of published genome-scale metabolic models (GEMs) for constraint-based modeling of the soil microbiome.
Repository Structure
Each organism has its own directory containing the original SBML files and processed versions. Models are sourced from various publications in different formats (primarily XML/SBML).
Processing Pipeline
The repository includes an automated R pipeline that standardizes metabolite annotations across all models:

Pattern Detection: Automatically identifies annotation formats (BiGG, RDF, string-based, or multi-database combinations)
Metabolite standardization: Converts annotations to MetanetX identifiers while preserving simulation outputs
Error recovery: Handles problematic SBML structures including Matrix errors and malformed group elements
Validation: Verifies growth rates match original publications using COBRApy

### Supported Annotation Patterns
The pipeline handles various annotation patterns and structures, including BiGG, ChEBI, KEGG, MetanetX, etc.

### File Structure
Each species directory contains:

*_input.xml - Original SBML file(s)
*_processed.xml - Standardized model with MetanetX annotations
processing_metadata.json - Processing log and conversion statistics
validation_results.json - COBRApy validation results
When available, files associated with original publication

### Community model integration
Processed models are compatible with standard constraint-based modeling frameworks:

    [![COBRApy](https://img.shields.io/badge/-COBRApy-028?&logo=GitHub)](https://github.com/opencobra/cobrapy)&nbsp; standard Flux-Balance Analysis (FBA), Flux-Variability Analysis (FVA) 
<br>
    [![COMETSpy](https://img.shields.io/badge/-COMETSpy-028?&logo=GitHub)](https://github.com/segrelab/cometspy)&nbsp; dynamic FBA, spatiotemporal simulations
<br>
    [![MICOM](https://img.shields.io/badge/-MICOM-028?&logo=GitHub)](https://github.com/micom-dev)&nbsp; community trade-offs, growth rate estimation


Processing Status
60+ models are currently being processed through the standardization pipeline. Models undergo validation to ensure growth rates match published values after annotation standardization.
Dependencies

R packages: sybilSBML, tidyverse, stringr, xml2, jsonlite
Python: cobra (for validation)
Reference data: MetanetX chemical cross-references and deprecated ID mappings (beta release, 2025)