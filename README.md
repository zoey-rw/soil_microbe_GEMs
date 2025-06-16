# Genome-scale Metabolic Models (GEMs) for the soil and rhizosphere microbiome

```mermaid
flowchart TD
    A1[40+ Manually-Curated GEMs:<br/>Heterogeneous publications] --> B
    A2[500+ Template-Based GEMs:<br/>Soil fungi, bacteria, and archaea] --> B[SBML Preprocessing:<br/>Error recovery for problematic files]
    B --> C[Metabolite Standardization:<br/>Pattern detection, extraction & conversion]
    C --> D[Growth Validation:<br/>Verify model equivalence using COBRApy]
    C --> E[Output for Multispecies Frameworks:<br/>Community modeling integration with MICOM and COMETS]
    
    style A1 fill:#e3f2fd,stroke:#1976d2,stroke-width:2px
    style A2 fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
    style B fill:#fff3e0,stroke:#f57c00,stroke-width:2px
    style C fill:#e8f5e8,stroke:#388e3c,stroke-width:2px
    style D fill:#e1f5fe,stroke:#0288d1,stroke-width:2px
    style E fill:#fafafa,stroke:#616161,stroke-width:2px
```

Genome-scale Metabolic Models (GEMs) represent the genomic basis of metabolism for individual microbial species. Soil microbiomes present unique modeling challenges due to extreme diversity, complex interactions, and the predominance of uncultured species. Most soil microbes cannot be grown in pure culture, making template-based modeling approaches essential for community-level simulations.

Published GEMs use different metabolite annotation formats that prevent integration into multi-species modeling frameworks. This repository provides an automated pipeline that standardizes metabolite annotations across both manually-curated and template-based GEMs, enabling constraint-based modeling of soil microbial communities.

## Processing Pipeline

The R-based pipeline handles diverse annotation formats and SBML structural issues:

**Pattern detection:** Automatically identifies annotation formats (RDF, string-based, or multi-database combinations)

**Metabolite standardization:** Converts annotations to MetanetX identifiers using cross-reference databases

**Error recovery:** Handles SBML validation errors, encoding issues, and malformed XML structures

**Growth validation:** Verifies that processed models maintain equivalent growth rates using COBRApy

## Current Status

36 manually-curated GEMs have been processed through the pipeline with the following results:

- 34/36 successful validations (94%)
- 32/34 input files readable by COBRApy (94%)
- 28/34 processed files readable by COBRApy (82%)
- 25/28 simulation equivalent when both files are readable (89%)

The collection includes nitrogen cycle bacteria (ammonia and nitrite oxidizers), rhizobia, soil decomposers, mycorrhizal fungi, yeasts, and methanogens. These curated models are sourced from publications and the BiGG database. 

Over 500 additional models are derived from template-based algorithms such as CarveFungi and COMMIT.

## File Structure

Each species directory for curated models contains:

- Original SBML file(s)
- `*_processed.xml` - Standardized model with MetanetX annotations
- `processing_metadata.json` - Conversion statistics and logs
- `validation_results.json` - COBRApy validation results
- Supplementary information from publication, when available

## Dependencies

**R packages:** sybilSBML, dplyr, stringr, xml2, jsonlite

**Python:** cobra (for validation)

**Reference data:** MetanetX chemical cross-references and deprecated ID mappings (beta release, 2025)

## Community Modeling Integration

Processed models are compatible with:

[![COBRApy](https://img.shields.io/badge/-COBRApy-028?&logo=GitHub)](https://github.com/opencobra/cobrapy) - Flux Balance Analysis (FBA), Flux Variability Analysis (FVA)

[![COMETSpy](https://img.shields.io/badge/-COMETSpy-028?&logo=GitHub)](https://github.com/segrelab/cometspy) - Dynamic FBA, spatiotemporal simulations

[![MICOM](https://img.shields.io/badge/-MICOM-028?&logo=GitHub)](https://github.com/micom-dev) - Community trade-offs, growth rate estimation
