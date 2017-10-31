# Armillaria_TE_annotations 

Supplementary methods and additional information for transposable element analyses in  [https://www.nature.com/articles/s41559-017-0347-8](https://www.nature.com/articles/s41559-017-0347-8)

Sipos, György, Arun N Prasanna, Mathias C Walter, Eoin O’Connor, Balázs Bálint, Krisztina Krizsán, Brigitta Kiss, et al. 2017. “Genome Expansion and Lineage-Specific Genetic Innovations in the Forest Pathogenic Fungi Armillaria.” Nature Ecology & Evolution 12 (October): 515. doi:10.1038/s41559-017-0347-8.

## REPET configuration files and scripts

Sample configuration files and scripts for _de novo_ discovery and annotation of transposable elements using [REPET](https://urgi.versailles.inra.fr/Tools/REPET). 

Provided is a representative example for a single species (_Agaricus bisporus_). All other species were processed using the same parameters & pipeline. 

### REPET _de novo_

-  	TEdenovo.cfg:			Configuration file for the TE _de novo_ pipeline
- 	TE_denovo_pipeline.sh:	_de novo_ pipeline script
-	TEdenovo_panSpecies.cfg:	Configuration file for clustering combined consensus sequences (step 8 of the pipeline)

### REPET anno

-	TEannot.cfg:				Configuration file for the TE anno pipeline
- 	TE_anno_pipeline.sh:		Anno pipeline script

## TE annotations

 *chr_allTEs_nr_noSSR_join_path.map.merge:	TE annotations for each genome
 
 Format: ElementID (see REPET website for naming scheme)	scaffold	start end
 
## Windowed analysis

-  run_windowed_analysis.sh:	Script to generate data for windowed analysis
-  Windowed_analysis.R:		Plot windowed analysis









