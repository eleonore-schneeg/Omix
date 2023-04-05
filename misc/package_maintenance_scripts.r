#install.packages(c("devtools", "roxygen2", "usethis", "spelling"))
library(devtools)
library(roxygen2)
library(usethis)
library(spelling)

## Ongoing package dev
use_r("batch_correction_protein")
use_r("circular_heatmap")
use_r("clusters_DE")
use_r("combat_correction")
use_r("custom_rna")
use_r("DE_formating")
use_r("DE_multiomics")
use_r("enrichment")
use_r("extract_multiomic_signature")
use_r("extract_weights")
use_r("fgsea")
use_r("filter_protein")
use_r("generate_multiassay")
use_r("get")
use_r("imputation_protein")
use_r("integration_visualisations")
use_r("integrative_models")
use_r("integrative_results_clustering")
use_r("integrative_results_supervised")
use_r("integrative_results_unsupervised")
use_r("median_centering_correction")
use_r("multiomics_clustering_heatmap")
use_r("multiomics_modules")
use_r("open_targets")
use_r("process_protein")
use_r("process_rna")
use_r("protein_DE_analysis")
use_r("pseudotime_inference")
use_r("pseudotime_ordering_modules")
use_r("rna_DE_analysis")
use_r("single_omic_comparison")
use_r("vertical_integration")

# add to DESCRIPTION
# Only add the pkg names which were directly used in any of the functions
# Additional dependencies can go to requirements-bioc.R for dockerfile generation.

usethis::use_github_action("check-standard")

usethis::use_pkgdown()
pkgdown::init_site() #create favicons from pkg logo
pkgdown::build_site()


### BEFORE EVERY COMMIT
#Restart R Session Cmd+Shift+F10 (Ctrl+Shift+F10 for Windows)
#Document Package Cmd+Shift+D (Ctrl+Shift+D for Windows)
#Check Package Cmd+Shift+E (Ctrl+Shift+E for Windows)

## before every release
# knit the readme.Rmd <<--
# update the site
usethis::use_version()
pkgdown::build_site()


options(rgl.useNULL = TRUE)

#########################
