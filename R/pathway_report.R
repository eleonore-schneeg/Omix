#' Generate a report for single omic analysis
#'
#' @param pathways list of pathways 
#' @param report_folder_path report_folder_path folder path to save the report.
#' @param report_file  filename for report (without an extension).
#' @param database Enrichment database `GO_Molecular_Function_2021`,`GO_Cellular_Component_2021`,
#'  `GO_Biological_Process_2021`, `Reactome_2016` , `KEGG_2021_Human` , `MSigDB_Hallmark_2020`
#' @family Report
#' @return Single omic analyses report html
#' @export
#'
pathway_report <- function(pathways,
                              report_folder_path = getwd(),
                              report_file = "single_omic_report_Omix",
                              database='Reactome_2016',
                              num_path=20){
  uniomic=list()
  uniomic$param$pathways=pathways
  uniomic$param$database=database
  uniomic$param$num_path=num_path

  
  report_file <- tools::file_path_sans_ext(report_file)
  
  cli::cli_h2("Generating report for pathways analyses")
  
  
  metadata_tmp_path <- file.path(tempdir(), "metadata.qs")
  
  cli::cli_text("Writing temp files for report...")
  qs::qsave(
    uniomic,
    metadata_tmp_path
  )
  
  krd <- file.path(tempdir(), "krdqc")
  intd <- file.path(tempdir(), "idqc")
  dir.create(krd, showWarnings = FALSE)
  dir.create(intd, showWarnings = FALSE)
  
  cli::cli_text("Generating Pathway report...")
  rmarkdown::render(
    system.file(
      "rmarkdown/templates/pathways/skeleton.Rmd",
      package = "Omix"
    ),
    params = list(
      metadata_path = metadata_tmp_path
    ),
    output_dir = report_folder_path,
    output_file = report_file,
    knit_root_dir = krd,
    intermediates_dir = intd,
    quiet = TRUE
  )
  
  report_file_name <- paste(report_file, ".html", sep = "")
  
  
  
}

