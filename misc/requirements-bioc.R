bioc_pkgs<-c(
  "MatrixGenerics",
  "DelayedArray",
  "ensembldb",
  "proBatch",
  "mixOmics",
  'basilisk',
  "ComplexHeatmap",
  "circlize",
  'BiocGenerics',
  'BiocStyle',
  'fgsea',
  'org.Hs.eg.db',
  'clusterProfiler',
  'enrichplot',
  'biomaRt',
  'edgeR',
  'GenomicRanges',
  'AnnotationDbi',
  'WGCNA',
  'limma',
  'S4Vectors',
  'DESeq2',
  'limma',
  'Biobase',
  'EnsDb.Hsapiens.v86',
  'msigdbr',
  "SummarizedExperiment",
  'SingleCellExperiment',
  "MultiAssayExperiment"
)


requireNamespace("BiocManager")
BiocManager::install(bioc_pkgs, ask=F)
