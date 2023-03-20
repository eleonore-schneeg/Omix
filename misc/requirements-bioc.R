bioc_pkgs<-c(
  "MultiAssayExperiment",
  "SummarizedExperiment",
  "proBatch",
  "mixOmics",
  'basilisk',
  "MOFA2",
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
  'EWCE',
  'GenomicRanges',
  'AnnotationDbi',
  'WGCNA',
  'limma',
  'S4Vectors',
  'slingshot',
  'DESeq2',
  'limma',
  'Biobase',
  'EnsDb.Hsapiens.v86',
  'msigdbr'
)


requireNamespace("BiocManager")
BiocManager::install(bioc_pkgs, ask=F)
