bioc_pkgs<-c(
  'basilisk',
  'Biobase',
  'BiocGenerics',
  'BiocStyle',
  'circlize',
  'clusterProfiler',
  'ComplexHeatmap',
  'DelayedArray',
  'DelayedMatrixStats',
  'MatrixGenerics',
  'edgeR',
  'enrichplot',
  'SummarizedExperiment',
  'fgsea',
  'GenomicRanges',
  'graph',
  'HDF5Array',
  'IRanges',
  'limma',
  'mixOmics',
  'msigdbr',
  'MultiAssayExperiment',
  'multtest',
  'AnnotationDbi',
  'ensembldb',
  'EnsDb.Hsapiens.v86',
  'org.Hs.eg.db',
  'preprocessCore',
  'proBatch',
  'rhdf5',
  'S4Vectors',
  'SingleCellExperiment',
  'EWCE',
  'slingshot',
  'WGCNA',
  'DESeq2',
  'MOFA2'
)


requireNamespace("BiocManager")
BiocManager::install(bioc_pkgs, ask=F)

