if(pipe=="Seurat"){
  # a = c('Seurat','survival','leiden','spatstat.sparse','BiocNeighbors','scater',
  #       'uwot','reticulate','scuttle','BiocSingular','mgcv','beachmat',
  #       'DelayedMatrixStats','rsvd','spatstat.explore','sctransform',
  #       'sparseMatrixStats','SummarizedExperiment','spatstat.data',
  #       'ScaledMatrix','DelayedArray','irlba','SeuratObject')
  pkgs =  c( 'scater', 'SingleCellExperiment', 'uwot', 'rtracklayer', # 'scran', 
             'GenomicFeatures', 'GenomicRanges', 'IRanges', 'clusterProfiler', 'rliger','SeuratWrappers',#'install',
             'ComplexHeatmap', 'ggsci', 'RColorBrewer', 'ggplotify', 'ggpubr', 'scales', 'cowplot', 'patchwork', 'tidyr', 'clustree',
             'VennDiagram', 'ggvenn', #'ggVennDiagram', 
             'tidyverse', 'plyr','data.table', 'dplyr', #'Seurat',
             'umap', 'reshape2' )
}else if( pipe=='scp' ){
  pkgs =  c( 'jsonlite', 'languageserver','IRanges','GenomicRanges','ComplexHeatmap',
             'pheatmap', 'data.table', 'tidyverse', 'magrittr', 'RColorBrewer', 'cowplot', 'patchwork', 'ggridges',
             'scales', 'textshape', 'stringr',"scuttle","SingleCellExperiment",
             'impute', 'scater', 'sva', 'GGally','limma',
             'scpdata', 'QFeatures', 'scp', 'Seurat')
}else if( pipe=='rna-seq'){
  pkgs = c( ###plotting
    'tidyverse','cowplot','pheatmap','ComplexHeatmap','RColorBrewer','ggsci',
    'ggsignif','ggpubr','ggplotify','patchwork','scales','ggrepel','reshape2','textshape',
    ### manipulate
    'dplyr','data.table','parallel','preprocessCore','readxl','WriteXLS',
    ## diff gene analysis
    'IsoformSwitchAnalyzeR','Rsubread','DESeq2','sva',
    ## Genomic ranges and genome
    'org.Mm.eg.db','ensembldb','FactoMineR','factoextra',
    ## GO related
    'topGO','clusterProfiler','enrichplot','pathview','Rgraphviz','grDevices')
  
  ### analysis functions ###
  Deseq2_dds = function(input_counts,index,designGroup = character()) {
    require(DESeq2)
    dds = DESeqDataSetFromMatrix(countData = input_counts,
                                 colData = index,
                                 design = ~ PhenoGroup )
    keep = rowSums(counts(dds))>10
    dds = dds[keep,]
    dds= DESeq( dds )
    return(dds)
  }
  cal_TPM = function(data,length_var="Length"){
    length = data[,which(colnames(data)==length_var)]
    kb = length/1000
    rpk = data[,-which(colnames(data)==length_var)]/kb
    tpm = apply(rpk,2,function(x){x/(sum(x)/1000000)})
    return( tpm )
  }
  
}



Loadpkgs = function(i){
  print(paste("Installing package:", i))
  if(!require(i, quietly = TRUE)){
    BiocManager::install(i, update = FALSE)
  }
  library(i, character.only = TRUE)
}

sapply( pkgs,Loadpkgs )

