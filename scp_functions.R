cellnames = NULL
sc.run = function(data){
  sc.runs = intersect(unique(data$Raw.file[data$celltype %in% cellnames]),names(data))
  # sc.runs = unique(data$Raw.file[data$celltype %in% c("LPS","untreated")])
  return(sc.runs)
}
scp.stats = function(data){
  sc.runs = intersect(unique(data$Raw.file[data$celltype %in% cellnames]),names(data))

  tmp = rbindRowData(data,i = sc.runs)
  
  all.runs = length(na.omit(intersect(unique(data$Raw.file),names(data))))
  n.peptide = length(na.omit(unique(tmp$Modified.sequence))) #tmp$Modified.sequence
  n.protein = length(na.omit(unique(tmp$Leading.razor.protein)))#tmp$Proteins
  n.cell = length( data$celltype[(data$celltype %in% cellnames)& data$Raw.file %in% sc.runs ] )
  sc.runs = length( sc.runs )
  print(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
               ", peptides = ",n.peptide,", protein = ", n.protein))
  return(paste0("all.runs = ",all.runs,", sc.runs = ",sc.runs,", cells = ",n.cell,
                ", peptides = ",n.peptide,", protein = ", n.protein))
}
keep.run = function(data){
  sc.runs = unique(data$Raw.file[data$SampleType %in% cellnames])
  # sc.runs = unique(data$Raw.file[data$celltype %in% c("Blank")])
  return(sc.runs)
}


# pep_mm = ev_update$Sequence #%>% grep(.,pattern = "K",invert = T,value = T)
# # pep = rownames(assay(scp[["peptides"]]))
# pep_hu = fread("/lustre/user/liclab/fangq/proj/scp/scope_protocol_data/Data_for_Fig_3-4_MaxQuant_txt_files/ev_updated.txt",header = T)
# pep_hu = pep_hu$Sequence
# 
# aa = toupper(c('a','r','n','d','c','q','e','g','h','i','l','k','m','f','p','s','t','w','y','v')) %>% as.list
# names(aa) = toupper(c('a','r','n','d','c','q','e','g','h','i','l','k','m','f','p','s','t','w','y','v')) 
# 
# aa = lapply( aa,function(x){
#   exist_1 = length(grep(pep_mm,pattern = x))/length(pep_mm)
#   exist_2 = length(grep(pep_hu,pattern = x))/length(pep_hu)
#   return(c(exist_1,exist_2))
# } )
# 
# library(plyr)
# 
# aa=ldply(aa,.id = "AA")
# 
# colnames(aa)[2:3] = c("MM_public","HU_public")
# aa$MMrank = rank(as.numeric(-aa$MM_public))
# aa$HUrank = rank(as.numeric(-aa$HU_public))
# length(unique(c(grep(pep,pattern = "K",value = F),grep(pep,pattern = "R",value= F))))/length(pep)
# 
# length(intersect(grep(pep,pattern = "K",value = F),grep(pep,pattern = "E",value= F)))/length(grep(pep,pattern = "K",value = F))
# 

Step1_generate_raw = function(ev_data = data.frame(),
                              meta_data = data.frame(),
                              batchCol = "Raw.file",channelCol = "tmt.label",
                              filt_with_min_psm = logical(),
                              psm_threshold = numeric(),
                              samplePattern = "Mono|Macro|CD4|CLP|CMP|MPP|MEP|Neutro",
                              uniq_protein_column = "Leading.razor.protein",
                              PEP_colname = 'dart_PEP',
                              carrierPattern_scr = 'Carrier',
                              sc_batch = sc.run(scp)
                              ){
    require(future)
    future::plan("multisession",workers = 10)
    scp = readSCP(featureData = ev_data,
                  colData = meta_data,
                  batchCol = batchCol,
                  channelCol = channelCol,
                  removeEmptyCols = TRUE)
    scp = zeroIsNA(scp,i = names(scp))
    n.psm = dims(scp)[1, ] ;summary(n.psm)       # psm counts
    n.runs = dims(scp)[2, ] ;summary(n.runs)      # run tmt counts
    if(filt_with_min_psm){
      keep = dims(scp)[1,][dims(scp)[1,]>psm_threshold] %>% names()
      scp = scp[,,keep]
    }
    scp = pep2qvalue(scp,
                      i = names(scp), ## need to be presnt with single cell runs
                      PEP = PEP_colname,
                      rowDataName = "qvalue_PSMs")
    scp = pep2qvalue(scp,
                      i = names(scp), ## need to be presnt with single cell runs
                      PEP = PEP_colname,
                      groupBy = uniq_protein_column,
                      rowDataName = "qvalue_proteins")
    scp = computeSCR(scp,
                       i = sc_batch,
                       colvar = "celltype",
                        carrierPattern = carrierPattern_scr,
                        samplePattern = samplePattern,
                        sampleFUN = "mean",
                        rowDataName = "MeanSCR")
    save(scp,file = "analysis/Step1_raw_scp.RData",compress = T)
    scp.stats(scp)
    print(dims(scp))
    assign("scp",scp,envir = .GlobalEnv)
    future::plan("sequential")
    }

Step2_psm_to_pep = function(scp = QFeatures(),
                            Raw.file_feature = "^\\d+",
                            medianCV_batches = character(),
                            medianCV_nobs = 5,
                            if_generate_psm = logical(),
                            uniq_protein_column = "Leading.razor.protein",
                            sc_batch = sc.run(scp)){
    require(future)
    future::plan("multisession",workers = 20)
    if(if_generate_psm){
      cat("Raw PSM aggregate Starting...",fill = T)       
      scp = joinAssays( scp,
                      i = grep(Raw.file_feature,names(scp),value = T),
                      name = "psms_raw" )
      cat("Raw PSM aggregate done! Now Starting normalization via Reference & aggregation of normlized PSM...",fill = T)
    }else(cat("Skip PSM aggregating ! Now Starting normalization via Reference & aggregation of normlized PSM...",fill = T))
    
    scp = divideByReference( scp,
                            i = sc_batch, ### must have Reference or will be error
                            colvar = "celltype",
                            samplePattern = ".",
                            refPattern = "Reference" )
    cat("Normalization Done! Now calculate medianRI and medianCV for each pep & cell...",fill = T)
    scp = medianCVperCell( scp, 
                           i = medianCV_batches,
                           groupBy = uniq_protein_column,
                           nobs = medianCV_nobs,  # least peptide numbers per protein,2 and 5
                           norm = "SCoPE2",
                           na.rm = T, 
                           colDataName = "MedianCV" )
    if(if_generate_psm){
      scp = joinAssays( scp,
                        i = grep(Raw.file_feature,names(scp),value = T),
                        name = "psms" )
    }
    cat("PSM aggregate to Peptides Starting......",fill = T)
    scp = aggregateFeaturesOverAssays(scp,
                                      i = grep(Raw.file_feature,names(scp),value = T),
                                      fcol = "Modified.sequence",
                                      name = paste0("peptides_",grep(Raw.file_feature,names(scp),value = T)),
                                      fun = matrixStats::colMedians, na.rm = TRUE)
    cat("Now joining into peptides assay, may take long time if dataset is big",fill = T)
    peptide.runs = grep("^peptides_",names(scp),value = T)
    scp = joinAssays( scp,
                      i = peptide.runs,
                      name = "peptides")  ### may take long time
    plot(scp, interactive = T)
    cat("Saving scp in pep_aggregated_scp.RData...",fill = T)
    save(scp,file = "analysis/Step2_pep_aggregated_scp.RData",compress = T)
    scp.stats(scp)
    cat(dims(scp))
    plot(scp)
    assign("scp",scp,envir = .GlobalEnv)
    future::plan("sequential")
        }

Step3_normalization_pep_aggre_protein = function(scp = QFeatures(),
                                                 filt_pep_NA = F,pep_pNA = 0.99,
                                                 filt_pro_NA = T,protein_pNA = 0.99,
                                                 uniq_protein_column = "Leading.razor.protein"){
    require(future)
    future::plan("multisession",workers = 20)

    cat("Now starting peptide normalization and log tranformation...",fill = T)
    scp = sweep(scp, 
              i = "peptides", ## assay to be normalized
              MARGIN = 2,
              FUN = "/",
              STATS = colMedians(assay(scp[["peptides"]]), na.rm = TRUE), ## *****Can be changed by different methods for normalization (eg, quantile?)
              name = "peptides_norm_col")

    scp = sweep(scp,
                i = "peptides_norm_col",  ## after col normalization do row normalization 
                MARGIN = 1,
                FUN = "/",
                STATS = rowMeans(assay(scp[["peptides_norm_col"]]),  na.rm = TRUE),
                name = "peptides_norm")
    if(filt_pep_NA){
      cat(paste0("Filtering out sparse peptide that over ",pep_pNA*100,"% is NA....."),fill= T)
      scp = filterNA(scp,
                     i = "peptides_norm",
                     pNA = pep_pNA)
    }else{
      cat("Skip filteringout sparse peptides",fill= T)
    }


    scp = logTransform(scp,
                        base = 2,
                        i = "peptides_norm",
                        name = "peptides_log")
    cat("Now aggregating peptide into protein, normalization and log tranformation...",fill = T)

    scp = aggregateFeatures(scp,
                            i = "peptides_log",
                            name = "proteins",
                            # fcol = "Leading.razor.protein",
                            fcol = uniq_protein_column,
                            fun = matrixStats::colMedians, 
                            na.rm = TRUE)

    scp = sweep(scp, 
                i = "proteins", ## assay to be normalized
                MARGIN = 2,
                FUN = "-",
                STATS = colMedians(assay(scp[["proteins"]]), na.rm = TRUE), ## *****Can be changed by different methods for normalization (eg, quantile?)
                name = "proteins_norm_col")

    scp = sweep(scp,
                i = "proteins_norm_col",  ## after col normalization do row normalization 
                MARGIN = 1,
                FUN = "-", ### initial : /
                STATS = rowMeans(assay(scp[["proteins_norm_col"]]),  na.rm = TRUE),
                name = "proteins_norm")
    if(filt_pro_NA){
      cat(paste0("Filtering out sparse protein that over ",protein_pNA*100,"% is NA....."),fill = T)
      scp = filterNA(scp,
                     i = "proteins_norm",
                     pNA = protein_pNA)
    }else{
      cat("Skip filteringout sparse proteins",fill= T)
    } 
    scp = sweep(scp,
                 i = "proteins_norm_col",  ## after col normalization do row normalization 
                 MARGIN = 1,
                 FUN = "-",
                 STATS = rowMeans(assay(scp[["proteins_norm_col"]]),  na.rm = TRUE),
                 name = "proteins_norm_unimputed")
    scp = filterNA(scp,
                i = "proteins_norm_unimputed",
                pNA = protein_pNA)
    cat("Now remove batch effect by Combat & limma...",fill = T)
    require(limma)
    sce = getWithColData(scp, "proteins_norm_unimputed")
    batch = colData(sce)$Raw.file
    assay(sce) = removeBatchEffect( assay(sce),
                                    batch = batch)
    scp = addAssay(scp,
                   y = sce,
                   name = "proteins_batchL")
    scp = addAssayLinkOneToOne(scp, 
                               from = "proteins_norm_unimputed",
                               to = "proteins_batchL")
    cat("Done! Saving files...",fill = T)
    plot(scp,interactive = T)
    save(scp,file = "analysis/Step3_protein_aggregated_normed_scp.RData",compress = T)
    scp.stats(scp)
    cat(dims(scp))
    assign("scp",scp,envir = .GlobalEnv)
    future::plan("sequential")
}



GSEA_Analysis = function( diff_data, gsea_pvalue = 0.05 ){
  
  gene_set = read.gmt("GSEA_imput_files/merge.gmt")
  
  gsea_data = diff_data[,which(colnames( diff_data )%in% colnames_needed)]
  gsea_data = rownames_to_column(gsea_data,"SYMBOL")
  # transfer_id = bitr( gsea_data$SYMBOL,
  #                      fromType = "SYMBOL",toType = "ENTREZID",
  #                      OrgDb = GO_database)
  # gsea_data = merge( gsea_data, transfer_id,by = 'SYMBOL')
  gsea_data = drop_na(gsea_data)
  
  gsea_input = gsea_data$log2FoldChange
  names(gsea_input) = gsea_data$SYMBOL
  gsea_input = sort( gsea_input,decreasing = T )
  
  gsea_output = GSEA( gsea_input,TERM2GENE = gene_set,verbose = F,pvalueCutoff = gsea_pvalue)
  head(gsea_output)
  return(gsea_output)
}

GO_Plot = function( GO_dataset, padj_threshold=0.05, show_num = 15, 
                    fill_color = "#e63f46",wid = 10, hgt = 8,
                    fontsize = 2,fontcolor = "black",
                    barwidth = 0.8,
                    GO_term = character(),
                    keywords = character(),
                    discard = character()){
  tmp = GO_dataset@result
  tmp = tmp[which(tmp$p.adjust<padj_threshold),] 
  tmp$logpval = (-log10( tmp$pvalue ))
  tmp  = tmp %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
  
  if(!is_empty( GO_term )){
    tmp = tmp %>% filter( ONTOLOGY %in% GO_term )
  }
  if(!is_empty( keywords )){
    seek = numeric()
    for( i in keywords){
      j = grep(i,tmp$Description,value = T)
      seek = c(seek,j)
    }
    if(length(seek)==0){return("No aim GO terms was found!")}
  }
  if(show_num==0){
    all_term = unique(seek)
  }else {
    all_term = unique(c(tmp$Description[1:min(length(tmp$Description),show_num)],seek))
  }
  if(length(discard)!=0){
    remove = numeric()
    for( i in discard){
      j = grep(i,tmp$Description,value = T)
      remove = c(remove,j)
    }
    if(length(remove)==0){return("No discarded GO terms was found!")}
    tmp = tmp[tmp$Description %in% setdiff(tmp$Description,remove),]
  }
  
  tmp = tmp[tmp$Description %in% all_term,] %>% arrange( dplyr::desc(logpval),dplyr::desc(Count)) 
  tmp$Description = factor(tmp$Description,levels = rev(tmp$Description))
  go = ggplot()+
    geom_bar( data=tmp,
              aes_string(x = "Description" ,y = "logpval"),
              fill = fill_color,stat = "identity",width = barwidth) + ylab( "-log(p value)" ) +
    geom_text( data=tmp,aes(x = Description ,y = 0.1, label = Description),
               hjust = 0,size = fontsize,color = fontcolor,family = "sans" )+
    xlab(NULL)+scale_y_continuous(expand = c(0,0)) + 
    coord_flip( )+
    # facet_grid( drop = T,rows = "ONTOLOGY",scales = "free",space = "free_y")+
    theme_cowplot()+
    theme( text = element_text( size = 30),
           axis.title = element_text(size = 30),
           axis.text.y = element_blank(),
           axis.text.x= element_text(size = 30,family = "sans"),
           axis.line = element_line(size = 1.5),
           axis.ticks = element_line(size = 1.5),
           axis.ticks.y = element_blank()
           
    )
  return(go)
}


# pep = assay(scp.blood[["peptides"]]) %>% as.data.frame() %>% 
#       rownames_to_column( "peptide" )
# na_tbl = apply(pep[,-c(1)],2,function(x){
#   na_num = length(x[is.na(x)])
#   num = length(x[!is.na(x)])
#   return(c(na_num,num))
# }) %>% t %>% as.data.frame() %>% rownames_to_column("Channel")
# colnames(na_tbl) = c("Channel","NA_number","Valid_number")


# pep = pep %>% gather( key = "Channel",value = "RI", -peptide )
# meta = metadata[,c("Raw.file","tmt.label","celltype")] %>% mutate(Channel = paste0(Raw.file,tmt.label)) %>%
#         select( Channel,celltype ) 
# meta = merge(meta,na_tbl,by = "Channel")
# pep = merge( pep, meta, by = c("Channel"),all.x = T )
# pep$batch = gsub("RI\\d+","",pep$Channel) %>% factor(.)
# pep$Channel = gsub("\\d_\\d+","",pep$Channel)%>% factor(.,levels = paste0("RI",1:16),order = T)
# pep$TMTlabel = ifelse(pep$Channel=="RI1","126",
#                       ifelse(pep$Channel=="RI2","127N",
#                       ifelse(pep$Channel=="RI3","127C",
#                       ifelse(pep$Channel=="RI4","128N",
#                       ifelse(pep$Channel=="RI5","128C",
#                       ifelse(pep$Channel=="RI6","129N",
#                       ifelse(pep$Channel=="RI7","129C",
#                       ifelse(pep$Channel=="RI8","130N",
#                       ifelse(pep$Channel=="RI9","130C",
#                       ifelse(pep$Channel=="RI10","131N",
#                       ifelse(pep$Channel=="RI11","131C",
#                       ifelse(pep$Channel=="RI12","132N",
#                       ifelse(pep$Channel=="RI13","132C",
#                       ifelse(pep$Channel=="RI14","133N",
#                       ifelse(pep$Channel=="RI15","133C","134N"))))))))))))))) %>%
#                       factor(.,levels =c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N"),ordered = T)
# pep$celltype = factor(pep$celltype,levels = c("Carrier","Reference","Unused","Mono","Macro","Blank"))

# p3 = ggplot( pep, aes( x = TMTlabel, y = log10(RI), fill = celltype ) ) +
#   geom_boxplot(outlier.size = 1) + scale_fill_brewer( palette = "Set3" ) + 
#   theme_cowplot(font_size = 26,font_family = "sans")   +
#   xlab(NULL)+ 
#   # ggtitle( label = scp.stats(scp.blood) ) +
#     theme(text = element_text(size = 26),
#         axis.text.x = element_text( angle = 45,size = 14, hjust = 1, vjust = 1  ),
#         title = element_text(size = 26),
#         axis.ticks = element_line(size = 1.5),
#         axis.title = element_text(size = 26),
#         axis.line = element_line(size = 1.5),
#         strip.background = element_blank(),strip.text = element_text(size = 18)) + facet_wrap(~batch, ncol = 1)
# p3

# meta = meta %>% gather( key = "Type",value = "Number",-celltype,-Channel )
# meta$celltype = factor(meta$celltype,levels = c("Carrier","Reference","Unused","Mono","Macro","Blank"))
# meta$batch = gsub("RI\\d+","",meta$Channel) %>% factor(.) %>% factor(.)
# meta$Channel =  gsub("\\d_\\d+","",meta$Channel) %>% 
#               factor(.,levels = paste0("RI",1:16),order = T)
# p4 = ggplot( meta[meta$batch!="1_12",],aes(x = Channel, y = Number/1000, fill = Type ) ) + 
#   geom_bar(position = "dodge",stat = "identity") + 
#   scale_fill_manual(values = c("#8f8e8e","#eb6059")) +
#   theme_cowplot(font_size = 26,font_family = "sans")  + 
#   xlab(NULL) + ylab(paste0("Pepetide counts (x1000)")) +
#   # ggtitle( label = scp.stats(scp.blood) ) +
#   theme(text = element_text(size = 26),
#         axis.text.x = element_text( angle = 45,size = 26, hjust = 1, vjust = 1  ),
#         title = element_text(size = 26),
#         axis.title = element_text(size = 26),
#         axis.line = element_line(size = 1.5),
#         axis.ticks = element_line(size = 1.5),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 18)) + facet_wrap(facets = ~batch,ncol = 1)
# p4



# library(VennDiagram)
# library(venneuler)
# library(eulerr)


# fit1 <- euler(c("Mustache"=5330,  ### 1 = mustache,2  = juicer, 3 =  homer
#                 "Juicer"=685,
#                 "Homer"=1426,
#                 "Mustache&Juicer" = 1888,
#                 "Mustache&Homer" = 546,
#                 "Juicer&Homer" = 99,
#                 "Mustache&Juicer&Homer" = 2263))
# pdf("/lustre/user/liclab/fangq/proj/embryo/hic/analysis/loop/cis/three_methods_comp/venn.ko.pdf",width = 8,height = 6)
# # draw.triple.venn(area1=10025,  ### 1 = mustache,2  = juicer, 3 =  homer
# #                  area2=4956,
# #                  area3=4319,
# #                  n12 = 1752+2070,
# #                  n23 = 129+2070,
# #                  n13 = 498+2070,
# #                  n123 = 2070,
# #                  category = c("Mustache", "Juicer", "Homer"),
# #                  fill=c("#259645", '#e53d30',"#f9c01c"),
# #                  alpha = c(0.5,0.5,0.5),
# #                  ext.dist = -0.3, # color for each circle,
# #                  lwd=c(0,0,0),cat.cex = 2,cex = 2,fontfamily = "sans",lty = "blank",
# #                  cat.fontfamily = "sans",scaled = T  # line width for each circle 
# # )
# plot(fit1,
#      quantities = list(cex = 2),# TRUE,
#      fill=alpha(c("#259645", '#e53d30',"#f9c01c"),alpha = 0.5),#"transparent",
#      # lty = 1:3,
#      shape = "ellipse",fontsize = 10,
#      legend = list(labels = c("Mustache", "Juicer", "Homer"),cex = 2.5),
#      # labels = list(font = 1,cex = 4,just = "bottom",hjust = 1,vjust = 1,check.overlap = T),
#      adjust_labels = TRUE)
# dev.off()