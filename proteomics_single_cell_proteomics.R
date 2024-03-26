setwd("/lustre/user/liclab/lisky/zhanggx/scp/sc_merge/")
source("/lustre/user/liclab/lisky/zhanggx/scp/sc_merge/R_pack/scp_functions.r")#This is the second R script, the content is at the bottom of this Rscript
#load packages needed
pkgs =  c( 'jsonlite', 'languageserver','IRanges','GenomicRanges','ComplexHeatmap',
           'pheatmap', 'data.table', 'tidyverse', 'magrittr', 'RColorBrewer', 'cowplot', 'patchwork', 'ggridges',
           'scales', 'textshape', 'stringr',"scuttle","SingleCellExperiment",
           'impute', 'scater', 'sva', 'GGally','limma','QFeatures', 'scp', 'Seurat',
           'scpdata','eulerr','ggplotify','GGally','ComplexHeatmap')
Loadpkgs = function(i){
  print(paste("Installing package:", i))
  if(!require(i, quietly = TRUE)){
    BiocManager::install(i, update = FALSE)
  }
  library(i, character.only = TRUE)
}
sapply( pkgs,Loadpkgs )

cellnames = c("HeLa","U937")
TMTlabel = factor(c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N"),
                  levels = c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N"),ordered = T)

#load the data faster
evidence = fread("/lustre/user/liclab/lisky/zhanggx/scp/sc_merge/dart_id/ev_updated.txt",nThread = 10,header = T)
colnames(evidence)=gsub(" ",".",colnames(evidence))
colnames(evidence) = gsub("Reporter.intensity.corrected.","RI",colnames(evidence))

evidence$Raw.file = gsub("2023","",evidence$Raw.file)

sc_batch = unique(evidence$Raw.file) 
evidence = evidence %>% filter( Raw.file %in% sc_batch )

nBatchs = sc_batch %>% length
metadata = data.frame( Raw.file = rep(sc_batch,each = 16),
                       channel = rep(factor(paste0("RI",1:16)),nBatchs),
                       TMTlabel = rep(TMTlabel,nBatchs)) %>% arrange(Raw.file,TMTlabel)

metadata$celltype = ifelse(metadata$Raw.file %in% grep(metadata$Raw.file,pattern = '[B|A]\\d+',value=T),
                           ifelse(as.numeric(str_extract(metadata$Raw.file,pattern = '(?<=[B|A])\\d+'))>=13,
                                  c("Carrier","Reference","Unused","Unused",rep(c("HeLa","U937"),5),"HeLa","Blank"),
                                  c("Carrier","Reference","Unused","Unused","Blank",rep(c("U937","HeLa"),5),"U937")),
                           ifelse(metadata$Raw.file %in% c("wAP329","wAP330","wAP331"),
                                  c("Carrier","Reference","Unused","Unused",rep(c("HeLa","U937"),6)),
                                  c("Carrier","Reference","Unused","Unused",rep(c("U937","HeLa"),6))))



metadata$Batch = ifelse( metadata$Raw.file %in% grep(sc_batch,pattern = "B",value = T),"B",
                         ifelse(metadata$Raw.file %in% grep(sc_batch,pattern = "A\\d+",value = T),"A","Public"))
save(evidence,metadata,file = "Step0_raw_meta.RData",compress = T)


####get the AA length####

ggplot(evidence) + 
  geom_histogram(aes(x = Length),binwidth = 1) +
  scale_x_continuous(limits = c(5,51),breaks = c(5,6,7,10,12,30,50),expand = c(0,0)) + 
  scale_y_continuous( limits = c(0,62000)) + 
  Themes(rotate = F)

######scp analysis#######

pif = 0.8
pvalue = 0.01
meanscr = 0.1
medianCV = 0.4

cols = brewer.pal(name = "Set3",n=6)
names(cols) = c("Carrier","Reference","Unused","HeLa","U937","Blank")

Step1_generate_raw(ev_data = evidence,meta_data = metadata,
                   batchCol = "Raw.file",channelCol = "channel",
                   filt_with_min_psm = T,psm_threshold = 100,
                   samplePattern = "HeLa|U937",sc_batch = sc_batch)
load("Step1_raw_scp.RData")
RI_plot = function(scp = scp, assay = character()){
  tmp = assay(scp,i = assay) %>% 
    as.data.frame() %>% #drop_Nas(.,dims = 1) %>% 
    rownames_to_column("Peptides") %>% 
    gather(key = Channels, value = RI, -Peptides )
  coldata = colData(scp) %>% as.data.frame() %>% rownames_to_column("Channels")
  tmp = merge(tmp,coldata,by = "Channels")
  p = ggplot(tmp,aes( x = TMTlabel, y = log10(RI), fill = celltype )) +
    geom_violin(size = 1,width = 1) + 
    scale_y_continuous( limits = c(0,6)) +

    stat_summary(fun.data = data_summary,linewidth = 0.7,
                 color = "black",fatten = 3,geom = "pointrange") +
    scale_fill_manual( values = cols) + 
    scale_color_manual( values = cols ) + 
    ylab("log10(medianRI)") + xlab(NULL) + ggtitle(assay) + 
    Themes(pt_size = 2,axis_x_fontsize = 16,axis_y_fontsize = 20,fontsize = 18,titlefontsize = 20,linesize = 1)
  return(p)
}
### Get the relative intensity among single cell channels to see if any batches are abnormal###
scb_B = grep(sc_batch,pattern="B",value = T)
scb_H = grep(sc_batch,pattern="A",value = T)
scb_pub = grep(sc_batch,pattern="^w",value = T)
p = RI_plot(scp,assay = scb_B[1])
for(i in c(scb_B[-1],scb_H,scb_pub)){
  p_tmp = RI_plot(scp,assay = i)
  p = p + p_tmp
  rm(p_tmp)
}

pdf("figures/P1_raw_RI_sc_batch.pdf",width = 35,height = 30)
p + plot_layout(guides = "collect",ncol = 5)
dev.off()

### Check the MeanSC###
Stats = data.frame(dims(scp))
scp.tmp = filterFeatures(scp, ~ Reverse != "+" & Potential.contaminant != "+" & !grepl("REV|CON", Leading.razor.protein))
Stats = rbind(Stats,data.frame(dims(scp.tmp)))
scp.tmp = filterFeatures(scp.tmp, ~ PIF > pif & !is.na(PIF))
Stats = rbind(Stats,data.frame(dims(scp.tmp)))
scp.tmp = filterFeatures(scp.tmp, ~ qvalue_proteins < pvalue &  qvalue_PSMs < pvalue  )
Stats = rbind(Stats,data.frame(dims(scp.tmp)))
scp = scp.tmp ; rm(scp.tmp)
tmp = rbindRowData(scp,sc_batch) %>% as.data.frame() 
options(scipen = 0)

p1 = ggplot( tmp, aes( x = MeanSCR,y = Raw.file, group = Raw.file ) ) +
  geom_density_ridges(fill = "#8db5de",bandwidth = 0.1,size = 1,scale = 2,rel_min_height = .05) +
  scale_x_log10(breaks = c(1e-3,1e-2,1e-1,1),limits = c(1e-4,1e0),expand = c(0,0)) +
  geom_vline( xintercept = meanscr,linewidth =1,lty = 2) +
  ylab("Batch")+ 
  Themes(pt_size = 14,titlefontsize = 20,
         axis_x_fontsize = 20,axis_y_fontsize = 20,
         fontsize =22 ,rotate = T,linesize = 1)  
pdf("figures/P2_MeanSCR.pdf",width = 6,height = 10)
p1
dev.off()

scp = QFeatures::filterFeatures( scp, ~ !is.na(MeanSCR) & MeanSCR < meanscr )

Stats = rbind(Stats,data.frame(dims(scp)))
save(scp,Stats,file = "Step1-1_pep_filt_scp.RData",compress = T)
load("Step1-1_pep_filt_scp.RData",verbose = T)


### Turn pSM into peptides###

Step2_psm_to_pep(scp=scp,Raw.file_feature = "^\\w+",
                 if_generate_psm = F,
                 medianCV_batches = sc_batch,
                 sc_batch = sc_batch) ### may take long time
# load("Step2_pep_aggregated_scp.RData",verbose = T)

### Get the Median_RI###

medianRI = colMeans( assay(scp[["peptides"]]),na.rm = T )
scp$medianRI = medianRI

k = getWithColData(scp, "peptides") %>%
  colData %>%
  data.frame
k$celltype = factor(k$celltype,levels = rev(c("Carrier","Reference","Unused","HeLa","U937","Blank")))
p2 = ggplot(k,aes( x = log10(medianRI), y = celltype, fill =  celltype )) +
  geom_boxplot(size = 0.5) + 
  scale_fill_manual( values = cols ) + 
  scale_x_continuous(limits = c(-1,7)) +
  ggtitle("MedianRI") +facet_grid(Batch~.,drop = T,scales = "free",space = "free") +
  ggtitle(NULL) + xlab("log10(medianRI)") + ylab(NULL) + 
  Themes(pt_size = 2,fontsize = 20,
         axis_x_fontsize = 18,axis_y_fontsize = 18,
         titlefontsize = 20,linesize = 1,rotate = F)+ theme(legend.position = "NULL")
pdf("figures/P3_MedianRI.pdf",width = 5,height = 6)
p2
dev.off()

### Get the MedianCV， and some strict QC###

a = summary(k$MedianCV[k$celltype %in% c("Carrier","Reference")])
a
p3 = k %>%
  ggplot(aes(x = MedianCV,y = celltype,fill = celltype,group = celltype)) +
  geom_density_ridges(size = 1,scale = 2,rel_min_height = 1e-3,bandwidth = 0.015) + 
  scale_fill_manual(values = alpha(cols,alpha = 0.6)) +
  geom_vline(xintercept = c(0.4),
             size = 1,linetype = 2) +
  xlab("MedianCV")+ ylab(NULL) + facet_grid(Batch~.,drop = T,scales = "free",space = "free") +
  Themes(fontsize = 20, 
         axis_x_fontsize = 20,axis_y_fontsize = 20,
         titlefontsize = 18,linesize = 1,rotate = F) + theme(legend.position = "NULL") 
pdf("figures/P4_medianCV.pdf",width = 6,height = 7)
p3
dev.off()
medianCV = 0.4

exclusion = setdiff(unique(k$Raw.file),unique(k$Raw.file[k$MedianCV<medianCV]))
if(!is_empty(exclusion)){
  tbl = matrix(NA,ncol = length(exclusion),nrow = 2,dimnames = list(NULL,exclusion))
  exclusion = c(exclusion,paste0("peptides_",exclusion))
  keepAssay = names(scp)[!(names(scp) %in% exclusion)]
  scp = scp[,,keepAssay]
  scp = scp[,!is.na(scp$MedianCV) & scp$MedianCV < medianCV & scp$celltype %in% cellnames, ]
  tbl = cbind(tbl,dims(scp) %>% as.data.frame() %>%  select(c(intersect(names(scp),unique(metadata$Raw.file)))))
}else{
  scp = scp[,!is.na(scp$MedianCV) & scp$MedianCV < medianCV & scp$celltype %in% cellnames, ]
  tbl = dims(scp) %>% as.data.frame() %>%  select(c(intersect(names(scp),unique(metadata$Raw.file))))
}
colnames(tbl) = ifelse(colnames(tbl) %in% grep(colnames(tbl),pattern = '^\\d',value = T),paste0("X",colnames(tbl)),colnames(tbl))
Stats = rbind(Stats,tbl)

Stats$Steps = factor(rep(c("Raw","Contaminant",paste0("PIF",pif),paste0("q_val",pvalue),paste0("MeanSCR",meanscr),paste0("MedianCV",medianCV)),each = 2),
                     levels = c("Raw","Contaminant",paste0("PIF",pif),paste0("q_val",pvalue),paste0("MeanSCR",meanscr),paste0("MedianCV",medianCV)))
Stats$Data = rep(c("PSMs","Cells"),times = 6)

### For the PSM in each steps ####
tmp = Stats %>% 
  gather( key = "Raw.file",value = "Number", -Steps,-Data )
tmp$Raw.file = gsub("X","",tmp$Raw.file)
tmp = tmp %>% filter(Raw.file %in% c(sc_batch)) %>% 
  merge(.,metadata[,c("Raw.file","Batch")], by = "Raw.file",all.x = T, all.y = F) %>% 
  unique() 
tmp$Number[tmp$Data=="Cells"&tmp$Steps!="Median0.4"&tmp$Batch!="Public"] = tmp$Number[tmp$Data=="Cells"&tmp$Steps!="Median0.4"&tmp$Batch!="Public"]  - 5
tmp$Number[tmp$Data=="Cells"&tmp$Steps!="Median0.4"&tmp$Batch=="Public"] = tmp$Number[tmp$Data=="Cells"&tmp$Steps!="Median0.4"&tmp$Batch=="Public"]  - 4
tmp = tmp %>% dplyr::group_by(Steps,Data,Batch) %>% dplyr::mutate(median = median(Number,na.rm = T))

p4 = ggplot( tmp ) + 
  geom_bar(aes(x = Steps, y = median, group = Steps,fill = Steps),width = 0.8,
           stat = "identity",position = "dodge") +
  scale_fill_manual( values = alpha(brewer.pal("YlOrRd",n=6)),1 ) + 
  stat_summary(aes(x = Steps, y = Number, group = Steps),
               fun.data = data_summary,linewidth = 0.7,
               color = "black",fatten = 2,
               geom = "pointrange") +
  Themes(pt_size = 8,fontsize = 18,linesize = 1,
         axis_x_fontsize = 16,axis_y_fontsize = 18,
         titlefontsize = 22,rotate = T) + 
  xlab(NULL) + ylab( "Number" ) + facet_grid(Data~Batch,scales = "free_y",space = "free_x") + theme(legend.position = "none")

pdf("figures/P5_statistics.pdf",width = 10,height = 6)
p4
dev.off()
save(scp,Stats,file = "Step2-1_cell_filt_scp.RData",compress = T)


### Turn peptide into protein ###

Step3_normalization_pep_aggre_protein(scp = scp, pep_pNA = 0.99,protein_pNA = 0.99,filt_pep_NA = F,filt_pro_NA = T)
plot(scp,interactive = T)
# load("Step3_protein_aggregated_normed_scp.RData")
scp[["proteins_norm"]] %>%
  assay %>% 
  is.na %>%  mean()
scp = impute(scp,
             i = "proteins_norm",
             method = "knn",
             k = 5, rowmax = 1, colmax= 1,
             maxp = Inf, rng.seed = 1234 )
scp[["imputedAssay"]] %>%
  assay %>%
  is.na %>%
  mean #### 0
save(scp,file = "Step3-1_imputed_scp.RData",compress = T)

sce = getWithColData(scp, "imputedAssay")
batch = colData(sce)$Raw.file 
model = model.matrix( ~ celltype, data = colData(sce))  # before combat, excluded the intersted variables

## Combat  bacth corrected
assay(sce) = ComBat(dat = as.matrix(assay(sce)),
                    batch = batch,
                    mod = model)
scp = addAssay(scp,
               y = sce,
               name = "proteins_batchC")
scp = addAssayLinkOneToOne(scp, 
                           from = "imputedAssay",
                           to = "proteins_batchC")
plot(scp,interactive = T)
save(scp,file = "Step3-2_batch_correct_scp.RData",compress = T)
load("Step3-2_batch_correct_scp.RData")

#### correlate with Nikolai dataset ### 

a = getWithColData(scp,"proteins_norm")
# 
# hl = a[rownames(a) %in% c(pros),a$celltype=="U937"]
# pearscor = cor(assay(hl),method = "pearson",use = "na.or.complete")
pearscor = cor(assay(a),method = "pearson",use = "na.or.complete")
cell_cols = alpha(c("#6097b0","#ed7057"),0.8)
names(cell_cols) =cellnames
cluster_cols = alpha(c("#d65240","#487fe7","#549f55"),0.9)
names(cluster_cols) = unique(a$Batch)

top = HeatmapAnnotation(Batch = a$Batch,
                        Cell = a$celltype,
                        col = list(Batch = cluster_cols,
                                   Cell = cell_cols),
                        show_legend = T,
                        show_annotation_name = F,
                        annotation_legend_param = list(labels_gp  = gpar(fontsize = 15),
                                                       title_gp = gpar(fontsize = 15),
                                                       legend_direction = "horizontal")
)
pearscor[pearscor>0.8] = 0.8
pearscor[pearscor<(-0.8)] = -0.8
pdf("figures/P6_cell_corr_heatmap.pdf",width = 10,height = 8)
Heatmap( pearscor, 
         heatmap_legend_param = list( title = "Pearson",
                                      at = c(-0.8,-0.4,0,0.4,0.8),
                                      labels = c(-0.8,-0.4,0,0.4,0.8),
                                      labels_gp  = gpar(fontsize = 15),
                                      direction = "vertical",
                                      title_position = "topcenter",
                                      title_gp = gpar(fontsize = 15)),
         show_heatmap_legend = T,
         col = circlize::colorRamp2(breaks = seq(-0.8,0.8,length.out = 9),
                                    colors = colorRampPalette(c("blue","white","red"))(9),#colorRampPalette(c("#00aeff","Black","#f9d31a"))(9),
                                    transparency = 0),
         cluster_columns = T,
         cluster_rows = T ,
         show_column_names = F,
         clustering_method_rows = "complete",#ward.D2
         clustering_distance_rows = "euclidean",
         show_row_names = F,
         top_annotation = top,
         column_names_rot = 45,
         width = unit(7,"in"),height = unit(7,"in"),
         column_names_gp = gpar(fontsize = 15) )
dev.off()

##### protein number per cell
tmp = getWithColData(scp,"proteins_norm")
pro_valid = apply( assay(tmp),2, function(x){
  Valid = length(x[!is.na(x)])
}) %>% as.data.frame() 

colnames(pro_valid) = c("pro_num")
pro_valid$Batchnew = ifelse(pro_valid$Batch=="Public",'Public',"This Study")
pro_valid = cbind(pro_valid,scp@colData)
p5 = ggplot(pro_valid) + 
  geom_density_ridges(aes(y = Batchnew,x = pro_num, fill = celltype,color = celltype),
                      jittered_points = TRUE, scale = 1, rel_min_height = .01,size = 0.5,
                      point_shape = "|", point_size = 4, 
                      position = position_points_jitter(height = 0),
                      bandwidth = 3) +
  scale_y_discrete(expand = c(0.1,0)) +
  scale_fill_manual(values = alpha(c("#6097b0","#ed7057"),0.5)) +
  scale_color_manual(values = alpha(c("#6097b0","#ed7057"),1)) +
  ylab("Batch") + xlab( "Protein Number" ) + 
  Themes(pt_size = 4,fontsize = 25,
         axis_x_fontsize = 26,axis_y_fontsize = 26,
         linesize = 1.5,rotate = F)
p5
pdf("figures/P7_protein_number.pdf",width = 8,height = 3)
p5
dev.off()


### Get the Venn to see protein overlap###

names(pro_density) = pro_density = unique(scp$Raw.file)
pro_density = lapply( as.list(pro_density), function(x){
  p = assay(tmp[,tmp$Raw.file==x]) %>% drop_Nas(.,dims = 1,partitial = F)
  Valid_density = apply(p,1,function(i){
    return(length(i[!is.na(i)])/length(i))
  }) %>% as.data.frame() %>% rownames_to_column( .,"Proteins" ) 
  colnames(Valid_density) = c("Proteins","pro_density")
  return(Valid_density)
}) %>% plyr::ldply(.,.id = "Raw.file")

pro_density = merge(pro_density,metadata[,c("Raw.file","Batch")],by = "Raw.file")

a_pros = unique(pro_density$Proteins[pro_density$Batch=="A"])
b_pros = unique(pro_density$Proteins[pro_density$Batch=="B"])
pub_pros = unique(pro_density$Proteins[pro_density$Batch=="Public"])


pdf("figures/P8_protein_overlap.pdf",width = 5,height = 5)

### the protein heat map ###

protein_detect =  assay(tmp) #%>%
protein_detect[!is.na(protein_detect)] = 1
protein_detect[is.na(protein_detect)] = 0

protein_detect = as.data.frame(protein_detect) %>% rownames_to_column(.,"Proteins") %>% gather(.,key = "Cells",value = "is_detect",-Proteins)
n = as.data.frame(tmp@colData) %>% rownames_to_column(.,"Cells") %>% select(Raw.file,Cells,celltype,Batch)
protein_detect = merge(n,protein_detect,by = "Cells") %>% 
  dplyr::group_by(Raw.file,Proteins,celltype) %>% 
  mutate( Detect_density = sum(is_detect)/length(is_detect) ) %>% 
  select(-Cells,-is_detect) %>% unique()
k = protein_detect %>%
  dplyr::mutate( tmp = paste0(Raw.file,"_",celltype) ) 
anno = k[,c(1:3,6)] %>% unique
k = k[,-c(1:3)] %>% spread( .,key = tmp,value = Detect_density) %>% column_to_rownames(.,"Proteins")

cell_cols = alpha(c("#6097b0","#ed7057"),0.8)
names(cell_cols) =cellnames
cluster_cols = alpha(c("#d65240","#487fe7","#549f55"),0.9)
names(cluster_cols) = unique(anno$Batch)

top = HeatmapAnnotation(Batch = anno$Batch,
                        Cell = anno$celltype,
                        col = list(Batch = cluster_cols,
                                   Cell = cell_cols),
                        # annotation_label = c("Cell","Batch"),
                        show_legend = T,
                        show_annotation_name = F,
                        annotation_legend_param = list(labels_gp  = gpar(fontsize = 15),
                                                       title_gp = gpar(fontsize = 15))
)
pdf("figures/P9_pro_detect_heatmap.pdf",width = 7,height = 8)
Heatmap( k, 
         heatmap_legend_param = list( title = "Detect Density",
                                      at = c(0,0.5,1),
                                      labels = c(0,0.5,1),
                                      labels_gp  = gpar(fontsize = 15),
                                      direction = "vertical",
                                      title_gp = gpar(fontsize = 15)),
         show_heatmap_legend = T,
         col = viridis_pal(option = "turbo",begin = 1,end = 0)(20),
         # col = rev(colorRampPalette(c("Red","Blue","Black"))(10)),
         cluster_columns = T,
         cluster_rows = T ,
         show_column_names = F,
         clustering_method_rows = "complete",#ward.D2
         clustering_distance_rows = "euclidean",
         #("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
         show_row_names = F,
         top_annotation = top,
         # left_annotation = left,
         column_names_rot = 45,
         width = unit(4,"in"),height = unit(6,"in"),
         column_names_gp = gpar(fontsize = 15) )
dev.off()

### bar plot of protein detected###

a = protein_detect %>% spread( .,key = celltype,value = Detect_density) 
a$Protrein_type = ifelse( (a$HeLa+a$U937)==2,"All detect",
                          ifelse( (a$HeLa+a$U937)==0,"Non detect",
                                  ifelse( a$HeLa>=0.5&a$U937<=0.4,"HeLa specific",
                                          ifelse( a$U937>=0.5&a$HeLa<=0.4,"U937 specific","Variable"))))
a$Protrein_type = factor(a$Protrein_type,levels = rev(c("All detect","Variable","HeLa specific","U937 specific","Non detect")),ordered = T)
a = a[a$Protrein_type!="Non detect",]
pdf("figures/P9_pro_detect_bar.pdf",width = 12,height = 5)
ggplot(a) +
  geom_bar( aes(x = Raw.file,fill = Protrein_type ),stat = "count",position = position_stack(),width = 0.9 ) + 
  scale_fill_manual(values = c("#d65240","#54a056","#f1be44","skyblue")) + 
  xlab(NULL) + ylab('Protein Number') + 
  annotate(geom = "segment",color = "#e97e36",x=0.55,xend = 13.45,y = -10,yend = -10,size = 3)+
  annotate(geom = "segment",color = "Navy",x=13.55,xend = 21.45,y = -10,yend = -10,size = 3)+
  annotate(geom = "segment",color = "Red",x=21.55,xend = 27.45,y = -10,yend = -10,size = 3)+
  Themes(axis_x_fontsize = 13,fontsize = 18,axis_y_fontsize = 20,linesize = 1)
dev.off()


##### PCA reduction analysis

scp.tmp = getWithColData(scp,"proteins_batchC")#proteins_batchC
scp.tmp = runPCA(scp.tmp,
                 ncomponents = 10, 
                 ntop = Inf,
                 scale = TRUE,
                 exprs_values = 1,
                 name = "PCA")
p = plotReducedDim(scp.tmp,
                   dimred = "PCA",
                   colour_by = "celltype",
                   point_alpha = 1, text_size = 20,point_size = 4 ) + 
  theme_cowplot() 
p
tmp = getWithColData(scp,"proteins_batchC")#imputedAssays
mat = assay(tmp) #%>% drop_Nas(dims = 1,partitial = T,valid_proportions = 0.9)
# mat = mat[,grep(colnames(mat),pattern = "H",value = T,invert = T)]rm 

pca = prcomp( t(mat %>% t %>% scale(.,center = T,scale = T) %>% t ))
deviation = (pca$sdev)^2/sum((pca$sdev)^2)*100
sum(deviation[1:2]) ;length(deviation)
anno = tmp@colData %>% as.data.frame() #%>% filter(Batch!="H")
pca_mat = cbind(pca$x,anno)
pca_mat$celltype = factor(pca_mat$celltype,levels = cellnames)
pca_mat$Raw.file = as.factor(pca_mat$Raw.file)

p8 = ggplot(pca_mat,aes(x = PC1,y = PC2)) + 
  geom_point(aes(color = celltype,shape = Batch),size =5) + 
  scale_shape_manual(values = c(15:21,8,13)) +
  scale_color_manual(values = alpha(c("#6097b0","#ed7057"),0.6)) +
  xlab(paste0("PC1 (",round(deviation[1],digits = 1),"%)")) +
  ylab(paste0("PC2 (",round(deviation[2],digits = 1),"%)")) +
  Themes(pt_size = 4,fontsize = 20,
         axis_x_fontsize = 25,axis_y_fontsize = 25,
         linesize = 1.5,rotate = F)
p8
pdf("figures/P10_raw_PCA.pdf",width = 7,height = 5)
p8
dev.off()



##### only single cells, get the protein list , feature plot and GO###
{
  normedassay0 = getWithColData(scp, "proteins_norm")
  seur = as.Seurat( normedassay0,
                    counts = "aggcounts", data = "assay" )
  seur = RenameAssays(seur,originalexp = "proteins_norm")
  
  normedassay1 = getWithColData(scp, "proteins_batchL")
  add = as.Seurat( normedassay1,
                   counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["unimputed_raw_limma"]] = CreateAssayObject( data = ttt )
  
  assay(normedassay1,"assay") = assay(normedassay1,"assay") %>% t() %>% apply( .,2,function(x){x = (x-mean(x,na.rm = T))/sd(x,na.rm = T)} ) %>% t()
  add = as.Seurat( normedassay1,
                   counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["unimputed_raw_limma_zscore"]] = CreateAssayObject( data = ttt )
  normedassay2 = getWithColData(scp, "imputedAssay")
  add = as.Seurat( normedassay2,
                   counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["imputed_raw"]] = CreateAssayObject( data = ttt )
  assay(normedassay2,"assay") = assay(normedassay2,"assay")  %>% t() %>% apply( .,2,function(x){x = (x-mean(x,na.rm = T))/sd(x,na.rm = T)} ) %>% t()
  add = as.Seurat( normedassay2,
                   counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["imputed_raw_zscore"]] = CreateAssayObject( data = ttt )
  
  normedassay3 = scp[["proteins_batchC"]]
  add = as.Seurat( normedassay3,
                   counts = "aggcounts", data = "assay" )
  ttt = GetAssayData( object = add,slot = "data" )
  seur[["imputed_raw_Combat"]] = CreateAssayObject( data = ttt )

  Idents(seur) = "celltype"
  save(seur,file = "seurat.RData",compress = T)
}
VlnPlot( seur,features = c("nFeature_originalexp"),group.by = "Raw.file")

DefaultAssay(seur)  = "proteins_norm"
seur = FindVariableFeatures(seur,nfeatures = 2000,selection.method = "vst")
seur = ScaleData( seur,do.scale = T,do.center = F )
seur = Seurat::RunPCA( seur,assay = "proteins_norm",npcs = 60)
ElbowPlot(seur,ndims = 60)
DimPlot( seur,reduction = "pca",group.by = c("celltype","Batch"),
         pt.size = 3) & 
  Themes(rotate = F) 

seur = RunUMAP( seur,dims = 1:30,reduction = "pca",n.neighbors = 30,min.dist = 0.05 )#n.neighbors = 30,min.dist = 0.05 
# seur = RunTSNE( seur,dims = 1:20,reduction = "pca" )
DimPlot( seur,reduction = "umap",group.by = c("celltype","Batch"),pt.size = 2 )

# seur = subset(seur,subset = Batch!="H")
DefaultAssay(seur)  = "imputed_raw_Combat"
VariableFeatures(seur)=rownames(seur)
seur = ScaleData( seur,do.scale = T,do.center = F )
seur = Seurat::RunPCA( seur,assay = "imputed_raw_Combat",npcs = 60)
seur = Seurat::RunUMAP( seur,dims = 1:30,reduction = "pca",n.neighbors = 40,min.dist = 0.2 )
seur = Seurat::RunTSNE( seur,dims = 1:30,reduction = "pca")

DimPlot( seur,reduction = "pca",group.by = c("celltype","Batch"),
         pt.size = 2) #& Themes(rotate = F) 

ElbowPlot(seur,ndims = 60)
DimPlot( seur,reduction = "umap",group.by = c("celltype","Batch"),
         # cols = alpha(c("#6097b0","#ed7057"),0.6),
         pt.size = 3) #& Themes(rotate = F) 

DefaultAssay(seur) = "imputed_raw_Combat"
seur = FindNeighbors(object = seur,assay = "imputed_raw_Combat")
cluster.test = FindClusters(object = seur,resolution = seq(0,1,by = 0.2))
clustree(cluster.test@meta.data, prefix = "imputed_raw_Combat_snn_res.") ### resolution = 0.4
seur = FindClusters(object = seur,resolution = 0.8)
Idents(seur) = "celltype"
table(seur$seurat_clusters,seur$celltype)

DefaultAssay(seur)  = 'imputed_raw'
Idents(seur) = "celltype"
seur = ScaleData(seur,do.scale = T,do.center = F)
markers = FindAllMarkers(seur,logfc.threshold = 0.25,min.pct = 0,only.pos = T,return.thresh = 0.05,slot = "data" )
markers = markers[markers$p_val_adj<0.05,]
write.table(rownames(seur),"proteins.txt",quote = F,col.names = F,row.names = F)
library(ensembldb)
pro_to_gene = readxl::read_xlsx("uniprot_2008human.xlsx")
library(clustree)
marker_gene = merge(markers,pro_to_gene,by.x = "gene",by.y = "Entry",all.x = T,all.y = F)

marker_gene
DefaultAssay(seur)  = 'proteins_norm'
seur = ScaleData(seur,do.scale = T,do.center = F)
pa_1 = DimPlot(seur,reduction = "pca",group.by = c('celltype'),pt.size = 2,cols = cell_cols)
pa_2 = FeaturePlot( seur,reduction = "pca",features = c('P05204','O60506','P04083','P21796'),
                    slot = "data",
                    cols = alpha(viridis_pal(option = 'magma')(6),0.9),
                    pt.size = 2)   #& Themes(rotate = F) 
pa_2[[1]]$labels$title = marker_gene$`Gene Names (primary)`[marker_gene$gene==pa_2[[1]]$labels$title]
pa_2[[2]]$labels$title = marker_gene$`Gene Names (primary)`[marker_gene$gene==pa_2[[2]]$labels$title]
pa_2[[3]]$labels$title = marker_gene$`Gene Names (primary)`[marker_gene$gene==pa_2[[3]]$labels$title]
pa_2[[4]]$labels$title = marker_gene$`Gene Names (primary)`[marker_gene$gene==pa_2[[4]]$labels$title]

pdf('figures/P13_featureplot.pdf',width = 6,height = 8)
pa_1 + pa_2 + plot_layout(design = '
                          #AA#
                          BBBB
                          BBBB
                          BBBB
                          ')
dev.off()

library(org.Hs.eg.db)
library(clusterProfiler)
species = "human"
database = function(species = character()){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(ensembldb)
  require(clusterProfiler)
  if(species=="human"){
    ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    orgdb = org.Hs.eg.db
  }else if(species=="mouse"){
    ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
    orgdb = org.Mm.eg.db
  }
  wg.gene = data.frame(genes(ensdb)) %>%  dplyr::filter( gene_biotype == "protein_coding" )
  wg.gene = wg.gene$symbol
  return(list(wg.gene=wg.gene,orgdb=orgdb))
}
go_hl = enrichGO(marker_gene$`Gene Names (primary)`[marker_gene$cluster=="HeLa"],
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",ont = "BP",minGSSize = 5,
                 pAdjustMethod = "BH",universe =database(species)[["wg.gene"]],
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
go_hl_simplified = clusterProfiler::simplify(x = go_hl ,cutoff=0.4,by="p.adjust",select_fun=min)  
go_u9 = enrichGO(marker_gene$`Gene Names (primary)`[marker_gene$cluster=="U937"],
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",ont = "BP",minGSSize = 5,
                 pAdjustMethod = "BH",universe = database(species)[["wg.gene"]],
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
go_u9_simplified = clusterProfiler::simplify(x = go_u9 ,cutoff=0.6,by="p.adjust",select_fun=min)  


p1 = GO_Plot(go_hl,fontsize = 7,fill_color = alpha("#6097b0",0.6),
             padj_threshold = 0.3,show_num = 15,keywords = c("cytokine" ),discard = "DNA")
p1
p2 = GO_Plot(go_u9,fontsize = 7,fill_color = alpha("#ed7057",0.6),
             padj_threshold = 0.3,show_num = 15,keywords = c("leu","immun","cytok" ),discard = "DNA")
pdf("figures/P12_GO.pdf",width = 15,height = 8)
p1|p2
dev.off()

### get the pca plot of protein, rna is just the same ###

library(data.table)
library(ggplot2)
library(Rtsne)
library(umap)
library(factoextra)

setwd("D:\\2023Autumn\\genomics\\project")
dir.create("./results")

protein_data <-read.csv(file="Proteins-processed.csv",sep = ",",header = T)
row.names(protein_data) <- protein_data[,1]
protein_data <- protein_data[,-1]

cell_annotation <- read.csv(file="Cells.csv",sep = ",",header = T)
row.names(cell_annotation) <- cell_annotation[,1]
cell_annotation <- cell_annotation[,-1]
cell_annotation <- t(cell_annotation)
cell_annotation <- as.data.frame(cell_annotation, stringsAsFactors = FALSE)

columns_to_remove <- grep("protein", names(protein_data), ignore.case = TRUE)
protein_data <- protein_data[ , -columns_to_remove]
protein_data <- t(protein_data)
pca_result <- prcomp(protein_data, center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)
pca_data <- cbind(pca_data, celltype = cell_annotation$celltype)

pca_data$celltype[pca_data$celltype == "sc_m0"] <- "Macrophage"
pca_data$celltype[pca_data$celltype == "sc_u"] <- "U937"


ggplot(pca_data, aes(x = PC1, y = PC2, color = celltype)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 16)) +
  labs(title = "PCA of Protein Expression Data")

ggsave("./results/PCA_plot.png", plot = last_plot(), width = 10, height = 8, dpi = 300)

tsne_result <- Rtsne(as.matrix(protein_data), dims = 2, perplexity = 30, verbose = TRUE)
tsne_data <- as.data.frame(tsne_result$Y)
tsne_data$celltype <- cell_annotation$celltype

ggplot(tsne_data, aes(x = V1, y = V2, color = celltype)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 16)) +
  labs(title = "t-SNE of Protein Expression Data")

ggsave("./results/tsne_plot.png", plot = last_plot(), width = 10, height = 8, dpi = 300)


umap_result <- umap(protein_data)

umap_data <- as.data.frame(umap_result$layout)

umap_data$celltype <- cell_annotation$celltype

ggplot(umap_data, aes(x = V1, y = V2, color = celltype)) +
  geom_point() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 16)) +
  labs(title = "UMAP of Protein Expression Data")

ggsave("./results/umap_plot.png", plot = last_plot(), width = 10, height = 8, dpi = 300)



elbowplot <- fviz_eig(pca_result)

ggsave("./results/elbowplot.png", plot = last_plot(), width = 10, height = 8, dpi = 300)
library("data.table")
library(ggplot2)

dir.create("./image")
data <-read.csv(file="joint.csv",sep = ",",header = T)
row.names(data) <- data[,1]
data <- data[,-1]

calculateCorrelations <- function(data) {
  count_cells_protein <- sum(substr(names(data), 1, 1) == "i")
  count_cells_RNA <- ncol(data) - count_cells_protein
  protein <- data[, 1:count_cells_protein]
  rna <- data[, (count_cells_protein + 1):(count_cells_protein + count_cells_RNA)]
  t_protein <- t(protein)
  t_rna <- t(rna)
  correlation_matrix_pi <- cor(t_protein)
  correlation_matrix_ri <- cor(t_rna)
  n_cols <- ncol(correlation_matrix_pi)
  correlations <- numeric(n_cols)
  for (i in 1:n_cols) {
    correlations[i] <- cor(correlation_matrix_ri[, i], correlation_matrix_pi[, i])
  }
  hist(correlations)
  correlation_data <- data.frame(correlations = correlations)
  return(correlation_data)
}

calculateCorrelations_null <- function(data) {
  # Count the number of cells measured for protein (assumed to start with 'i')
  count_cells_protein <- sum(substr(names(data), 1, 1) == "i")
  count_cells_RNA <- ncol(data) - count_cells_protein
  protein <- data[, 1:count_cells_protein]
  rna <- data[, (count_cells_protein + 1):(count_cells_protein + count_cells_RNA)]
  t_protein <- t(protein)
  t_rna <- t(rna)
  correlation_matrix_pi <- cor(t_protein)
  correlation_matrix_ri <- cor(t_rna)
  n_cols <- ncol(correlation_matrix_pi)
  n_rows <- nrow(correlation_matrix_pi)
  correlation_matrix_ri <- correlation_matrix_ri[sample(n_rows), sample(n_cols)]
  correlation_matrix_pi <- correlation_matrix_pi[sample(n_rows), sample(n_cols)]
  correlations <- numeric(n_cols)
  for (i in 1:n_cols) {
    correlations[i] <- cor(correlation_matrix_pi[, i], correlation_matrix_ri[, i])
  }
  correlation_data <- data.frame(correlations = correlations)
  return(correlation_data)
}

correlation_data <- calculateCorrelations(data)
head(correlation_data)
correlation_data$gene_id <- row.names(data)
positive_control_genes <- correlation_data[correlation_data$correlations > 0.2, ]
negative_control_genes <- correlation_data[correlation_data$correlations < -0.2, ]

positive_control <- data[positive_control_genes$gene_id, ]

negative_control <- data[negative_control_genes$gene_id, ]

### Correlation statistics of positively and negatively regulated genes ###

positive_control_correlation <- calculateCorrelations(positive_control)

negative_control_correlation <- calculateCorrelations(negative_control)

#### NULL as control ###
null_correlation <- calculateCorrelations_null(data)
correlation_data <- subset(correlation_data, select = -gene_id)

positive_control_correlation$Group <- "Cluster I"
negative_control_correlation$Group <- "Cluster II"
null_correlation$Group <- 'Null'
correlation_data$Group <- 'All genes'

combined_data <- rbind(null_correlation,
                       correlation_data,
                       positive_control_correlation, 
                       negative_control_correlation)

combined_data$Group <- factor(combined_data$Group, levels = c("Null", "All genes", "Cluster I", "Cluster II"))

ggplot(combined_data, aes(x = Group, y = correlations, fill = Group)) + 
  geom_violin(trim = FALSE, alpha = 0.6) +  
  geom_segment(aes(x = 1.8, xend = 2.2, y = 0.2, yend = 0.2), linetype = "dashed", color = "blue") +  # Short line at y = 0.2
  geom_segment(aes(x = 1.8, xend = 2.2, y = -0.2, yend = -0.2), linetype = "dashed", color = "red") +  # Short line at y = -0.15
  annotate("text", x = 2, y = 0.25, label = "Cluster I", size = 5, color = "blue") +
  annotate("text", x = 2, y = -0.25, label = "Cluster II", size = 5, color = "red") +
  theme_minimal(base_size = 18) + 
  theme(text = element_text(size = 16),  
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5)) +  
  labs(title = "Violin Plot of Different Correlation Groups",
       x = "",
       y = "Correlation of r_i and p_i") +
  scale_y_continuous(breaks = seq(-0.4, 1, by = 0.2)) +  
  scale_fill_brewer(palette = "Set2") 

ggsave("./results/difference.png", plot = last_plot(), width = 12, height = 8, dpi = 300)

negative_control_normalized <- t(scale(t(negative_control)))

positive_control_normalized <- t(scale(t(positive_control)))
pheatmap(negative_control_normalized,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(255),
         breaks = seq(-0.5, 0.5, length.out = 256)) # Setting color breaks from -0.5 to 0.5


ggsave("./results/heatmap1.png", plot = last_plot(), width = 10, height = 8, dpi = 300)


pheatmap(positive_control_normalized,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(255),
         breaks = seq(-0.5, 0.5, length.out = 256)) # Setting color breaks from -0.5 to 0.5
ggsave("./results/heatmap2.png", plot = last_plot(), width = 10, height = 8, dpi = 300)



#====================== GO ========================
  positive_control_genes_sorted <- positive_control_genes[order(positive_control_genes$correlations, decreasing = TRUE), ]

negative_control_genes_sorted <- negative_control_genes[order(negative_control_genes$correlations, decreasing = FALSE), ]

top_genes_up <- head(positive_control_genes_sorted$gene_id, 10)
top_genes_down <- head(negative_control_genes_sorted$gene_id, 10)

write(top_genes_up, file = "topup10.txt")
write(top_genes_down, file = "topdown10.txt")


### CITE seq ###

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())


dir = c('/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample1/', 
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample2/',
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample3/',
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample4/',
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample7/',
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample8/',
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample9/',
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample10',
        '/Users/zhangguangxin/Desktop/Final Project/Cite-seq/GSE201048/Sample11')

names(dir) = c('Sample1','Sample2','Sample3','Sample4','Sample7','Sample8','Sample9','Sample10','Sample11') #多样本整合时命名barcode-ID    
counts <- Read10X(data.dir =dir)


scRNA <- CreateSeuratObject(counts[["Gene Expression"]])
adt <- CreateAssayObject(counts[["Antibody Capture"]])
all.equal(colnames(counts[["Gene Expression"]]), colnames(counts[["Antibody Capture"]]))
scRNA[["ADT"]] <- adt
Assays(scRNA)#[1] "RNA" "ADT"
rownames(scRNA[["ADT"]])#what protein do I have? 

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
scRNA <- subset(scRNA, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)

scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA,nfeatures = 2000)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA,dims = 1:20)

scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.9, verbose = FALSE)
scRNA <- RunUMAP(scRNA, dims = 1:20)
DimPlot(scRNA, label = TRUE)

#Finding differential genes
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scRNA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# Normalize ADT data
DefaultAssay(scRNA) <- "ADT"
scRNA <- NormalizeData(scRNA, normalization.method = "CLR", margin = 2)

#visualize
rownames(scRNA[["ADT"]])
p1 <- FeaturePlot(scRNA, "CD45-TotalSeqB", cols = c("lightgrey", "darkgreen")) + ggtitle("ADT")
DefaultAssay(scRNA) <- "RNA"
p2 <- DimPlot(scRNA, label = TRUE)+ggtitle("RNA")
p1|p2

#Section_5: Integrating
#The following statement sets the parallel speedup 
library(future)
plan("multiprocess", workers =4)
options(future.globals.maxSize = 2000 * 1024^2)

scRNA_WNN <- CreateSeuratObject(counts[["Gene Expression"]])
adt_WNN <- CreateAssayObject(counts[["Antibody Capture"]])
# add this assay to the previously created Seurat object
scRNA_WNN[["ADT"]] <- adt_WNN

DefaultAssay(scRNA_WNN) <- 'RNA'
scRNA_WNN[["percent.mt"]] <- PercentageFeatureSet(scRNA_WNN, pattern = "^MT-")
scRNA_WNN <- subset(scRNA_WNN, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)
scRNA_WNN <- NormalizeData(scRNA_WNN) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(scRNA_WNN) <- 'ADT'
VariableFeatures(scRNA_WNN) <- rownames(scRNA_WNN[["ADT"]])
scRNA_WNN <- NormalizeData(scRNA_WNN, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

scRNA_WNN <- FindMultiModalNeighbors(
  scRNA_WNN, reduction.list = list("pca", "apca"),
  dims.list = list(1:20, 1:10), modality.weight.name = "RNA.weight"
)
# UMAP visualization of data created based on weighted combination of RNA and pro
scRNA_WNN <- RunUMAP(scRNA_WNN, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
scRNA_WNN <- FindClusters(scRNA_WNN, graph.name = "wsnn", algorithm = 1, resolution = 0.6, verbose = FALSE)
p1 <- DimPlot(scRNA_WNN, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5)

# compare the graphs of cells formed by mRNA and WNN alone
p1 <- DimPlot(scRNA_WNN, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5)+ggtitle("WNN")
p2 <- DimPlot(scRNA,label = TRUE) + ggtitle("RNA")
p1|p2


### cite seq, the same process as data from SCoPE2###


library("data.table")
library(ggplot2)
data <-read.csv(file="Macro_Joint_LogNorm.csv",sep = ",",header = T)
row.names(data) <- data[,1]
data <- data[,-1]
data[is.na(data)] <- 0

calculateCorrelations <- function(data) {
  count_cells_protein <- 1329
  count_cells_RNA <- 1329
  rna <- data[, 1:count_cells_protein]
  protein <- data[, (count_cells_protein + 1):(count_cells_protein + count_cells_RNA)]
  t_protein <- t(protein)
  t_rna <- t(rna)
  correlation_matrix_pi <- cor(t_protein)
  correlation_matrix_ri <- cor(t_rna)
  n_cols <- ncol(correlation_matrix_pi)
  correlations <- numeric(n_cols)
  
  for (i in 1:n_cols) {
    correlations[i] <- cor(correlation_matrix_ri[, i], correlation_matrix_pi[, i])
  }
  hist(correlations)
  
  correlation_data <- data.frame(correlations = correlations)
  
  return(correlation_data)
}

calculateCorrelations_null <- function(data) {
  count_cells_protein <- 1329
  count_cells_RNA <- 1329
  protein <- data[, 1:count_cells_protein]
  rna <- data[, (count_cells_protein + 1):(count_cells_protein + count_cells_RNA)]
  t_protein <- t(protein)
  t_rna <- t(rna)
  correlation_matrix_pi <- cor(t_protein)
  correlation_matrix_ri <- cor(t_rna)
  n_cols <- ncol(correlation_matrix_pi)
  n_rows <- nrow(correlation_matrix_pi)
  correlation_matrix_ri <- correlation_matrix_ri[sample(n_rows), sample(n_cols)]
  correlation_matrix_pi <- correlation_matrix_pi[sample(n_rows), sample(n_cols)]
  correlations <- numeric(n_cols)
  for (i in 1:n_cols) {
    correlations[i] <- cor(correlation_matrix_pi[, i], correlation_matrix_ri[, i])
  }
  hist(correlations)
  correlation_data <- data.frame(correlations = correlations)
  return(correlation_data)
}

correlation_data <- calculateCorrelations(data)
null_correlation <- calculateCorrelations_null(data)

# Combine the data frames and create a new column to distinguish them
correlation_data$dataset <- "membrane gene"
null_correlation$dataset <- "Null"
combined_data <- rbind(correlation_data, null_correlation)

# Plotting the combined histogram
ggplot(combined_data, aes(x = correlations, fill = dataset)) +
  geom_histogram(position = "identity", alpha = 0.4, binwidth = 0.1, color = "black") +
  labs(title = "Correlations of RNA and Protein", x = "Correlations", y = "Frequency") +
  scale_fill_manual(values = c("blue", "red")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)) 
library(pheatmap)

correlation_data <- correlation_data[, "correlations", drop = FALSE]
numeric_data <- scale(correlation_data)

pheatmap(numeric_data,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(255),
         breaks = seq(-0.5, 0.5, length.out = 256))



### functions ###

sc.run = function(data){
  sc.runs = intersect(unique(data$Raw.file[data$celltype %in% cellnames]),names(data))
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
  return(sc.runs)
}

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