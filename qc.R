setwd("D://LILab//SC蛋白组//231122//final_result")
getwd()
dir.create(file.path(getwd(),"analysis"),recursive = F)
dir.create(file.path(getwd(),"figures"),recursive = F)
pipe = 'scp'
source("D://LILab//SC蛋白组//231122//code//packages.R")
source("D://LILab//SC蛋白组//231122//code//scp_functions.R")
source("D://LILab//SC蛋白组//231122//code//plot_functions.R")
cellnames = c("naiveT","actT") 
TMTlabel = factor(c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N"),
                  levels = c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N"),ordered = T)

#####load evidence_updated data and generate matadata(****meta need to be arranged by RI1-16*****)
evidence = fread(paste0(getwd(),"/ev_updated1.txt"),nThread = 10,header = T)
colnames(evidence)=gsub(" ",".",colnames(evidence))
colnames(evidence) = gsub("Reporter.intensity.corrected.","RI",colnames(evidence))

evidence$Raw.file = gsub("MS20230366_SCP_","",evidence$Raw.file)
evidence$Raw.file = gsub("_20231120235906","",evidence$Raw.file)
evidence$Raw.file = gsub("_20231120180441","",evidence$Raw.file)
evidence$Raw.file = gsub("_20231120194430","",evidence$Raw.file)
evidence$Raw.file = gsub("_20231120162451","",evidence$Raw.file)


sc_batch = unique(evidence$Raw.file) %>% grep(.,pattern ="B",value = T)
evidence = evidence %>% dplyr::filter( Raw.file %in% sc_batch )

nBatchs = sc_batch %>% length
metadata = data.frame( Raw.file = rep(sc_batch,each = 16),
                       channel = rep(factor(paste0("RI",1:16)),nBatchs),
                       TMTlabel = rep(TMTlabel,nBatchs)) %>% arrange(Raw.file,TMTlabel)

metadata$celltype = ifelse(metadata$Raw.file %in% grep(metadata$Raw.file,pattern = '[B|A]\\d+',value=T),
                           ifelse(as.numeric(str_extract(metadata$Raw.file,pattern = '(?<=[B|A])\\d+'))>=13,
                                  c("Carrier","Reference","Unused","Unused",rep(c("naiveT","actT"),5),"naiveT","Blank"),
                                  c("Carrier","Reference","Unused","Unused","Blank",rep(c("actT","naiveT"),5),"actT")),
                           c("Carrier","Reference","Unused","Unused",rep(c("actT","naiveT"),6)))
save(evidence,metadata,file = "analysis/Step0_raw_meta.RData",compress = T)


######scp analysis#######
pif = 0.8
pvalue = 0.01
meanscr = 0.1
medianCV = 0.4

gradient_cols = ggsci::pal_npg(palette = c("nrc"), alpha = 1)(6)
names(gradient_cols) = c("Carrier","Reference","Unused","Blank",cellnames)
show_col(gradient_cols)

Step1_generate_raw(ev_data = evidence,meta_data = metadata,
                   batchCol = "Raw.file",channelCol = "channel",
                   filt_with_min_psm = T,psm_threshold = 100,
                   samplePattern = "naiveT|actT",sc_batch = sc_batch)
load("analysis/Step1_raw_scp.RData")

### observ relativ intensity among single cell channels to see if any batches are abnormal

p = RI_plot(scp,assay = sc_batch[1],relative_to_Carrier = T)
for(i in c(sc_batch[-1])){
  p_tmp = RI_plot(scp,assay = i,relative_to_Carrier = T)
  p = p + p_tmp
  rm(p_tmp)
}

pdf("figures/P1_raw_RI_sc_batch.pdf",width = 7*2,height = 4*2)
p + plot_layout(guides = "collect",ncol = 2)
dev.off()

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
  # geom_histogram(binwidth = 0.1,bins = 30) +
  # scale_fill_manual( values = c("navy","red") ) +
  scale_x_log10(breaks = c(1e-3,1e-2,1e-1,1),limits = c(1e-4,1e0),expand = c(0,0)) +
  geom_vline( xintercept = meanscr,linewidth =1,lty = 2) +
  ylab("Batch")+ 
  Themes(pt_size = 14,titlefontsize = 20,
         axis_x_fontsize = 20,axis_y_fontsize = 20,
         fontsize =22 ,rotate = T,linesize = 1)  
pdf("figures/P2_MeanSCR.pdf",width = 6,height = 4)
p1 
dev.off()

scp = QFeatures::filterFeatures( scp, ~ !is.na(MeanSCR) & MeanSCR < meanscr )

Stats = rbind(Stats,data.frame(dims(scp)))
save(scp,Stats,file = "analysis/Step1-1_pep_filt_scp.RData",compress = T)
# load("analysis/Step1-1_pep_filt_scp.RData",verbose = T)

Step2_psm_to_pep(scp=scp,Raw.file_feature = "^\\w+",
                 if_generate_psm = F,
                 medianCV_nobs = 5,
                 medianCV_batches = sc_batch,
                 sc_batch = sc_batch) ### may take long time
# load("analysis/Step2_pep_aggregated_scp.RData",verbose = T)

medianRI = colMeans( assay(scp[["peptides"]]),na.rm = T )
scp$medianRI = medianRI

k = getWithColData(scp, "peptides") %>%
  colData %>%
  data.frame
k$celltype = factor(k$celltype,levels = rev(c("Carrier","Reference","Unused","Blank",cellnames)))
p2 = ggplot(k,aes( x = log10(medianRI), y = celltype, fill =  celltype )) +
  geom_boxplot(size = 0.5) + 
  # stat_summary(fun = mean,size = 4,color = "Black",
  #              geom = "point" ) + 
  # stat_summary(fun.data = data_summary,linewidth = 1,color = "Black",
  #              geom = "errorbar",width = 0 ) +  
  scale_fill_manual( values = gradient_cols ) + 
  scale_x_continuous(limits = c(-1,7)) +
  ggtitle("MedianRI") +#facet_grid(Batch~.,drop = T,scales = "free",space = "free") +
  ggtitle(NULL) + xlab("log10(medianRI)") + ylab(NULL) + 
  Themes(pt_size = 2,fontsize = 20,
         axis_x_fontsize = 18,axis_y_fontsize = 18,
         titlefontsize = 20,linesize = 1,rotate = F)+ theme(legend.position = "NULL")
pdf("figures/P3_MedianRI.pdf",width = 5,height = 4)
p2
dev.off()

a = summary(k$MedianCV[k$celltype %in% c("Carrier","Reference")])
a
p3 = k %>%
  ggplot(aes(x = MedianCV,y = celltype,fill = celltype,group = celltype)) +
  geom_density_ridges(size = 1,scale = 2,rel_min_height = 1e-3,bandwidth = 0.015) + 
  scale_fill_manual(values = alpha(gradient_cols,alpha = 0.6)) +
  geom_vline(xintercept = c(0.4),
             size = 1,linetype = 2) +
  xlab("MedianCV")+ ylab(NULL) + #facet_grid(Batch~.,drop = T,scales = "free",space = "free") +
  Themes(fontsize = 20, 
         axis_x_fontsize = 20,axis_y_fontsize = 20,
         titlefontsize = 18,linesize = 1,rotate = F) + theme(legend.position = "NULL") 
pdf("figures/P4_medianCV.pdf",width = 6,height = 4)
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
  tbl = cbind(tbl,dims(scp) %>% as.data.frame() %>%  dplyr::select(c(intersect(names(scp),unique(metadata$Raw.file)))))
}else{
  scp = scp[,!is.na(scp$MedianCV) & scp$MedianCV < medianCV & scp$celltype %in% cellnames, ]
  tbl = dims(scp) %>% as.data.frame() %>%  dplyr::select(c(intersect(names(scp),unique(metadata$Raw.file))))
}
colnames(tbl) = ifelse(colnames(tbl) %in% grep(colnames(tbl),pattern = '^\\d',value = T),paste0("X",colnames(tbl)),colnames(tbl))
Stats = rbind(Stats,tbl)

Stats$Steps = factor(rep(c("Raw","Contaminant",paste0("PIF",pif),paste0("q_val",pvalue),paste0("MeanSCR",meanscr),paste0("MedianCV",medianCV)),each = 2),
                     levels = c("Raw","Contaminant",paste0("PIF",pif),paste0("q_val",pvalue),paste0("MeanSCR",meanscr),paste0("MedianCV",medianCV)))
Stats$Data = rep(c("PSMs","Cells"),times = 6)

### For the PSM in each steps ####
tmp = Stats %>% 
  # filter( Data=="PSMs" ) %>% 
  gather( key = "Raw.file",value = "Number", -Steps,-Data )
tmp$Raw.file = gsub("X","",tmp$Raw.file)
tmp = tmp %>% dplyr::filter(Raw.file %in% c(sc_batch)) %>% 
  unique() 
tmp$Number[tmp$Data=="Cells"&tmp$Steps!="Median0.4"] = tmp$Number[tmp$Data=="Cells"&tmp$Steps!="Median0.4"]  - 5
tmp = tmp %>% dplyr::group_by(Steps,Data) %>% dplyr::mutate(median = median(Number,na.rm = T))

p4 = ggplot( tmp[tmp$Data=="PSMs",] ) + 
  geom_bar(aes(x = Steps, y = median, group = Steps,fill = Steps),width = 0.8,
           stat = "identity",position = "dodge") +
  # geom_boxplot(aes(x = Steps, y = Number, group = Steps,fill = Steps),size = 0.7,varwidth = 0.8,color ="black",notch = F) +
  # geom_density_ridges(size = 1,scale = 2,rel_min_height = 1e-3,bandwidth = 0.015) + 
  # scale_fill_manual( values = alpha(brewer.pal("YlOrRd",n=6)),1 ) + 
  scale_fill_manual( values = alpha(ggsci::pal_npg(palette = c('nrc'))(6),1 )) + 
  stat_summary(aes(x = Steps, y = Number, group = Steps),
               fun.data = data_summary,linewidth = 0.7,
               color = "black",fatten = 2,
               geom = "pointrange") +
  Themes(pt_size = 8,fontsize = 18,linesize = 1,
         axis_x_fontsize = 16,axis_y_fontsize = 18,
         titlefontsize = 22,rotate = T) + 
  xlab(NULL) + ylab( "Number" ) + facet_grid(Data~.,scales = "free_y",space = "free_x") + theme(legend.position = "none")

pdf("figures/P5_statistics.pdf",width = 6,height = 4)
p4
dev.off()
save(scp,Stats,file = "analysis/Step2-1_cell_filt_scp.RData",compress = T)

Step3_normalization_pep_aggre_protein(scp = scp, pep_pNA = 0.99,protein_pNA = 0.99,filt_pep_NA = F,filt_pro_NA = T)
plot(scp,interactive = T)
# load("analysis/Step3_protein_aggregated_normed_scp.RData")
scp[["proteins_norm"]] %>%
  assay %>% 
  is.na %>%  mean() #[1]0.2938052

scp = impute(scp,
             i = "proteins_norm",
             method = "knn",
             k = 5, rowmax = 1, colmax= 1,
             maxp = Inf, rng.seed = 1234 )

scp[["imputedAssay"]] %>%
  assay %>%
  is.na %>%
  mean #### 0
save(scp,file = "analysis/Step3-1_imputed_scp.RData",compress = T)

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
save(scp,file = "analysis/Step3-2_batch_correct_scp.RData",compress = T)
#load("analysis/Step3-2_batch_correct_scp.RData")

#### correlate with Nikolai dataset
a = getWithColData(scp,"proteins_batchC") #proteins_norm
# 
# hl = a[rownames(a) %in% c(pros),a$celltype=="actT"]
# pearscor = cor(assay(hl),method = "pearson",use = "na.or.complete")
pearscor = cor(assay(a),method = "pearson",use = "na.or.complete")
summary(pearscor)
quantile(pearscor, probs = c(0.05, 0.95), na.rm = T, names = TRUE)
pearscor[pearscor>0.5] = 0.5
pearscor[pearscor<(-0.5)] = -0.5

cell_cols = alpha(c("#4d82e9","#d65341"),0.8)
names(cell_cols) =cellnames
cluster_cols = alpha(brewer.pal('Paired',n=3),0.9)
names(cluster_cols) = unique(a$Raw.file)
top = HeatmapAnnotation(Batch = a$Raw.file,
                        Cell = a$celltype,
                        col = list(Batch = cluster_cols,
                                   Cell = cell_cols),
                        # annotation_label = c("Cell","Batch"),
                        show_legend = T,
                        show_annotation_name = F,
                        annotation_legend_param = list(labels_gp  = gpar(fontsize = 15),
                                                       title_gp = gpar(fontsize = 15),
                                                       legend_direction = "horizontal")
)
pdf("figures/P6_cell_corr_heatmap.pdf",width = 10,height = 8)
Heatmap( pearscor, 
         heatmap_legend_param = list( title = "Pearson",
                                      at = c(-0.5,-0.25,0,0.25,0.5),
                                      labels = c(-0.5,-0.25,0,0.25,0.5),
                                      labels_gp  = gpar(fontsize = 15),
                                      direction = "vertical",
                                      title_position = "topcenter",
                                      title_gp = gpar(fontsize = 15)),
         show_heatmap_legend = T,
         # col = viridis_pal(option = "magma",begin = 1,end = 0)(20),
         col = circlize::colorRamp2(breaks = seq(-0.5,0.5,length.out = 9),
                                    colors = colorRampPalette(c("blue","white","red"))(9),#colorRampPalette(c("#00aeff","Black","#f9d31a"))(9),
                                    transparency = 0),
         cluster_columns = T,
         # row_km = 4,
         cluster_rows = T ,
         show_column_names = F,
         clustering_method_rows = "complete",#ward.D2
         clustering_distance_rows = "euclidean",
         #("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
         show_row_names = F,
         top_annotation = top,
         # left_annotation = left,
         column_names_rot = 45,
         width = unit(3,"in"),height = unit(3,"in"),
         column_names_gp = gpar(fontsize = 15) )
dev.off()

##### protein number per cell
tmp = getWithColData(scp,"proteins_norm")
pro_valid = apply( assay(tmp),2, function(x){
  Valid = length(x[!is.na(x)])
}) %>% as.data.frame() 

colnames(pro_valid) = c("pro_num")
pro_valid = cbind(pro_valid,scp@colData)
p5 = ggplot(pro_valid) + 
  # geom_boxplot(aes(x = Batch,y = pro_num, fill = celltype), color = "Black",size = 1) + 
  # geom_violin(aes(x = Batch,y = pro_num, fill = celltype), color = "Black",linewidth = 1) +
  geom_density_ridges(aes(y = Raw.file,x = pro_num, fill = celltype,color = celltype),
                      jittered_points = TRUE, scale = 1, rel_min_height = .01,size = 0.5,
                      point_shape = "|", point_size = 4, 
                      position = position_points_jitter(height = 0),
                      # position = "points_sina",stat = "density_ridges",
                      bandwidth = 3) +
  scale_y_discrete(expand = c(0.1,0)) +
  # geom_text(aes(x = PC1,y = PC2,label = Raw.file),size = 6,check_overlap = T)+
  # scale_color_manual(values = alpha(c("#6097b0","#ed7057"),0.6)) +
  scale_fill_manual(values = alpha(c("#6097b0","#ed7057"),0.5)) +
  scale_color_manual(values = alpha(c("#6097b0","#ed7057"),1)) +
  ylab("Batch") + xlab( "Protein Number" ) + #facet_grid(Set~.,scales = "free",space = "free") + theme(legend.position = "none")+
  Themes(pt_size = 4,fontsize = 25,
         axis_x_fontsize = 26,axis_y_fontsize = 26,
         linesize = 1.5,rotate = F)
p5
pdf("figures/P7_protein_number_new.pdf",width = 8,height = 3)
p5
dev.off()

names(pro_density) = pro_density = unique(scp$Raw.file)
pro_density = lapply( as.list(pro_density), function(x){
  p = assay(tmp[,tmp$Raw.file==x]) %>% drop_Nas(.,dims = 1,partitial = F)
  Valid_density = apply(p,1,function(i){
    return(length(i[!is.na(i)])/length(i))
  }) %>% as.data.frame() %>% rownames_to_column( .,"Proteins" ) 
  colnames(Valid_density) = c("Proteins","pro_density")
  return(Valid_density)
}) %>% plyr::ldply(.,.id = "Raw.file")

# a_pros = unique(pro_density$Proteins[pro_density$Batch=="A"])
# b_pros = unique(pro_density$Proteins[pro_density$Batch=="B"])
# pub_pros = unique(pro_density$Proteins[pro_density$Batch=="Public"])
# pdf("figures/P8_protein_overlap.pdf",width = 5,height = 5)
# library(eulerr)
# require(ggplotify)
# fit_up = euler(c("A"=length(setdiff(a_pros,c(pub_pros,b_pros))),
#                  "B"=length(setdiff(b_pros,c(pub_pros,a_pros))),
#                  "Public"=length(setdiff(pub_pros,c(b_pros,a_pros))),
#                  "A&B&Public"=length(intersect(a_pros,intersect(pub_pros,b_pros))),
#                  "A&B"=length(setdiff(intersect(a_pros,b_pros),pub_pros)),
#                  "B&Public"=length(setdiff(intersect(pub_pros,b_pros),a_pros)),
#                  "A&Public"=length(setdiff(intersect(pub_pros,a_pros),b_pros))
# ))
# plot(fit_up,quantities = list(cex = 1.8,hjust = 1,fontfamily = "sans"),# TRUE
#      fill=alpha(brewer.pal("Set3",n=3),alpha = 0.5),#"transparent",
#      edges = F,
#      cat.fontface = 8,
#      main = list(label="Proteins",cex = 2,fontfamily = "sans",just = "center"),
#      shape = "ellipse",fontfamily = "sans",
#      legend = F,#list(cex = 2,side = "right"),
#      label = list(fontfamily = "sans",cex = 2,check.overlap = T,hjust = 0),
#      adjust_labels = F)  %>% as.ggplot()
# dev.off()

# pros = intersect(a_pros,intersect(pub_pros,b_pros))
# a = getWithColData(scp,"proteins_norm")
# a = a[rownames(a) %in% c(pros),]
# 
# protein_expr = assay(a) %>% drop_Nas(.,dims=1,partitial = F)
# protein_expr = as.data.frame(protein_expr) %>% rownames_to_column(.,"Proteins") %>% gather(.,key = "Cells",value = "expr",-Proteins)
# n = as.data.frame(a@colData) %>% rownames_to_column(.,"Cells") %>% select(Raw.file,Cells,celltype,Batch)
# # n$Batchnew = ifelse(n$Batch=="Public","Public","OurStudy")
# protein_expr = merge(n,protein_expr,by = "Cells") %>% 
#   dplyr::group_by(Batch,Proteins,celltype) %>% 
#   dplyr::mutate( MedianExpr = median(expr,na.rm = T) ) %>% 
#   select(-Cells,-expr,-Raw.file,-Batch) %>% unique()
# 
# a_cells = n$Cells[n$Batchnew=="A"]
# pub_cells = n$Cells[n$Batchnew=="Public"]
# 
# k = as.data.frame(assay(hl)[,c("wAP329RI9",'wAP332RI12')]) %>% rownames_to_column(.,"Proteins")

# library(GGally)
# protein_expr$id = paste0(protein_expr$Batch,"_",protein_expr$celltype)
# k = protein_expr[,-c(1,2)] %>% spread( .,key = id,value = MedianExpr) 
# pdf("figures/cor_point.pdf",width = 12,height = 12)
# limitRange <- function(data, mapping, ...){ 
#   ggplot(data = data, mapping = mapping, ...) + 
#     geom_point(...) + 
#     geom_smooth(method = "lm", se = FALSE,color = alpha("grey",0.9),linewidth = 1) +
#     scale_y_continuous(limits = c(-2.1, 3)) +
#     scale_x_continuous(limits = c(-2.1, 3)) 
# }
# ggpairs(k,columns = 2:7,lower = list(continuous = wrap(limitRange, size = 0.8, color = alpha("red",0.5))),
#         upper = list(continuous = wrap(ggally_cor, size = 6, color = "black",use = "na.or.complete")) ) + 
#   theme_cowplot() + Themes(axis_x_fontsize = 18,linesize = 0.8,axis_y_fontsize = 18,fontsize = 6,titlefontsize = 18)
# dev.off()
# 
# pdf("figures/cor_point2.pdf",width = 12,height = 12)
# limitRange2 <- function(data, mapping, ...){ 
#   ggplot(data = data, mapping = mapping, ...) + 
#     geom_point(...) + 
#     geom_abline(intercept = 0,slope = 1,linewidth=1,color = alpha("red",0.8),linetype = 2) +
#     scale_y_continuous(limits = c(-2.1, 3)) +
#     scale_x_continuous(limits = c(-2.1, 3)) 
# }
# ggpairs(k,columns = 2:7,lower = list(continuous = wrap(limitRange2, size = 1, color = alpha("black",0.4))),
#         upper = NULL ) + 
#   theme_cowplot() + Themes(axis_x_fontsize = 18,linesize = 0.8,axis_y_fontsize = 18,fontsize = 6,titlefontsize = 18)
# dev.off()

library(ComplexHeatmap)
tmp = getWithColData(scp,"proteins")

protein_detect =  assay(tmp) #%>%
# drop_Nas(.,dims = 1,partitial = F)
protein_detect[!is.na(protein_detect)] = 1
protein_detect[is.na(protein_detect)] = 0

protein_detect = as.data.frame(protein_detect) %>% rownames_to_column(.,"Proteins") %>% gather(.,key = "Cells",value = "is_detect",-Proteins)
n = as.data.frame(tmp@colData) %>% rownames_to_column(.,"Cells") %>% select(Raw.file,Cells,celltype)
protein_detect = merge(n,protein_detect,by = "Cells") %>% 
  dplyr::group_by(Raw.file,Proteins,celltype) %>% 
  mutate( Detect_density = sum(is_detect)/length(is_detect) ) %>% 
  select(-Cells,-is_detect) %>% unique()
k = protein_detect %>%
  dplyr::mutate( tmp = paste0(Raw.file,"_",celltype) ) 
anno = k[,c(1:2,5)] %>% unique
k = k[,-c(1:2)] %>% spread( .,key = tmp,value = Detect_density) %>% column_to_rownames(.,"Proteins")

top = HeatmapAnnotation(Batch = anno$Raw.file,
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
         width = unit(2,"in"),height = unit(4,"in"),
         column_names_gp = gpar(fontsize = 15) )
dev.off()
a = protein_detect %>% spread( .,key = celltype,value = Detect_density) 
a$Protrein_type = ifelse( (a$naiveT+a$actT)==2,"All detect",
                          ifelse( (a$naiveT+a$actT)==0,"Non detect",
                                  ifelse( a$naiveT>=0.5&a$actT<=0.4,"naiveT specific",
                                          ifelse( a$actT>=0.5&a$naiveT<=0.4,"actT specific","Variable"))))
a$Protrein_type = factor(a$Protrein_type,levels = rev(c("All detect","Variable","naiveT specific","actT specific","Non detect")),ordered = T)
a = a[a$Protrein_type!="Non detect",]
pdf("figures/P9_pro_detect_bar.pdf",width = 6,height = 4)
ggplot(a) +
  geom_bar( aes(x = Raw.file,fill = Protrein_type ),stat = "count",position = position_stack(),width = 0.9 ) + 
  scale_fill_manual(values = c("#d65240","#54a056","#f1be44","skyblue")) + 
  xlab(NULL) + ylab('Protein Number') + 
  annotate(geom = "segment",color = gradient_cols[1],x=0.55,xend = 1.45,y = -10,yend = -10,size = 3)+
  annotate(geom = "segment",color = gradient_cols[3],x=1.55,xend = 2.45,y = -10,yend = -10,size = 3)+
  annotate(geom = "segment",color = gradient_cols[4],x=2.55,xend = 3.45,y = -10,yend = -10,size = 3)+
  # annotate(geom = "segment",color = gradient_cols[5],x=11.55,xend = 14.45,y = -10,yend = -10,size = 3)+
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
  geom_point(aes(color = celltype),size =5) + 
  # geom_text(aes(x = PC1,y = PC2,label = Raw.file),size = 6,check_overlap = T)+
  scale_shape_manual(values = c(15:21,8,13)) +
  scale_color_manual(values = cell_cols) +
  xlab(paste0("PC1 (",round(deviation[1],digits = 1),"%)")) +
  ylab(paste0("PC2 (",round(deviation[2],digits = 1),"%)")) +
  Themes(pt_size = 4,fontsize = 20,
         axis_x_fontsize = 25,axis_y_fontsize = 25,
         linesize = 1.5,rotate = F)
p8
pdf("figures/P10_raw_PCA.pdf",width = 7,height = 5)
p8
dev.off()


#### general pca reduction
mat = assay(tmp)
pearson = cor(t(mat))
# pearson = cor(t(mat %>% t %>% scale(.,center = T,scale = T) %>% t ))
head(pearson[1:10,1:10])
rsum = rowSums(pearson^2)

# Calculate the weighted data matrix:
mat = diag(rsum) %*%  mat
pca.imp.cor <- cor(mat)

# pca = prcomp(t(mat))
pca = prcomp(pca.imp.cor)


# pca = prcomp( t(assay(scp.tmp) %>% t %>% scale(.,center = T,scale = T) %>% t ))
deviation = (pca$sdev)^2/sum((pca$sdev)^2)*100
sum(deviation[1:2]) ;length(deviation)
# [1] 21.45363
# [1] 34
anno = tmp@colData %>% as.data.frame()
pca_mat = cbind(pca$x,anno)

# pearson = cor(scp.tmp)
pca_mat$celltype = factor(pca_mat$celltype,levels = cellnames)
pca_mat$Raw.file = as.factor(pca_mat$Raw.file)
p9 = ggplot(pca_mat,aes(x = PC1,y = PC2)) + 
  geom_point(aes(color = celltype,shape = Raw.file),size =5) + 
  # geom_text(aes(x = PC1,y = PC2,label = Raw.file),size = 6,check_overlap = T)+
  scale_shape_manual(values = c(15:21,8,13)) +
  scale_color_manual(values = cell_cols) +
  xlab(paste0("PC1 (",round(deviation[1],digits = 1),"%)")) +
  ylab(paste0("PC2 (",round(deviation[2],digits = 1),"%)")) +
  Themes(pt_size = 4,fontsize = 20,
         axis_x_fontsize = 25,axis_y_fontsize = 25,
         linesize = 1.5,rotate = F)
p9
pdf("figures/P11_weighted_PCA_new.pdf",width = 7,height = 5)
p9
dev.off()


#### only single cells
{
  normedassay0 = getWithColData(scp, "proteins_norm")
  # assay(normedassay0,"assay") = assay(normedassay1,"assay") %>% t() %>% apply( .,2,function(x){x = (x-mean(x,na.rm = T))/sd(x,na.rm = T)} ) %>% t()
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
  
  # pearson = cor(t(assay(normedassay3,"assay")))
  # rsum = rowSums(pearson^2)
  # assay(normedassay3,"assay",withDimnames = F) = diag(rsum) %*%  assay(normedassay3,"assay") 
  # add = as.Seurat( normedassay3,
  #                  counts = "aggcounts", data = "assay" )
  # ttt = GetAssayData( object = add,slot = "data" )
  # seur[["weighted_imputed_batchC"]] = CreateAssayObject( data = ttt )
  # 
  # tmp = assay(normedassay3,"assay") %>% cor()
  # tmp = SingleCellExperiment(assays = list(assay = tmp,aggcounts = tmp))
  # tmp@colData = normedassay3@colData
  # 
  # add = as.Seurat( tmp,counts = "aggcounts", data = "assay" )
  # ttt = GetAssayData( object = add,slot = "data" )
  # seur[["weighted_cor_imputed_batchC"]] = CreateAssayObject( data = ttt )
  
  #   load("mergeGO_median_mat.RData",verbose = T)
  #   merge = column_to_rownames(merge,"GOterms")
  #   add = CreateSeuratObject( merge,assay = "GOterm",meta.data = seur@meta.data )
  #   ttt = GetAssayData( object = add,slot = "data" )
  #   seur[["GOterm"]] =  CreateAssayObject( data = ttt )
  Idents(seur) = "celltype"
  save(seur,file = "seurat.RData",compress = T)
}


seur
# 
set.seed(1234)

VlnPlot( seur,features = c("nFeature_originalexp"),group.by = "Raw.file")

DefaultAssay(seur)  = "proteins_norm"
seur = FindVariableFeatures(seur,nfeatures = 2000,selection.method = "vst")
seur = ScaleData( seur,do.scale = T,do.center = F )
seur = Seurat::RunPCA( seur,assay = "proteins_norm",npcs = 60)
ElbowPlot(seur,ndims = 60)
DimPlot( seur,reduction = "pca",group.by = c("celltype","Batch"),
         # cols = alpha(c("#6097b0","#ed7057"),0.6),
         pt.size = 3) & 
  Themes(rotate = F) 

seur = RunUMAP( seur,dims = 1:30,reduction = "pca",n.neighbors = 30,min.dist = 0.05 )#n.neighbors = 30,min.dist = 0.05 
# seur = RunaiveTSNE( seur,dims = 1:20,reduction = "pca" )
DimPlot( seur,reduction = "umap",group.by = c("celltype","Batch"),pt.size = 2 )

# seur = subset(seur,subset = Batch!="H")
DefaultAssay(seur)  = "imputed_raw_Combat"
VariableFeatures(seur)=rownames(seur)
seur = ScaleData( seur,do.scale = T,do.center = F )
seur = Seurat::RunPCA( seur,assay = "imputed_raw_Combat",npcs = 30)
seur = Seurat::RunUMAP( seur,dims = 1:30,reduction = "pca",n.neighbors = 10,min.dist = 0.2 )
seur = Seurat::RunaiveTSNE( seur,dims = 1:30,reduction = "pca")

DimPlot( seur,reduction = "pca",group.by = c("celltype","Batch"),
         # cols = alpha(c("#6097b0","#ed7057"),0.6),
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
# seur = FindNeighbors(seur, dims = 1:20)
# seur = FindClusters(seur, resolution = 0.8)
Idents(seur) = "celltype"
table(seur$seurat_clusters,seur$celltype)

DefaultAssay(seur)  = 'imputed_raw_Combat'
Idents(seur) = "celltype"
seur = ScaleData(seur,do.scale = T,do.center = F)
markers = FindAllMarkers(seur,logfc.threshold = 0,min.pct = 0,only.pos = T,return.thresh = 0.05,slot = "data" )
markers = markers[markers$p_val<0.05,]
write.table(rownames(seur),"proteins.txt",quote = F,col.names = F,row.names = F)
library(ensembldb)
# ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
# genes = ensembldb::proteins(ensdb,columns = c("tx_id","uniprot_id","protein_id","uniprot_db"))
# genes 
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
go_hl = enrichGO(marker_gene$`Gene Names (primary)`[marker_gene$cluster=="naiveT"],
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",ont = "BP",minGSSize = 5,
                 pAdjustMethod = "BH",universe =database(species)[["wg.gene"]],
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
go_hl_simplified = clusterProfiler::simplify(x = go_hl ,cutoff=0.4,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值。
go_u9 = enrichGO(marker_gene$`Gene Names (primary)`[marker_gene$cluster=="actT"],
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",ont = "BP",minGSSize = 5,
                 pAdjustMethod = "BH",universe = database(species)[["wg.gene"]],
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
go_u9_simplified = clusterProfiler::simplify(x = go_u9 ,cutoff=0.6,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值。


p1 = GO_Plot(go_hl,fontsize = 7,fill_color = alpha("#6097b0",0.6),
             padj_threshold = 0.3,show_num = 15,keywords = c("cytokine" ),discard = "DNA")
p1
p2 = GO_Plot(go_u9,fontsize = 7,fill_color = alpha("#ed7057",0.6),
             padj_threshold = 0.3,show_num = 15,keywords = c("leu","immun","cytok" ),discard = "DNA")
pdf("figures/P12_GO.pdf",width = 15,height = 8)
p1|p2
dev.off()


library(ComplexHeatmap)
cell_cols = alpha(c("#6097b0","#ed7057"),0.6)
names(cell_cols) =cellnames

cluster_cols = alpha(c("#DF4E50","#BAacd1","#EDE5a1"),0.7)
names(cluster_cols) = unique(seur$Batch)
# seur = ScaleData(seur,do.center = F,do.scale = T)
k = seur[["imputed_raw"]]@data[markers$gene,] %>% t %>% scale(.) %>% t()
top = HeatmapAnnotation(Cells = seur@meta.data$celltype,
                        Batchs = seur@meta.data$Batch,
                        col = list(Cells = cell_cols,
                                   Batchs = cluster_cols),
                        show_legend = T,
                        show_annotation_name = F,
                        annotation_legend_param = list(labels_gp  = gpar(fontsize = 15),
                                                       title_gp = gpar(fontsize = 15))
)
summary(k[,1:10])
quantile(k,0.95,na.rm = T)
quantile(k,0.05,na.rm = T)
k[k>1.5]=1.5
k[k<(-1.5)] = -1.5

pdf("figures/heatmap.pdf",width = 8,height = 4)
Heatmap( k, 
         heatmap_legend_param = list( title = "Expr",
                                      at = seq(-1.5,1.5,by = 0.75),
                                      labels =seq(-1.5,1.5,by = 0.75),
                                      labels_gp  = gpar(fontsize = 15),
                                      direction = "vertical",
                                      title_gp = gpar(fontsize = 15)),
         show_heatmap_legend = T,
         col = colorRampPalette(c("#00aeff","Black","#f9d31a"))(9),
         # col = colorRampPalette(brewer.pal("RdYlBu",n=9))(20),
         cluster_columns = T,row_km = 2,
         cluster_rows = T ,
         show_column_names = F,
         clustering_method_rows = "complete",#ward.D2
         clustering_distance_rows = "euclidean",
         #("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
         show_row_names = F,
         top_annotation = top,
         # left_annotation = left,
         column_names_rot = 45,
         column_names_gp = gpar(fontsize = 15) )
dev.off()
