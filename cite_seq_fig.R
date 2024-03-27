
library(Seurat)
library(RColorBrewer)
library(gridExtra)

library(ggplot2)
library(Rmagic)
library(phateR)

library(readr)
library(viridis)
library(visNetwork)


###############################################
################################ Figure 1 Code  ###########################
#setwd("~/NatNeu2022/test_fig")
imm_all<-readRDS("./data.rds")

## Figure 1a
cols<-c(brewer.pal(12, "Set3"),brewer.pal(n = 8, name = "Dark2"),"#FDB462",
        "#0000CD","#7B68EE","#4682B4","#808000","#7CFC00") 
f1a<-DimPlot(imm_all,label=TRUE,reduction = "tsne",cols= cols)+NoLegend()

#以下操作可以实现从Clusters中抓取特定的一类或几类细胞;
#归一化数据
imm_all <- NormalizeData(imm_all, normalization.method = "LogNormalize", scale.factor = 10000)
object_17 <- imm_all[,imm_all$seurat_clusters %in% c(17)]
##构建稀疏矩阵
sparse <- Matrix(object_17@assays$ADT@counts,sparse = T)
sparse <- NormalizeData(sparse, normalization.method = "LogNormalize", scale.factor = 10000)
##输出基因feature
features <- data.frame(ID = sparse@Dimnames[[1]],Name = sparse@Dimnames[[1]],EX = "Antibody capture")
write.table(x = features,file = "/Users/houshuyang/Desktop/features.tsv",sep = "\t",quote = F,col.names = F,row.names = F)
##输出barcode
write(x = sparse@Dimnames[[2]],file = "/Users/houshuyang/Desktop/barcodes.tsv")
##输出表达矩阵
writeMM(sparse,file = "/Users/houshuyang/Desktop/matrix.mtx")

library(Matrix)
mtx_data <- readMM("./Macrophage_RNA/matrix.mtx")
writeMM(pro,file = "./matrix.csv")

## figure 1b
plt2<-c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080")
features<- c("CD45-Ab","CD3-Ab","CD19-Ab","CD20-Ab","CD25-Ab","CD27-Ab",
             "CD14-Ab","CD11b-Ab","CD56-Ab","CD16-Ab","CD69-Ab","HLA-DR-Ab")
p<-FeaturePlot(imm_all, features = features,reduction = "tsne", 
            min.cutoff = "q05", max.cutoff = "q95", ncol = 6,cols = plt2,combine = FALSE)



for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}


### Figure 1c
clus_freq<-data.frame(clus=Idents(imm_all),sname=imm_all$sampleName)

head(clus_freq)
table(clus_freq$clus)
ctype<-as.character(clus_freq$clus)
ctype[grep("15|8|16|17|19|21",ctype)]<-"Imm"
ctype[grep("13|20|22",ctype)]<-"NVU"
ctype[grep("18",ctype)]<-"Oligo"
ctype[grep("24|23|25",ctype)]<-"Oth"
ctype[!(ctype %in% c("Imm","NVU","Oligo","Oth"))]<-"Micr"


####
clus_freq$clus<-ctype
df<-as.data.frame(xtabs(~clus+sname,data=clus_freq))
df$sname<-factor(df$sname,levels = c("P6.C","P6.B","P6.A","P5.B", "P5.A","P4",
                                     "P3.B","P3.A","P2", "P1.B", "P1.A"))
cplt<-c("#E97171","#FBD46D","#4F8A8B","#222831","#AD9D9D","#005086","#318FB5","#F7D6BF","#B0CAC7")
f1c<-ggplot(df, aes(fill=clus, y=Freq, x=sname)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = cplt[1:5])+xlab("")+ylab("frequency")+
  theme_classic()+theme(axis.text=element_text(size=20,color = "black"))+coord_flip()

lay<-rbind(c(1,2),
           c(3,3))
figure1<-grid.arrange(grobs= list(f1a,f1c,f1b),layout_matrix=lay)
figure1

# setwd("~/NatNeu2022/test_fig")
# ggsave(filename="figure1.pdf", 
#        plot = figure1, 
#        device = cairo_pdf, 
#        width = 240, 
#        height = 270, 
#        units = "mm")

################################ Figure 1 Code ends ###########################

############################# Figure 2 ########################################
# read velmeshev et al data from seurat R object
velmeshev<-readRDS("./tsne_umap_rawMatrix_microglia_seurat_object_autism_data.rds")
f2a<-DimPlot(velmeshev, reduction = 'tsne',group.by = "diagnosis",label = FALSE)+NoAxes()

plt2<-c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080")
features <- c("P2RY12","CX3CR1","AIF1","CSF1R","IL18")
p<-FeaturePlot(velmeshev, features = features,reduction = "tsne", 
            min.cutoff = "q05", max.cutoff = "q95", ncol = 5,col=plt2,combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

f2b<-cowplot::plot_grid(plotlist = p,nrow = 1)


masuda<-readRDS("./Masuda_seurat_object_tsne_umap.rds")
plt2<-c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080")
#DimPlot(masuda, reduction = 'tsne')
features <- c("P2RY12","CX3CR1","AIF1","CSF1R","IL18")
p<-FeaturePlot(masuda, features = features,reduction = "tsne", 
               min.cutoff = "q05", max.cutoff = "q95", ncol = 5,col=plt2,pt.size = 0.2,combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
f2c<-cowplot::plot_grid(plotlist = p,nrow = 1)



features <- c("IL1B", "IL1A","TNF","CCL4","CCL2")

p<-FeaturePlot(imm_all, features = features,reduction = "tsne", 
               min.cutoff = "q05", max.cutoff = "q95", ncol = 5,col=plt2,pt.size = 0.2,combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
f2d<-cowplot::plot_grid(plotlist = p,nrow = 1)
#####
p<-FeaturePlot(velmeshev, features = features,reduction = "tsne", 
               min.cutoff = "q05", max.cutoff = "q95", ncol = 5,col=plt2,pt.size = 0.2,combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
f2e<-cowplot::plot_grid(plotlist = p,nrow = 1)

p<-FeaturePlot(masuda, features = features,reduction = "tsne", 
               min.cutoff = "q05", max.cutoff = "q95", ncol = 5,col=plt2,pt.size = 0.2,combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
f2f<-cowplot::plot_grid(plotlist = p,nrow = 1)

lay <- rbind(c(1,2),
             c(NA,3),
             c(4,4),
             c(5,5),
             c(6,6))
figure2<-grid.arrange(grobs=list(f2a,f2b,f2c,f2d,f2e,f2f), ncol=2, widths = 1:2, 
             layout_matrix=lay)

figure2


################################ Figure 2 Code ends ###########################


############################# Figure 3 ########################################
                      #### Figure 3a #####
imm<-readRDS("~/NatNeu2022/data/immune_cells_16_cluster_figure3.rds")

## Figure 3a
cols<-c(brewer.pal(12, "Set3"),brewer.pal(n = 8, name = "Dark2"),"#FDB462",
        "#0000CD","#7B68EE","#4682B4","#808000","#7CFC00") 
f3a<-DimPlot(imm,label=TRUE,reduction = "tsne",cols= cols)+NoLegend()

                      #### Figure 3b #####
plt2<-c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080")
features<- c("CD3-Ab","CD4-Ab","CD11b-Ab","CD14-Ab","CD8a-Ab","CD19-Ab",
             "HLA-DR-Ab","CD20-Ab","CD16-Ab","CD56-Ab")
p<-FeaturePlot(imm, features = features,reduction = "tsne", 
               min.cutoff = "q05", max.cutoff = "q95", ncol = 2,cols = plt2,combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

f3b<-cowplot::plot_grid(plotlist = p,ncol = 2)


#### Figure 3c #####
pairFreqDF_to_network<-function(pairFreqDF,collist = c("gray","gold","tomato"),cutpoints = c(0,1,5)){
  
  c_to_c<-as.character(pairFreqDF$pair_list)
  enL<-rep("F",length(c_to_c))
  enR<-rep("T",length(c_to_c))
  
  for(i in 1:length(enL)){
    enL[i]<-strsplit(c_to_c[i],split = "_")[[1]][1]
    enR[i]<-strsplit(c_to_c[i],split = "_")[[1]][2]
  }
  edgeDF<-data.frame(from=enL,to=enR)
  edgeDF$arrows<-rep("to",nrow(edgeDF))
  
  pair_freq<-pairFreqDF$Freq
  b<-c(cutpoints,max(pair_freq))
  names<-collist
  weightCol<-cut(pair_freq, breaks = b, labels = names)
  
  edgeDF$color<-weightCol
  edgeDF$value<-sqrt(pair_freq)
  id<-unique(c(enL,enR))
  group<-rep("g",length(id))
  group[id %in% unique(enL)]<-"Lig"
  group[id %in% unique(enR)]<-"Rec"
  nodesDF<-data.frame(id=id,group=group,label=id,title=id)
  return(list(node=nodesDF,edge=edgeDF))
}
######

nDF<-readRDS("./NVU_to_imm_network.rds")
visNetwork(nodes = nDF$nodeList,edges = nDF$edgeList) %>%
  visIgraphLayout(randomSeed = 1234)%>%
  visNodes(font=list(size=26,color="teal"),size=30)%>%
  visLegend(main = "network")

#### Figure 3d #####
pDF<-pairFreqDF_to_network(pairFreqDF = nDF$pairFreqDF,
                           collist = c("gray","gold","tomato","teal"),
                           cutpoints = c(0,1,5,10))


visNetwork(nodes = pDF$node,edges = pDF$edge) %>%
  visIgraphLayout(randomSeed = 1234)%>%
  visNodes(font=list(size=30,color="teal"),size=30)%>%
  visLegend(main = "network")


################################ Figure 3 Code ends ###########################

###############################################################################
################################ Figure 4 Code  ###############################

library(pheatmap)
library(RColorBrewer)

sampleInfo<-readRDS("./sampleInfo.rds")
heatmap_data<-readRDS("./figure4_heatmap_data.rds")
#mmTohg<-readRDS("~/neuro10x/PublicDataSets/mouseData/MouseEpilepsy/mmToHg.rds")
logcpm_comm_l<-heatmap_data$logcpm_ligands
rmdupl_index<-which(duplicated(logcpm_comm_l$gene_name))
logcpm_comm_l$gene_name[which(duplicated(logcpm_comm_l$gene_name))]
logcpm_comm_l<-logcpm_comm_l[- rmdupl_index,]
rownames(logcpm_comm_l)<-logcpm_comm_l$gene_name # for Receptor it did not had any duplicates

all(sampleInfo$accNum==colnames(logcpm_comm_l)[4:ncol(logcpm_comm_l)])

logcpm_comm_l<-logcpm_comm_l[order(logcpm_comm_l$HGNC),]
HGNC<-paste0(as.character(logcpm_comm_l$gene_name),"::",as.character(logcpm_comm_l$HGNC))
logcpm_mat<-logcpm_comm_l[,4:ncol(logcpm_comm_l)]


ligand_heatmap<-pheatmap(as.matrix(logcpm_mat),colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(7),
                         scale = "row",cluster_rows = FALSE,cluster_cols = FALSE,labels_row = HGNC,
                         labels_col =sampleInfo$group )


logcpm_comm_r<-heatmap_data$logcpm_r
rownames(logcpm_comm_r)<-logcpm_comm_r$gene_name 

all(sampleInfo$accNum==colnames(logcpm_comm_r)[4:ncol(logcpm_comm_r)])

logcpm_comm_r<-logcpm_comm_r[order(logcpm_comm_r$HGNC),]
HGNC<-paste0(as.character(logcpm_comm_r$gene_name),"::",as.character(logcpm_comm_r$HGNC))
logcpm_mat<-logcpm_comm_r[,4:ncol(logcpm_comm_r)]

receptor_heatmap<-pheatmap(as.matrix(logcpm_mat),colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(7),
                           scale = "row",cluster_rows = FALSE,cluster_cols = FALSE,labels_row = HGNC,
                           labels_col =sampleInfo$group )


################################ Figure 4 Code ends  ###############################

################################ Figure 5 Code  ###############################

### Fig 5a
c3_clusters <-readRDS("./C3_reclusters.rds")

DimPlot(c3_clusters,reduction = "tsne_adt",label = T)

### Fig 5b
plt2<-c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080")
features<-c("CD4-Ab","CD8a-Ab","CD3-Ab","CD19-Ab","CD11b-Ab","CD16-Ab")
p<-FeaturePlot(c3_clusters,features = features,reduction = "tsne_adt",
               cols = plt2,pt.size = 0.5,ncol = 3,combine = F)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

f5b<-cowplot::plot_grid(plotlist = p,ncol = 3)

f5b
# Fig f5c
features<-c("NKG7","GZMB","GNLY")
p<-FeaturePlot(c3_clusters,features = features,reduction = "tsne_adt",
               cols = plt2,pt.size = 0.5,ncol = 2,combine = F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

f5c<-cowplot::plot_grid(plotlist = p,ncol = 2)

f5c
##### Fig 5d
##### plot gene contribution figures ######
mle_res<-readRDS("./mle_res_cd4.rds")
sig_gene_id <-readRDS("./cd4_sig_gene_id_dataframe.rds")
pbmc.PIC<-readRDS("./cells_for_pic_seq_analysis.rds")

##############
# Compute the expected transcription of PICs
good_pics <-rownames(mle_res)
alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
split_count = cbind(alpha[good_pics], 1 - alpha[good_pics])# * umicount[good_pics]
split_count = na.omit(split_count)
#split_count = split_count / rowSums(split_count)
t_col<-'gold'; dc_col<-'darkblue'

split_count_ord <- split_count[order(split_count[,1]),] # order by T cell mixing
#barplot(t(split_count_ord), col = c(t_col, dc_col), border = NA, xaxs = "i", space = 0, names.arg = rep("", nrow(split_count)), las = 2)

set.seed(5566187)

cells_CD4 <- sample(names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="CD4T"]),250)
#cells_CD8 <- sample(names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="CD8T"]),250)
cells_Micr <- sample(names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="Micr"]),250)
sample_cells <- c(cells_CD4,cells_Micr)
#sample_cells <- c(cells_CD8,cells_Micr)

d_CD4<- names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="d_CD4"])
#d_CD4<- names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="d_CD8"])
d_CD4 <- intersect(rownames(split_count),d_CD4)

all(rownames(split_count)==d_CD4)

features <- c(head(sig_gene_id,10)$genes, tail(sig_gene_id,10)$genes)

pbmc.PIC<-ScaleData(pbmc.PIC,features = features)

## fig 5d
DoHeatmap(pbmc.PIC,features = features,cells = rownames(split_count_ord),draw.lines = FALSE)+
  scale_fill_gradientn(colors = c("gray", "white", "red"))

barplot(t(split_count_ord), col = c(t_col, dc_col), border = NA, xaxs = "i", space = 0, names.arg = rep("", nrow(split_count)), las = 2)

#### Fig 5e
DoHeatmap(pbmc.PIC,features = features,cells = sample_cells,draw.lines = FALSE)+
  scale_fill_gradientn(colors = c("gray", "white", "red"))


########################## Fig5f, fig5g

mle_res<-readRDS("./mle_res_cd8.rds")
sig_gene_id <-readRDS("./cd8_sig_gene_id_dataframe.rds")
pbmc.PIC<-readRDS("./cells_for_pic_seq_analysis.rds")

##############
# Compute the expected transcription of PICs
good_pics <-rownames(mle_res)
alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
split_count = cbind(alpha[good_pics], 1 - alpha[good_pics])# * umicount[good_pics]
split_count = na.omit(split_count)
#split_count = split_count / rowSums(split_count)
t_col<-'gold'; dc_col<-'darkblue'

split_count_ord <- split_count[order(split_count[,1]),] # order by T cell mixing

set.seed(5566187)

cells_CD8 <- sample(names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="CD8T"]),250)
cells_Micr <- sample(names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="Micr"]),250)
sample_cells <- c(cells_CD8,cells_Micr)

d_CD4<- names(Idents(pbmc.PIC)[Idents(pbmc.PIC)=="d_CD8"])
d_CD4 <- intersect(rownames(split_count),d_CD4)

all(rownames(split_count)==d_CD4)

features <- c(head(sig_gene_id,10)$genes, tail(sig_gene_id,10)$genes)

pbmc.PIC<-ScaleData(pbmc.PIC,features = features)

## fig 5f
DoHeatmap(pbmc.PIC,features = features,cells = rownames(split_count_ord),draw.lines = FALSE)+
  scale_fill_gradientn(colors = c("gray", "white", "red"))

barplot(t(split_count_ord), col = c(t_col, dc_col), border = NA, xaxs = "i", space = 0, names.arg = rep("", nrow(split_count)), las = 2)

#### Fig 5g
DoHeatmap(pbmc.PIC,features = features,cells = sample_cells,draw.lines = FALSE)+
  scale_fill_gradientn(colors = c("gray", "white", "red"))

################################ Figure 5 Code ends  ###############################

