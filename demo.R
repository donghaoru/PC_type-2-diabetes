#
# GEO database
#
library(GEOquery)
# GSE76894   GPL570   hgu133plus2
# GSE76895   GPL570
eSet <- getGEO('GSE118139', destdir=".",
               AnnotGPL = F,
               getGPL = F)
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
raw_exprSet <- exprs(eSet[[1]])
phe <- pData(eSet[[1]])

library(stringr) 
Group <- as.data.frame(str_split(phe$source_name_ch1,',',simplify = T)[,1]) 
Group <- as.data.frame(c(rep('ND',50),rep('T2D',50)))
colnames(Group) <- 'Group' 
rownames(Group) <- rownames(phe)

BiocManager::install('hgu133plus2.db') 
library(hgu133plus2.db)
ids <- toTable(hgu133plus2SYMBOL)

dat <- raw_exprSet
dat_id <- dat[rownames(dat) %in% ids$probe_id,] 
ids2 <- ids[match(rownames(dat_id),ids$probe_id),] 
dat_id <- as.matrix(dat_id)
rownames(dat_id) <- ids2$symbol

library(limma) 
dimnames <- list(rownames(dat_id),colnames(dat_id))
d <- matrix(as.numeric(as.matrix(dat_id)),nrow=nrow(dat_id),dimnames=dimnames) 
dat_final <- as.data.frame(avereps(d)) 

library(reshape2)
melt_dat<-melt(dat_final)
colnames(melt_dat)=c('Sample','Value')
group_list <- Group$Group
melt_dat$Group=rep(group_list,each=nrow(dat_final))
ggplot(melt_dat,aes(x=Sample,y=Value,fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(3,14))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

library(limma)
dat <- dat_final
design <- model.matrix(~0+factor(Group$Group)) 
colnames(design)=levels(factor(Group$Group)) 
rownames(design)=colnames(dat) 
contrast.matrix <- makeContrasts(T2D - ND,
                                 levels = design)
fit1 <- lmFit(dat,design)
fit2 <- contrasts.fit(fit1, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf) 
allDEG = na.omit(tempOutput)

gene <- c('PC','PDHA1','PDHB','PDHX','PDK1','CS','ACO2','IDH1','IDH2','OGDH','SUCLG2','SUCLG1','SUCLA2','SDHB','FH','MDH2','SLC2A1','HK1','HK2','GPI','ALDOA','TPI1','GAPDH','PGK1','PGM1','ENO1','ENO2','PKM2','PFKP','PFKL','PFKM','MDH1','ME1','ME2','ME3','SLC25A11','PCCA','PCCB','ECH1','GLS','LDHA','LDHB','ACLY')
a <- intersect(gene,rownames(allDEG))
dat1 <- allDEG[a,]
library(ggpubr)
library(ggthemes)
dat1$Group = as.factor(ifelse(dat1$adj.P.Val < 0.05 & dat1$logFC <= -log2(2), 'Down',
                              ifelse(dat1$adj.P.Val < 0.05 & dat1$logFC >= log2(2), 'Up', 'Not')))
dat1$v=-log2(dat1$adj.P.Val)
dat1$X <- rownames(dat1)
ggscatter(dat1, legend = "right",ylim = c(0, 10), xlim=c(-3,6),
          x = "logFC", y = "v",
          xlab = 'LogFC', ylab="-Log2(FDR)",
          color = "Group", size = 1, 
          palette = c("#00AFBB", "#999999", "#FC4E07"),
          label = 'X',
          repel = T,
          label.select =rownames(dat1[abs(dat1$logFC) > log2(1.2) & dat1$adj.P.Val < 0.05,])
)+
  theme_base()+
  geom_hline(yintercept = -log2(0.05),linetype='dashed')+
  geom_vline(xintercept = c(-log2(2),log2(2)),linetype='dashed')

save(allDEG1,allDEG2,file = 'DEG.Rdata')


#
Group$y <- ifelse(Group$Group=='T2D',NA,'a')
Group <- na.omit(Group)
data <- dat_final[,rownames(Group)]
Group$pc <- as.data.frame(t(data['PC',]))$PC
Group$Group <- ifelse(Group$pc < median(Group$pc), 'Low','High')
Group <- Group[order(Group$Group),]
data <- data[,rownames(Group)]

g <- as.data.frame(Group[,1])
colnames(g) <- 'Group'
rownames(g) <- rownames(Group)

library(limma)
dat <- data
Group <- g
design <- model.matrix(~0+factor(Group$Group)) 
colnames(design)=levels(factor(Group$Group)) 
rownames(design)=colnames(dat) 
contrast.matrix <- makeContrasts(High - Low,
                                 levels = design)
fit1 <- lmFit(dat,design)
fit2 <- contrasts.fit(fit1, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf) 
allDEG = na.omit(tempOutput)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GSEABase)
d='D:/app/geneset'
geneset=read.gmt(file.path(d,'h.all.v7.5.1.symbols.gmt')) #c2.cp.kegg.v7.5.1.symbols.gmt c5.go.bp.v7.5.1.symbols.gmt
genelist=allDEG$logFC
names(allDEG)=rownames(allDEG)
genelist=sort(genelist,decreasing = T)

gsea.hall1=GSEA(genelist,TERM2GENE = geneset,verbose = F)
gseaplot2(gsea.hall, geneSetID = 23,pvalue_table = T,color = '#00AFBB')

save(gsea.hall1,gsea.hall2,file = 'pc_analysis.Rdata')

##  single cell
library(Seurat)
library(tidyverse)
srt <- CreateSeuratObject(counts = data,
                          meta.data = group,
                          min.cells = 3, 
                          min.features = 200)
srt[["percent.mt"]] <- PercentageFeatureSet(srt, 
                                            pattern = "^MT-")

VlnPlot(object = srt, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "Group",
        log = T,
        pt.size = 0.1)

FeatureScatter(object = srt, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = "Group")

srt <- subset(srt,
              subset = nFeature_RNA > 1000 &
                percent.mt > -Inf & # 极端值
                percent.mt < 10)
srt <- CellCycleScoring(
  srt, 
  s.features = cc.genes$s.genes, 
  g2m.features = cc.genes$g2m.genes)
head(srt2@meta.data[ , 5:10])

save(dat_final,Group,srt2,file = 'smart2.Rdata')

srt2 <- NormalizeData(object = srt2,
                      normalization.method = "LogNormalize")

srt2 <- FindVariableFeatures(object = srt2, 
                             selection.method = "vst", 
                             mean.function = ExpMean,
                             dispersion.function = LogVMR,
                             nfeatures = 2000 # 选择???变基因的数???，1000-5000
)
srt2<- ScaleData(object = srt2,
                 features = rownames(srt2))
HVFInfo(object = srt2)[1:15,1:3]
top10 <- head(VariableFeatures(srt2),
              10)
plot1 <- VariableFeaturePlot(srt2)
LabelPoints(plot = plot1,
            points = top10,
            repel = TRUE)
srt2 <- RunPCA(srt2)
ElbowPlot(srt2,
          ndims = 50) 
srt2 <- RunTSNE(srt2, 
                dims = 1:8) 
DimPlot(object = srt2, 
        reduction = "tsne", 
        group.by = "Group",
        dims = c(1,2),
        shuffle = TRUE,
        label = TRUE,
        label.size = 4,
        label.color = "black",
        label.box = F,
        sizes.highlight = 1)
# UMAP分析，t-SNE的完美替代品，速度快，结果解释多
srt2 <- RunUMAP(srt2, 
                dims = 1:8)
DimPlot(object = srt2, 
        reduction = "umap", 
        group.by = "Group",
        dims = c(1,2),
        shuffle = TRUE,
        label = TRUE,
        label.size = 4,
        label.color = "black",
        label.box = TRUE,
        sizes.highlight = 1)

FeaturePlot(srt2, 
            reduction = "umap", 
            features = "PC")

VlnPlot(srt2,features = 'PC',group.by = "Group")

srt2 <- FindNeighbors(srt2,
                      k.param = 20,
                      dims = 1:8)
srt2 <- FindClusters(srt2,
                     resolution = 0.2, # 最重要参数，该值越???，cluster 越多
                     method = "igraph", # 根据结果调整
                     algorithm = 1, # 根据结果调整'
                     random.seed = 2022)
UMAPPlot(object = srt2,
         group.by = "seurat_clusters",
         pt.size = 2, 
         label = TRUE)
TSNEPlot(object = srt2,
         group.by = "seurat_clusters",
         pt.size = 2,
         label=T)

library(SingleR)
library(celldex)
library(BiocParallel)
ref <- celldex::HumanPrimaryCellAtlasData()
load('HumanPrimaryCellAtlasData.Rda')
singler <- SingleR(test = srt2@assays$RNA@data,
                   ref = ref,
                   labels = ref$label.main,
                   clusters = srt2@meta.data$seurat_clusters,
                   fine.tune = TRUE, 
                   BPPARAM = SnowParam(8))
singler$pruned.labels
cell_ids <- singler$pruned.labels
names(cell_ids) <- levels(srt2)

srt2 <- RenameIdents(srt2, cell_ids)
UMAPPlot(object = srt2, pt.size = 0.5, label = TRUE)

cell_marker <- read.table('cellmarker.txt',header = T,sep = '\t')
DotPlot(srt2, 
        features = unique(cell_marker$Cell.Marker),group.by = "seurat_clusters") +
  theme(axis.text = element_text(size = 8,
                                 angle = 90,
                                 hjust = 1))
DimPlot(srt2, 
        reduction = "umap", # pca, umap, tsne
        group.by = "Cell.tytle",
        label = T)

Cell_type <- c("0" = "Alpha cell",
               "1" = "Beta cell",
               "2" = "Alpha cell1",
               "3" = "Ductal cell",
               "4" = "Acinar cell",
               "5" = "Alpha cell2",
               '6'='Stellate cells',
               '7'='Alpha cell3')
srt2[['cell_type']] <- unname(Cell_type[srt2@meta.data$seurat_clusters])
DimPlot(srt2, 
        reduction = "umap", 
        group.by = "cell_type",
        label = TRUE, 
        pt.size = 2) +
  NoLegend()
FeaturePlot(srt2, 
            reduction = "umap", 
            features = c("PC"),
            label = TRUE)
VlnPlot(srt2,
        features = c("PC"),
        group.by = 'cell_type')
target.cell <- row.names(srt@meta.data)[which(srt@meta.data$cell_type == "Mature beta cell")]
B.data <- subset(srt, cells = target.cell)
Idents(B.data) <- B.data@meta.data$PC_state
levels(B.data)
Idents(mbc) <- mbc@meta.data$PC_state
DEG <- FindMarkers(mbc,
                   ident.1 = 'Positive',
                   ident.2 = "Negative",
                   logfc.threshold =0 )

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GSEABase)
d='D:/app/geneset'
geneset1=read.gmt(file.path(d,'c2.cp.kegg.v7.5.1.symbols.gmt')) #c2.cp.kegg.v7.5.1.symbols.gmt c5.go.bp.v7.5.1.symbols.gmt  h.all.v7.5.1.symbols.gmt
genelist=DEG$avg_log2FC
names(genelist)=rownames(DEG)
genelist=sort(genelist,decreasing = T)

gsea.KEGG=GSEA(genelist,TERM2GENE = geneset1,verbose = F,pvalueCutoff = 1)

a <- as.data.frame(gsea.KEGG)
b <- a[1:10,]
b$type <- ifelse(b$NES < 0, 'Up','Down')
ggplot(a,aes(x = NES,y = ID,fill=type))+
  geom_bar(stat="identity", width=0.5)+
  NoLegend()+
  ylab('')+
  labs(title = "The Most Enriched KEGG Terms")+
  scale_y_discrete(limits=c("KEGG_ERBB_SIGNALING_PATHWAY" ,"KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION","KEGG_CARDIAC_MUSCLE_CONTRACTION","KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_CHONDROITIN_SULFATE","KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY" ,"KEGG_RNA_DEGRADATION","KEGG_CELL_CYCLE","KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450","KEGG_P53_SIGNALING_PATHWAY"))
ggplot(b,aes(x = NES,y = ID,fill=type))+
  geom_bar(stat="identity", width=0.5)+
  NoLegend()+
  ylab('')+
  labs(title = "Top 10 Enriched KEGG Terms")+
  scale_y_discrete(limits=c("KEGG_CARDIAC_MUSCLE_CONTRACTION","KEGG_AMINOACYL_TRNA_BIOSYNTHESIS" ,"KEGG_NICOTINATE_AND_NICOTINAMIDE_METABOLISM" ,"KEGG_SNARE_INTERACTIONS_IN_VESICULAR_TRANSPORT","KEGG_CELL_ADHESION_MOLECULES_CAMS","KEGG_JAK_STAT_SIGNALING_PATHWAY","KEGG_INSULIN_SIGNALING_PATHWAY","KEGG_RETINOL_METABOLISM","KEGG_TYPE_II_DIABETES_MELLITUS","KEGG_P53_SIGNALING_PATHWAY"))

save(srt,DEG,gsea.KEGG,file = 'YM_sc.Rdata')


# PC-KO

data <- read.table('DATA.txt',sep = '\t',header = T)
library(limma)
rt <- as.matrix(data)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
dat <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dat <- avereps(dat)
dat <- as.data.frame(dat)

library(limma)
library(ggplot2)
library(reshape2)
melt_dat<-melt(data)
colnames(melt_dat)=c('Gene','Sample','Value')
group_list <- Group$Group
melt_dat$Group=rep(group_list,each=nrow(dat))
ggplot(melt_dat,aes(x=Sample,y=Value,fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(0,14))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

Group=as.data.frame(c(rep("WT",3),rep("KO",3)))
colnames(Group) <- 'Group'
rownames(Group) <- colnames(dat)

design <- model.matrix(~0+factor(Group$Group))
colnames(design) <- levels(factor(Group$Group))
fit <- lmFit(dat,design)
cont.matrix<-makeContrasts(KO-WT,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
diff <- topTable(fit2,p.value = 0.05,lfc = 1,number = 10000000000 )
DEG.id <- mapIds(x = org.Mm.eg.db,
                 keys = rownames(diff),
                 keytype = 'SYMBOL',
                 column = 'ENTREZID')
DEG.id <- na.omit(DEG.id)
KEGG <- enrichKEGG(gene = DEG.id,
                   organism = "mmu",
                   pvalueCutoff =0.05,
                   qvalueCutoff = 1)
aa <- as.data.frame(KEGG)
diff_kegg <- aa[c('mmu04668', 'mmu04210', 'mmu04350', 'mmu04115', 'mmu04915', 'mmu05330', 'mmu04610', 'mmu00980'),]
p = ggplot(diff_kegg,aes(x = GeneRatio,y = Description))
p=p + geom_point()  
p=p + geom_point(aes(size=Count))
pbubble = p+ geom_point(aes(size=Count,color=-1*log10(pvalue)))
pr = pbubble+scale_color_gradient(low="#00AFBB",high = "#FC4E07")
pr = pr+labs(color=expression(-log[10](Pvalue)),size="Count",  
             x="GeneRatio",y="Pathway name")
pr + theme_bw()+scale_y_discrete(limits=c("Hippo signaling pathway - multiple species","Protein export" ,"Citrate cycle (TCA cycle)","DNA replication","Ferroptosis","Fatty acid metabolism","Fc epsilon RI signaling pathway","GnRH secretion","Adipocytokine signaling pathway","p53 signaling pathway"))



