options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages('installr')
library(installr)
install.Rtools()

options(download.file.method ='libcurl')#BiocManager下载出错的时候用这个，每次下载都要重新运行这个代码
options(url.method='libcurl')
BiocManager::install("annotate")
BiocManager::install("IRanges")
BiocManager::install("Biobase")
BiocManager::install("BiocGenerics")
BiocManager::install("annotate")
BiocManager::install("IRanges")
BiocManager::install("annotate")
BiocManager::install("IRanges")
BiocManager::install("limma")
BiocManager::install("Seurat")

SummarizedExperiment
library(scater)
library(GSEABase)

#首先，从GSE131685中下载数据：里面的文件名要分别改为“barcodes.tsv”、“genes.tsv”和“matrix.mtx”，在Read10X
#导入10X和SmartSeq2数据Tabula Muris）时才不会报错。。。

rm(list = ls())  ## 魔幻操作，一键清空~
Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor
memory.limit(10000000) #改变R内存配置上限

#################################0.打开R包并载入数据###################################################
library(devtools)
install_github("immunogenomics/harmony")
library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)

setwd("C:\\Users\\muyumuge\\Desktop\\normal lung\\all normal")   
Lun1 <- readRDS("C:/Users/muyumuge/Desktop/normal lung/all normal/Lun1.rds")

#Lung data loading 并构建seurat object

L1.data <- Read10X(data.dir = "C:\\Users\\muyumuge\\Desktop\\normal lung\\all normal\\Lung1\\")
L1 <- CreateSeuratObject(counts = L1.data, project = "Lung1", min.cells = 8, min.features = 200)
L2.data <- Read10X(data.dir = "C:\\Users\\muyumuge\\Desktop\\normal lung\\all normal\\Lung2\\")
L2 <- CreateSeuratObject(counts = L2.data, project = "Lung2", min.cells = 6, min.features = 200)
L3.data <- Read10X(data.dir = "C:\\Users\\muyumuge\\Desktop\\normal lung\\all normal\\Lung3\\")
L3 <- CreateSeuratObject(counts = L3.data, project = "Lung3", min.cells = 8, min.features = 200)
Luntest <- merge(x = L1, y = list(L2,L3,L4,L5,L6,L7)) #读取文件并用merge函数进行合并

L7 <- CreateSeuratObject(counts = lung7, project = "Lung7", min.cells = 6, min.features = 200)
Lun=Lun1
saveRDS(Luntest,file="C:\\Users\\muyumuge\\Desktop\\normal lung\\all normal\\Luntest.rds")

############################################1.数据前期处理和矫正##############################################
Lun=Lun1#重新命名一个名称避免被修改后无法返回重做

# quality control,使用PercentageFeatureSet函数计算线粒体基因的百分比
Lun[["percent.mt"]] <- PercentageFeatureSet(Lun, pattern = "^MT-") #提取有关线粒体的基因
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 人类血液常见红细胞基因
HB_m <- match(HB.genes_total,rownames(Lun@assays$RNA))
HB.genes <- rownames(Lun@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
Lun[["percent.HB"]]<-PercentageFeatureSet(Lun,features=HB.genes)

pdf(file="1.beforeQC.featureViolin.pdf",width=10,height=6)
VlnPlot(Lun, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), ncol = 4) #由图可以看出分布还可以,基因特征的小提琴图
dev.off()

#测序深度的相关性绘图
pdf(file="1.beforeQC.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "percent.HB")
CombinePlots(plots = list(plot1, plot2,plot3))
dev.off()

#筛选条件,对数据进行过滤
Lun <- subset(Lun, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20& percent.HB <10) 

#过滤后基因特征的小提琴图
pdf(file="1.afterQC.featureViolin.pdf",width=10,height=6)
VlnPlot(Lun, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), ncol = 4) #由图可以看出分布还可以,基因特征的小提琴图
dev.off()

#过滤后测序深度的相关性绘图
pdf(file="1.afterQC.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "percent.HB")
CombinePlots(plots = list(plot1, plot2,plot3))
dev.off()

#对数据进行标准化
Lun <- NormalizeData(Lun, normalization.method = "LogNormalize", scale.factor = 10000)
Lun <- NormalizeData(Lun) #标准化

#提取那些在细胞间变异系数较大的基因
Lun <- FindVariableFeatures(Lun, selection.method = "vst", nfeatures = 2000) #查找高变基因
top10 <- head(VariableFeatures(Lun), 10)

#标准化后测序深度的相关性绘图
pdf(file="4.NormalizeQC.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Lun, feature1 = "nCount_RNA", feature2 = "percent.HB")
CombinePlots(plots = list(plot1, plot2,plot3))
dev.off()

pdf(file="1.after.normalizeQC.featureViolin.pdf",width=10,height=6)
VlnPlot(Lun, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), ncol = 4) #由图可以看出分布还可以,基因特征的小提琴图
dev.off()

#输出特征方差图
pdf(file="5.NormalizeQC.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(Lun)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

###################################2.分析细胞周期对分群的影响并去掉细胞周期相关基因###################################
#这里我想叨叨几句，据我看到的文献，多数是在进行降维后将细胞周期方面对分群的影响作为一个单独模块去叙述，作者在先期不管
#细胞周期对聚类是否有影响的情况下就对细胞周期相关基因进行去除也是比较明智的，因为作者并不想让该因素混杂其中影响分群
# 计算细胞周期
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
Lun <- CellCycleScoring(Lun, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all.genes <- rownames(Lun)
Lun <- ScaleData(Lun, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes) #等候大约1小时,这一步需要很大内存，memory limit(100000)

#当然我们还是要看是否细胞周期真的有影响，感兴趣的小伙伴可以看一下，确实是有一定影响的！
#Lun <- ScaleData(Lun, features = rownames(Lun))
#Lun  <- RunPCA(Lun , features = c(s.genes, g2m.genes))
#DimPlot(Lun)

##############################################3.PCA主成分分析################################################
##PCA分析，这一步是为了提取所有细胞的特征，线性降维，和后面的非线性降维不一样。非线性降维是为了可视化，更好看。
Lun=ScaleData(Lun)                     #PCA降维之前的标准预处理步骤
Lun=RunPCA(object= Lun,npcs = 20,pc.genes=VariableFeatures(object = Lun), verbose = F)     #PCA分析,经过上面一步的
#分析之后，Lun数据的active.ident已经变为了G1,G2M,S，renameident函数可以更改active.ident。

#绘制每个PCA成分的相关基因
pdf(file="03.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = Lun, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#主成分分析图形
pdf(file="03.PCATSNE1.pdf",width=6.5,height=6)
DimPlot(object = Lun08, reduction = "pca")
dev.off()

#主成分分析热图
pdf(file="03.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = Lun, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#每个PC的p值分布和均匀分布
Lun <- JackStraw(object = Lun, num.replicate = 100)
Lun <- ScoreJackStraw(object = Lun, dims = 1:20)
pdf(file="03.pcaJackStrawTSNE.pdf",width=10,height=8)
plot1 <- JackStrawPlot(object = Lun, reduction="pca", dims = 1:20,ymax = 0.8)
plot2 <- ElbowPlot(Lun)
CombinePlots(plots = list(plot1,plot2),legend = "bottom")
dev.off()
#JackStrawPlot虚线以上的为可用维度，你也可以调整 dims 参数，画出所有 pca 查看
#ElbowPlot(Lun),这是另外一种鉴定手段是绘制所有PC的分布点图,大多数软件都是通过拾取拐点处的pc作为选定数目
#这一步如果所有维度都有意义就选择全部，如果，   



#####################################4.利用harmony算法去除批次效应并细胞分簇#############################
#Eliminate batch effects with harmony and cell classification,在harmony之前需要对seruat对象进行常规的PCA分析
options(repr.plot.height = 2.5, repr.plot.width = 6)

Lun <- Lun %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)　#等候时间大约7分钟
harmony_embeddings <- Embeddings(Lun, 'harmony')
harmony_embeddings[1:5, 1:5]
#第一种UMAP降维聚类：
Lun <- Lun %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.25) %>%
  identity()        #Number of communities: 25，这一步会计算出分为了多少个cluster，下面的数字相应的改正
new.cluster.ids <- c(0,1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16) #这里根据上面的多少个cluster做相应的改正
names(new.cluster.ids) <- levels(Lun)
Lun <- RenameIdents(Lun, new.cluster.ids)
#第二种TSNE降维聚类：
Luntsne<- Lun %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.25) %>%
  identity()    #两种降维方法得出的分簇是不一样的，Number of communities: 17，这一步会计算出分为了多少个cluster，下面的数字相应的改正
new.cluster.ids0 <- c(0,1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16) #这里根据上面的多少个cluster做相应的改正
names(new.cluster.ids0) <- levels(Luntsne)
Luntsne <- RenameIdents(Luntsne, new.cluster.ids0) #从这里开始形成两套seurat对象。

Lun=Lun08 #重命名为了方便下面的运算


##########################################5.鉴定Marker基因############################################
#Calculating differentially expressed genes (DEGs) and Save rds file
logFCfilter=0.5
adjPvalFilter=0.05
Lun.markers <- FindAllMarkers(object = Lun,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter) #寻找高变基因，等候时间大约20分钟
sig.markers=Lun.markers[(abs(as.numeric(as.vector(Lun.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(Lun.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="05.markersTSNE.xls",sep="\t",row.names=F,quote=F)
saveRDS(Lun,file="05LunTSNE.rds")
#找出每个簇的前十个基因
top10 <- Lun.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#绘制marker在各个cluster的热图
pdf(file="06.tsneHeatmapTSNE.pdf",width=12,height=9)
DoHeatmap(object = Lun, features = top10$gene) + NoLegend() #这个可以用其他热图函数重新画，添加分组信息https://www.jianshu.com/p/b2816a72e79e
dev.off()

#寻找每个簇中显著表达的基因,3种函数画图表示FeaturePlot，VlnPlot，DotPlot
cluster0.markers <- FindMarkers(object = Lun, ident.1 = 0, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster0.markers, n = 9)

      #1绘制marker在各个cluster的散点图
pdf(file="05.markerScatterEpiU.pdf",width=10,height=6)
FeaturePlot(object = Lun, reduction = "umap", features = c("SLPI", "SFTPC", "SFTPB", "PGC", "AGER", "NAPSA","SCGB3A2", "SFTA2", "SFTPD")) 
dev.off()
pdf(file="05.markerScatterEpiT.pdf",width=10,height=6)
FeaturePlot(object = Lun, reduction = "tsne", features = c("SLPI", "SFTPC", "SFTPB", "PGC", "AGER", "NAPSA","SCGB3A2", "SFTA2", "SFTPD")) 
dev.off()

      #2绘制marker的小提琴图
pdf(file="05.markerViolinEpi.pdf",width=10,height=6)
VlnPlot(object = Lun, pt.size=0,combine=T,features = c("SLPI", "SFTPC", "SFTPB", "PGC", "AGER", "NAPSA","SCGB3A2", "SFTA2", "SFTPD")) 
dev.off()

      #3绘制marker在各个cluster的气泡图
pdf(file="05.markerBubbleEpi.pdf",width=12,height=6)
DotPlot(object =Lun, features = c("SLPI", "SFTPC", "SFTPB", "PGC", "AGER", "NAPSA","SCGB3A2", "SFTA2", "SFTPD")) #考虑如何将多表达和少表达的放一起??
dev.off()

cluster1.markers <- FindMarkers(object = Lun, ident.1 = 1, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster22.markers, n = 9)

cluster2.markers <- FindMarkers(object = Lun, ident.1 = 2, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster2.markers, n = 9)

cluster3.markers <- FindMarkers(object = Lun, ident.1 = 3, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster3.markers, n = 9)

cluster4.markers <- FindMarkers(object = Lun, ident.1 = 4, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster4.markers, n = 9)

cluster5.markers <- FindMarkers(object = Lun, ident.1 = 5, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster5.markers, n = 9)

cluster6.markers <- FindMarkers(object = Lun, ident.1 = 6, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster6.markers, n = 9)

cluster7.markers <- FindMarkers(object = Lun, ident.1 = 7, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster7.markers, n = 9)

cluster7.1.markers <- FindMarkers(object = Lun, ident.1 = 7,ident.2 =6, min.pct = 0.25) #这里需要根据每个Custer进行相应的参数设置
head(x = cluster7.1.markers, n = 9)


#这样是寻找单个聚类中的显著基因，ident.2 = c(0, 3)代表以0和3两个簇作为参照选择maker，只有ident.1则代表以全部簇作为参照
cluster5.markers <- FindMarkers(object = Lun, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#ident.2:A second identity class for comparison; if NULL, use all other cells for comparison; if an object of class phylo or 'clustertree' is passed to ident.1, must pass a node to find markers for
head(x = cluster5.markers, n = 5)
#这样寻找所有聚类中显著基因，计算速度很慢，需要等待，另外，我们有多种方法统计基因的显著性
FeaturePlot(object = Lun, features = c("NKG7", "GZMB", "GNLY", "TYROBP", "CCL3", "KLRD1")) #这里就是根据上面得到的top5基因改名字就行，n可以随意改变


#############################################06.umap和tsne可视化###########################################
#Some visual figure generation
pcSelect=20
Lun <- FindNeighbors(object = Lun, dims = 1:pcSelect)                #计算邻接距离
Lun <- FindClusters(object = Lun, resolution = 0.5)                  #对细胞分组,优化标准模块化
Lun <- RunTSNE(object = Lun, dims = 1:pcSelect)                      #TSNE聚类
pdf(file="06.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = Lun, pt.size = 1, label = TRUE)    #TSNE可视化
#完成聚类后，一定要记住保存数据，不然重新计算可要头疼了
saveRDS(Lun, file = "Lun.rds")
dev.off()
write.table(Lun$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)

##tSNE Plot
Lun <-RunTSNE(Lun, reduction = "harmony", dims = 1:20)  #等候大约

#可视化
pdf(file="06.umap.orig.pdf",width=10,height=6)
DimPlot(Lun, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
dev.off()
pdf(file="06.umap.orig1.pdf",width=10,height=6)
DimPlot(Lun, reduction = "umap", group.by = "orig.ident", pt.size = .1)
dev.off()

pdf(file="06.umap.phase.pdf",width=10,height=6)
DimPlot(Lun, reduction = "umap", group.by = "Phase", pt.size = .1)　#按照细胞周期进行划分
dev.off()

pdf(file="06.umap.pdf",width=10,height=6)
DimPlot(Lun, reduction = "umap", label = TRUE, pt.size = .1)　#注意作者在用同样参数设置后分为10个clusters，其实无关紧要，都需要通过marker重新贴现。
dev.off()

pdf(file="06.tsne.orig.pdf",width=10,height=6)
DimPlot(Lun, reduction = "tsne", group.by = "orig.ident", pt.size = 1, split.by = 'orig.ident')
dev.off()
pdf(file="06.tsne.orig1.pdf",width=10,height=6)
DimPlot(Lun, reduction = "tsne", group.by = "orig.ident", pt.size = 1)
dev.off()

pdf(file="06.tsne.phase.pdf",width=10,height=6)
DimPlot(Lun, reduction = "tsne", group.by = "Phase", pt.size = 1)　#按照细胞周期进行划分
dev.off()

pdf(file="06.tsne.pdf",width=10,height=6)
DimPlot(Lun, reduction = "tsne", label = TRUE, pt.size = 1)　#注意作者在用同样参数设置后分为10个clusters，其实无关紧要，都需要通过marker重新贴现。
dev.off()

#############################################07.注视细胞类型##################################################
setwd("C:\\Users\\muyumuge\\Desktop\\normal lung\\all normal")
rm(list = ls())  ## 魔幻操作，一键清空~
Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor
memory.limit(10000000) #改变R内存配置上限


library(SingleR) #查询R包的版本packageVersion("ExperimentHub")；查询所有包所在库的路径 ：.libPaths()
library(tibble)
library(dplyr)
library(ExperimentHub)

#新版本分型方法：singleR 1.1.14#
#因为我们刚刚从Seurat过来的，所以我们应该很想知道Seurat cluster 的细胞注释结果，因此，对Seurat的结果进行注释
#我们这里采用两个人类的参考集去做细胞注释，也可以选用这个参考集里的其他参考函数
hpca.se <- HumanPrimaryCellAtlasData() #ref.se <- HumanPrimaryCellAtlasData() 说明书给出的示例
bpe.se <- BlueprintEncodeData() #ref.se <- BlueprintEncodeData(rm.NA = "rows") 说明书给出的示例

#读入Seurat对象转换为SingleCell支持的对象
seurat.obj <- readRDS("../output/pbmc_tutorial.rds") #这一步因为Seurat对象还没有关闭，所以不需要导入，直接如下命名就可以
seurat.obj <-Lun
seurat.obj@meta.data$cell.type <- Idents(seurat.obj)
test <- as.SingleCellExperiment(seurat.obj)

#采用两个参考集一起进行注释
Anno <- SingleR(test = test,
                ref = list(HP = hpca.se , BP = bpe.se),
                labels = list(hpca.se$label.main , bpe.se$label.main),
                method = "cluster",
                cluster = test$cell.type)

#提取需要的细胞分类信息
Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)


#也可以将细胞注释信息重新添加到Seurat对象中去
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj, new.cluster.ids)



#旧版本分型方法：singleR1.0.0#
#提取需要备注的信息
counts<-Lun@assays$RNA@counts
clusters<-Lun@meta.data$seurat_clusters
ann=Lun@meta.data$orig.ident
#细胞自动注释
memory.limit(10000000) #这一步需要大内存支持分析运算，重新分配内存
singler = CreateSinglerObject(counts, annot = ann, "Lun", min.genes = 0,
                              species = "Human", citation = "",
                              ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
                              reduce.file.size = T, numCores = 1)
#提取细胞分型信息
singler$seurat = Lun
singler$meta.data$xy = Lun@reductions$tsne@cell.embeddings #reduction这里可更改tsne和umap，试试看哪个好
clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
#保存细胞分型信息
write.table(clusterAnn,file="07.clusterAnnTSNE.txt",quote=F,sep="\t",col.names=F)
write.table(singler$other,file="07.cellAnnTSNE.txt",quote=F,sep="\t",col.names=F)

#也可以将细胞注释信息重新添加到Seurat对象中去
new.cluster.ids <- clusterAnn
names(new.cluster.ids) <- levels(Lun)
Lun08 <- RenameIdents(Lun, new.cluster.ids) #重要认识：RenameIdents这个函数是更改Lun的active ident，即图片以这个作为分类依据


#标注后重新TSNE和UMAP及PCA作图
pdf(file="07.tsneTSNE.pdf",width=10,height=6)
TSNEPlot(object = Lun08, pt.size = 1, label = TRUE) 
dev.off()

pdf(file="07.umapTSNE.pdf",width=10,height=6)
UMAPPlot(object = Lun08, pt.size = 1, label = TRUE)
dev.off()




########################################08.提取某一类细胞进行后续分析###############################################
#Select a subset of Epithelial cells（近端小管），提取PT cells进行后续分析
PT <- SubsetData(Lun07, ident.use = "Epithelial cells", subset.raw = T) #等候大约2分钟，这里需要看文章注释细胞后选择哪个cluster来分析
saveRDS(PT,file="Epi.rds")



#install.packages(c("VGAM", "irlba", "matrixStats", "igraph", "combinat", "fastICA", "grid", "ggplot2", "reshape2", "plyr", "parallel", "methods"))

library(monocle)

#monocle输入文件类型有3种类型的格式：

#表达量文件：exprs基因在所有细胞中的count值矩阵。

#表型文件：phenoData。

#featureData

#3个特征文件转换成CellDataSet格式：


data <- as(as.matrix(PT@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = PT@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
my_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())


#对monocle对象进行归一化：
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1) ##过滤掉低质量的基因。等候大约5min。


#查看数据：
head(pData(my_cds))
head(fData(my_cds))

pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))
disp_table <- dispersionTable(my_cds)
head(disp_table)

#在进行降维聚类之前，先获得高表达的基因集作为每个聚类用来order的Feature基因。也可以使用所有的基因，但一些表达量特别低的
#基因提供的聚类信号往往会制造分析噪音，Feature基因的选择性很多，一种是可以根据基因的平均表达水平来进行筛选，另外我们也
#可以选择细胞间异常变异的基因。这些基因往往能较好地反映不同细胞的状态

table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)#细胞平均表达量大于0.1
my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)
pdf(file="08.dispersion.pdf",width=10,height=8)
plot_ordering_genes(my_cds)
dev.off()

expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10)) #细胞表达个数大于10
my_cds_subset <- my_cds[expressed_genes, ]
my_cds_subset
head(pData(my_cds_subset))
my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)
table(fData(my_cds_subset)$use_for_ordering)
pdf(file="08.variance.pdf",width=10,height=8)
plot_pc_variance_explained(my_cds_subset, return_all = FALSE)
dev.off()

#下面进行降维与聚类
my_cds_subset <- reduceDimension(my_cds_subset,max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)
my_cds_subset <- clusterCells(my_cds_subset, verbose = FALSE)
pdf(file="08.rho.pdf",width=10,height=8)
plot_rho_delta(my_cds_subset, rho_threshold = 2, delta_threshold = 10)
dev.off()

my_cds_subset <- clusterCells(my_cds_subset,rho_threshold = 2,delta_threshold = 10,skip_rho_sigma = T,verbose = FALSE)
table(pData(my_cds_subset)$Cluster)
pdf(file="08.tsneclusters.pdf",width=10,height=8)
plot_cell_clusters(my_cds_subset)
dev.off()

#选择定义细胞发展的基因
head(pData(my_cds_subset))
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~Cluster',cores = 22)
dim(clustering_DEG_genes)

#将细胞按照伪时间排序
library(dplyr)
clustering_DEG_genes %>% arrange(qval) %>% head()
my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree') #降维
my_cds_subset <- orderCells(my_cds_subset) #将细胞按照伪时间排序

#伪时序轨迹按不同方面绘图
pdf(file="08.cell_trajectory_seurat_clusters.pdf",width=10,height=8)
plot_cell_trajectory(my_cds_subset, color_by = "seurat_clusters")
dev.off()

pdf(file="08.cell_trajectory_RNA0.25.pdf",width=10,height=8)
plot_cell_trajectory(my_cds_subset, color_by = "RNA_snn_res.0.25")
dev.off()

pdf(file="08.cell_trajectory_RNA0.5.pdf",width=10,height=8)
plot_cell_trajectory(my_cds_subset, color_by = "RNA_snn_res.0.5")
dev.off()

pdf(file="08.cell_trajectory_old.ident.pdf",width=10,height=8)
plot_cell_trajectory(my_cds_subset, color_by = "old.ident")
dev.off()

pdf(file="08.cell_trajectory_Phase.pdf",width=10,height=8)
plot_cell_trajectory(my_cds_subset, color_by = "Phase")
dev.off()

head(pData(my_cds_subset))
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 22)
my_pseudotime_de %>% arrange(qval) %>% head()     #从这里挑选影响分化命运的top6基因
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
pdf(file="08.cell_trajectory_Pseudotime.pdf",width=10,height=8)
plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime")
dev.off()

#"A" stand for top 6 genes of affecting the fate decisions，影响fate decision的gene变化
A=c("RPS27","SFTPB","RPL9","AGER","RPS18","RPL10")
my_pseudotime_gene <-A
pdf(file="08.fate_decision_Top6gene.pdf",width=10,height=8)
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])
dev.off()

#Calculate the heat map of the top 50 genes
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$gene_short_name
pdf(file="08.pseudotime_heatmap.pdf",width=10,height=8)
my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],num_clusters = 5,cores = 22,show_rownames = TRUE,return_heatmap = TRUE)
dev.off()



























