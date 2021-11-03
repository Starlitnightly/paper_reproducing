# Nature复现：下游分析-Figure3b

在前面的教程里，我们已经得到了一个MOFA模型，很多分析都可以开始进行了，只不过下面的教程又将回到R语言中进行，在这里，我没用wsl上面的R，而是本地的R环境

## 1、多组学因子分析

### 1.1 加载依赖

```R
library(MOFA2)
library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)
```

### 1.2 基础设置

```R
## Define I/O ##

io <- list()
io$basedir <- "D:/drive/scnmt"
io$gene_metadata <- "D:/drive/scnmt/gene/Mmusculus_genes_BioMart.87.txt"
io$outdir <- "D:/drive/scnmt/result"

io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$rna.file <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$annos_dir  <- paste0(io$basedir, "/features/genomic_contexts")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")

## Define options ##
opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  # "genebody",
  "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)

opts$acc.annos <- c(
  # "genebody",
  "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)


# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  
  # E4.5
  "E4.5_Epiblast",
  
  # E5.5
  "E5.5_Epiblast",
  
  # E6.5
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Mesoderm",
  
  # E7.5
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 25      # minimum number of cells per feature
opts$met_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 25      # minimum number of cells per feature
opts$acc_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

# Define colors
opts$colors <- c(
  "Epiblast"="grey70",
  "Mesoderm"="#CD3278",
  "Primitive Streak"="sandybrown",
  "Endoderm"="#43CD80",
  "Ectoderm"="steelblue"
)

# window length for the overlap between genes and features
opts$overlapGenes  <- FALSE
opts$gene_window  <- 5e4

# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[stage_lineage%in%opts$stage_lineage]
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
opts$rna_cells <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells | id_acc %in% opts$acc_cells ]

```

### 1.3 模型加载

我们将上面生成的model_1模型放到了hdf5目录下

```R
model <- load_model(paste0(io$outdir,"/hdf5/model_1.hdf5"))
model
```

查看一下模型信息

```R
Trained MOFA with the following characteristics: 
 Number of views: 9 
 Views names: RNA Ectoderm enhancers (acc) Endoderm enhancers (acc) Mesoderm enhancers (acc) Acc CGI Promoters Ectoderm enhancers (met) Endoderm enhancers (met) Mesoderm enhancers (met) Met Promoters 
 Number of features (per view): 2500 1000 1000 1000 1000 1000 1000 1000 1000 
 Number of groups: 1 
 Groups names: gastrulation 
 Number of samples (per group): 2177 
 Number of factors: 9 
```

观察单细胞在多组学中的分布

```R
plot_data_overview(model)
```



<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/download.png" alt="download" style="zoom:50%;" />

### 1.4 模型重定义

在这里，我们需要对模型的一些参数跟变量进行重新定义

```R
sample_metadata_filt <- sample_metadata %>% 
  setkey(sample) %>% .[samples_names(model)] %>%
  .[,lineage10x_2:=stringr::str_replace_all(lineage10x_2,"_"," ")]
  
#重新定义组名
opts$views_names <- c(
  "acc_H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers (acc)",
  "acc_H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers (acc)",
  "acc_H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers (acc)",
  "acc_prom_2000_2000"="Acc CGI Promoters",
  "met_H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers (met)",
  "met_H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers (met)",
  "met_H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers (met)",
  "met_prom_2000_2000"="Met Promoters",
  "rna" = "RNA expression"
)

views_names(model) <- stringr::str_replace_all(views_names(model), opts$views_names)
```

### 1.5 计算多组学因子

我们计算MOFA模型中对于每一个样本的MOFA因子，在不同的组学数据中的方差比例

```R
r2 <- calculate_variance_explained_per_sample(model)$gastrulation
factors <- r2[,"RNA"]>0.005
#model <- subset_factors(model, which(factors))
factors_names(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")
```

因子可视化

```
plot_variance_explained(model)
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(2).png" alt="下载 (2)" style="zoom:50%;" />

可以看到，因子1，2，5，8可能是造成组学中样本数据变化的潜在因素

我们进一步观察这四个因子在不同细胞（样本）中的表现

```R
p <- plot_factor(model, 
  factors=c(1,2,5,8), 
  color_by=sample_metadata_filt$lineage10x_2,
  # shape_by=sample_metadata_filt$stage,
  dot_size = 1
)
p <- p + scale_color_manual(values=opts$colors)

# pdf(paste0(io$outdir,"/pdf/scatterplot_F1vsF2.pdf"), useDingbats = F, onefile = F, width=7, height=4)
print(p)
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(3).png" alt="下载 (3)" style="zoom:50%;" />

至此，我们可以确定因子1与Mesoderm（中胚层）及Ectoderm（外胚层）相关，因子2与Endoderm（内胚层）相关

## 2、细胞命运可视化

### 2.1 加载依赖

```R
library(umap)
library(Rtsne)
library(irlba)
```

### 2.2 定义函数

```R
theme_pub <- function() {
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.2)),
    axis.title = element_text(size=rel(1.2), color="black"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}
```

### 2.3 umap参数设置

```R
algorithms <- c("umap")

# umap.defaults$n_neighbors <- 20
# umap.defaults$min_dist <- 0.6

umap.defaults$n_neighbors <- 25
umap.defaults$min_dist <- 0.55
```

 ### 2.4 绘制不同种类的细胞分布图

```R
# Fetch factors
factors <- c(1,2,4,8)
Z <- get_factors(model,factors=c(1,2,4,8))$gastrulation 

for (algorithm in algorithms) {

  set.seed(42)
  if (algorithm=="tsne") {
    tsne <- Rtsne(Z, check_duplicates=FALSE, pca=FALSE, theta=0.5, dims=2)
    Z.out <- tsne$Y
  } else if (algorithm=="umap") {
    umap.out <- umap(Z, config = umap.defaults)
    Z.out <- umap.out$layout
  }
  
  # Flip a factor  
  Z.out[,2] <- -Z.out[,2]
  
  to.plot <- Z.out %>% as.data.table %>% .[,sample:=rownames(Z)] %>%
      merge(sample_metadata_filt, by="sample")

  p1 <- ggplot(to.plot, aes(x=V1, y=V2, color=lineage10x_2)) +
    geom_point(alpha=0.7, size=2.0) +
    scale_color_manual(values=opts$colors) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    theme_pub()
  
  # pdf(sprintf("%s/MOFA_%s.pdf",io$outdir,algorithm), width=5, height=7, useDingbats = F)
  print(p1)
  # dev.off()
  
  # Save coordinates
  # fwrite(to.plot[,c("sample","V1","V2")], sprintf("%s/%s_coordinates.txt",io$outdir,algorithm))
}
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(5).png" alt="下载 (5)" style="zoom:50%;" />

### 2.5 推断缺失值

我们使用impute推断在单细胞中目标值的缺失情况

```R
model <- impute(model)
```

### 2.6 绘制甲基化图谱

#### 2.6.1 外胚层增强子

```R
#
view <- "Ectoderm enhancers (met)"
factor <- 1

# Select top 50 enhancers with largest loading (Factor 1 for mesoderm enhancers, Factor 2 for endoderm enhancers)
tmp <- names(tail(sort(abs(get_weights(model, views=view, factor=factor)[[1]][,1])), n=50))
if (length(model@imputed_data)>0) {
  met <- colMeans(model@imputed_data[[view]]$gastrulation[tmp,], na.rm=T)
} else {
  met <- colMeans(model@data[[view]]$gastrulation[tmp,], na.rm=T)
}

# for better visualisation
met[met<((-6) )] <- (-6) 
met[met>(5)] <- 5

# Convert M-values to B-values
met <- 100*2**met/(1+2**met)

foo <- to.plot %>% merge(
  data.table(sample = samples_names(model)$gastrulation, met = met), by="sample"
)

p <- ggplot(foo, aes(x=V1, y=V2, color=met)) +
  geom_point(alpha=0.7, size=2.0) +
  scale_colour_gradientn(colours = brewer.pal(9, "OrRd")) +
  labs(x="UMAP Dimension 1", y=" UMAP Dimension 2") +
  theme_pub() 
# print(p)

# pdf(sprintf("%s/MOFA_%s_%d.pdf",io$outdir,view,factor), width=6, height=6, useDingbats = F)
print(p)
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(6).png" alt="下载 (6)" style="zoom:50%;" />

#### 2.6.2 中胚层增强子

```R
view <- "Mesoderm enhancers (met)"
factor <- 1

# Select top 50 enhancers with largest loading (Factor 1 for mesoderm enhancers, Factor 2 for endoderm enhancers)
tmp <- names(tail(sort(abs(get_weights(model, views=view, factor=factor)[[1]][,1])), n=50))
if (length(model@imputed_data)>0) {
  met <- colMeans(model@imputed_data[[view]]$gastrulation[tmp,], na.rm=T)
} else {
  met <- colMeans(model@data[[view]]$gastrulation[tmp,], na.rm=T)
}

# for better visualisation
met[met<((-6) )] <- (-6) 
met[met>(5)] <- 5

# Convert M-values to B-values
met <- 100*2**met/(1+2**met)

foo <- to.plot %>% merge(
  data.table(sample = samples_names(model)$gastrulation, met = met), by="sample"
)

p <- ggplot(foo, aes(x=V1, y=V2, color=met)) +
  geom_point(alpha=0.7, size=2.0) +
  scale_colour_gradientn(colours = brewer.pal(9, "OrRd")) +
  labs(x="UMAP Dimension 1", y=" UMAP Dimension 2") +
  theme_pub() 
# print(p)

# pdf(sprintf("%s/MOFA_%s_%d.pdf",io$outdir,view,factor), width=6, height=6, useDingbats = F)
print(p)
# dev.off()
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(7).png" alt="下载 (7)" style="zoom:50%;" />

#### 2.6.3 内胚层增强子

```R
#
view <- "Endoderm enhancers (met)"
factor <- 2

# Select top 50 enhancers with largest loading (Factor 1 for mesoderm enhancers, Factor 2 for endoderm enhancers)
tmp <- names(tail(sort(abs(get_weights(model, views=view, factor=factor)[[1]][,1])), n=50))
if (length(model@imputed_data)>0) {
  met <- colMeans(model@imputed_data[[view]]$gastrulation[tmp,], na.rm=T)
} else {
  met <- colMeans(model@data[[view]]$gastrulation[tmp,], na.rm=T)
}

# for better visualisation
met[met<((-6) )] <- (-6) 
met[met>(5)] <- 5

# Convert M-values to B-values
met <- 100*2**met/(1+2**met)

foo <- to.plot %>% merge(
  data.table(sample = samples_names(model)$gastrulation, met = met), by="sample"
)

p <- ggplot(foo, aes(x=V1, y=V2, color=met)) +
  geom_point(alpha=0.7, size=2.0) +
  scale_colour_gradientn(colours = brewer.pal(9, "OrRd")) +
  labs(x="UMAP Dimension 1", y=" UMAP Dimension 2") +
  theme_pub() 
# print(p)

# pdf(sprintf("%s/MOFA_%s_%d.pdf",io$outdir,view,factor), width=6, height=6, useDingbats = F)
print(p)
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(8).png" alt="下载 (8)" style="zoom:50%;" />

### 2.7 绘制染色质开发性图谱

#### 2.7.1 外胚层增强子

```R
#  "acc_H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers (acc)",
#"acc_H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers (acc)",
#  "acc_H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers (acc)",
#

view <- "Ectoderm enhancers (acc)"
factor <- 1

# Select top 50 enhancers with largest loading (Factor 1 for mesoderm enhancers, Factor 2 for endoderm enhancers)
tmp <- names(tail(sort(abs(get_weights(model, views=view, factor=factor)[[1]][,1])), n=50))
if (length(model@imputed_data)>0) {
  acc <- colMeans(model@imputed_data[[view]]$gastrulation[tmp,], na.rm=T)
} else {
  acc <- colMeans(model@data[[view]]$gastrulation[tmp,], na.rm=T)
}



# For better visualisation
acc[acc<(-1.5)] <- (-1.5)
acc[acc>(0.5)] <- 0.5

# Convert M-values to B-values
acc <- 100*2**acc/(1+2**acc)

foo <- to.plot %>% merge(
  data.table(sample = samples_names(model)$gastrulation, acc = acc), by="sample"
)

p <- ggplot(foo, aes(x=V1, y=V2, color=acc)) +
  geom_point(alpha=0.7, size=1.5) +
  scale_colour_gradientn(colours = rev(brewer.pal(9, "Blues"))) +
  labs(x="UMAP Dimension 1", y=" UMAP Dimension 2") +
  theme_pub() 

# pdf(sprintf("%s/MOFA_%s_%d.pdf",io$outdir,view,factor), width=6, height=6, useDingbats = F)
print(p)
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(9).png" alt="下载 (9)" style="zoom:50%;" />

#### 2.7.2 中胚层增强子

```R
#  "acc_H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers (acc)",
#"acc_H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers (acc)",
#  "acc_H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers (acc)",
#

view <- "Mesoderm enhancers (acc)"
factor <- 1

# Select top 50 enhancers with largest loading (Factor 1 for mesoderm enhancers, Factor 2 for endoderm enhancers)
tmp <- names(tail(sort(abs(get_weights(model, views=view, factor=factor)[[1]][,1])), n=50))
if (length(model@imputed_data)>0) {
  acc <- colMeans(model@imputed_data[[view]]$gastrulation[tmp,], na.rm=T)
} else {
  acc <- colMeans(model@data[[view]]$gastrulation[tmp,], na.rm=T)
}



# For better visualisation
acc[acc<(-1.5)] <- (-1.5)
acc[acc>(0.5)] <- 0.5

# Convert M-values to B-values
acc <- 100*2**acc/(1+2**acc)

foo <- to.plot %>% merge(
  data.table(sample = samples_names(model)$gastrulation, acc = acc), by="sample"
)

p <- ggplot(foo, aes(x=V1, y=V2, color=acc)) +
  geom_point(alpha=0.7, size=1.5) +
  scale_colour_gradientn(colours = rev(brewer.pal(9, "Blues"))) +
  labs(x="UMAP Dimension 1", y=" UMAP Dimension 2") +
  theme_pub() 

# pdf(sprintf("%s/MOFA_%s_%d.pdf",io$outdir,view,factor), width=6, height=6, useDingbats = F)
print(p)
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(10).png" alt="下载 (10)" style="zoom:50%;" />

#### 2.7.3 内胚层增强子

```R
#  "acc_H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers (acc)",
#"acc_H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers (acc)",
#  "acc_H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers (acc)",
#

view <- "Endoderm enhancers (acc)"
factor <- 2

# Select top 50 enhancers with largest loading (Factor 1 for mesoderm enhancers, Factor 2 for endoderm enhancers)
tmp <- names(tail(sort(abs(get_weights(model, views=view, factor=factor)[[1]][,1])), n=50))
if (length(model@imputed_data)>0) {
  acc <- colMeans(model@imputed_data[[view]]$gastrulation[tmp,], na.rm=T)
} else {
  acc <- colMeans(model@data[[view]]$gastrulation[tmp,], na.rm=T)
}



# For better visualisation
acc[acc<(-1.5)] <- (-1.5)
acc[acc>(0.5)] <- 0.5

# Convert M-values to B-values
acc <- 100*2**acc/(1+2**acc)

foo <- to.plot %>% merge(
  data.table(sample = samples_names(model)$gastrulation, acc = acc), by="sample"
)

p <- ggplot(foo, aes(x=V1, y=V2, color=acc)) +
  geom_point(alpha=0.7, size=1.5) +
  scale_colour_gradientn(colours = rev(brewer.pal(9, "Blues"))) +
  labs(x="UMAP Dimension 1", y=" UMAP Dimension 2") +
  theme_pub() 

# pdf(sprintf("%s/MOFA_%s_%d.pdf",io$outdir,view,factor), width=6, height=6, useDingbats = F)
print(p)
```

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/%E4%B8%8B%E8%BD%BD%20(11).png" alt="下载 (11)" style="zoom:50%;" />



## 3、后话

到这里，Figure3b就完全复现完了，中间还是有很多很有意思的点。由于我这里只有一组，所以MOFA计算出来的因子其实是不同组学数据可变性的来源，而并不是组间差异，进而借助组学数据的差异来确定

注：本研究运用MOFA分析的原因是，每个组学层的数据都有一大堆细胞，作者的目的是找出导致这一大堆细胞不同的原因，也就是所谓的方差最大。MOFA应用在找有病没病的人的hub基因可能会有些麻烦，因为可能hub基因在有病的人体内表达都很高，在没病的人体内表达都很低，这种情况下，组内组学数据的方差就会很小。不过也不是不可以，把有病没病的数据混成一组就好了。