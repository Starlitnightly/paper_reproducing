# Nature复现：MOFA模型-Figure3b

呃，其实最开始是想复现Figure3c的来着，等做到最后一步的时候才发现，作者在github上的标注出错了，没办法，只好将错就错了

## 1、数据下载

首先我们需要把作者的处理后的数据给下载下来，总共33g，当然分析的时候不会一次载入33g的数据量，只是其中一小部分，我将下载好的内容放到文件夹`scnmt`下

```
ftp://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/scnmt_gastrulation.tar.gz
```

作者其实每一幅图都给了代码，但是在这个教程里面，不会完全照搬里面的，因为MOFA这个包的缘故，你大概率，是运行不成功的，但这里还是放出作者代码的仓库，建议一并下载下来，可以供参考

```
https://github.com/rargelaguet/scnmt_gastrulation
```

## 2、MOFA模型数据准备

### 2.1 环境准备

我们在刚刚安装好的wsl系统里面，输入jupyter lab，打开R环境，可以看到右上角有个R的标注

![image-20211006231647297](https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/image-20211006231647297.png)

### 2.2 载入R包

接下来就是正式的跑代码环节，为什么是wsl呢，因为有一些代码我暂时还不知道怎么修改，只能保持跟作者一致的代码

以下#代表源注释，##代表我的注释

```R
suppressMessages(library(MOFA2))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(scater))
suppressMessages(library(reticulate))
suppressMessages(library(argparse))
```

### 2.3 配置I/O

在这一步，只需要注意`io$basedir <- "/mnt/d/drive/scnmt"`的内容就好了，把这个改成你上面数据下载下来的安装路径

```R
## Define I/O ##

io$basedir <- "/mnt/d/drive/scnmt"
io$gene_metadata <- "/mnt/d/drive/scnmt/gene/Mmusculus_genes_BioMart.87.txt"
io$outdir <- "/mnt/d/drive/scnmt/result"

io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$rna.file <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$annos_dir  <- paste0(io$basedir, "/features/genomic_contexts")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")
```

### 2.4 定义变量

在这里，我们会把本次数据分析所需要的变量以及环境都定义清楚

```R
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
##从sample.metadata里面读取了样本数据，并且做了以下处理
##1、定义了新列stage_lineage，类型为factor，值为stage与lineage10x_2列的值拼接而成
##2、移除了不在opts$stage_lineage里的stage_lineage的行
tmp <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[stage_lineage%in%opts$stage_lineage]
##选择了pass_metQC等于True的行，并单独取出id_met列存进met_cell列表里面
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
##选择了pass_rnaQC等于True的行，并单独取出id_rna列
opts$rna_cells <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]
```

### 2.5 定义单细胞信息

每个单细胞的特征作者已经处理好了，我们这里只需要直接导入就好

```R
##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells | id_acc %in% opts$acc_cells ]
```

### 2.6 定义矩阵处理函数

```R
matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}
```

### 2.7 加载数据

这一步要用到linux命令zcat，我不是很会改成windows下的命令，不然我就直接改了，也就不用安wsl了

```R
# Load Methylation data
met_dt <- lapply(opts$met.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir,n), showProgress=F, stringsAsFactors=F, quote="") %>% .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Load Accessibility data
acc_dt <- lapply(opts$acc.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[V1%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","N","rate"))

# Load RNA data
sce <- readRDS(io$rna.file) %>% .[,opts$rna_cells]

# Load annotation metadata
feature_metadata <- lapply(unique(c(opts$met.annos,opts$acc.annos)), function(i) 
  fread(sprintf("%s/%s.bed",io$annos_dir,i), stringsAsFactors=T)[,c(1,2,3,4,5,6)]) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno"))


# Load gene metadata 
gene_metadata <- fread(io$gene_metadata,stringsAsFactors=T) %>% 
  setnames(c("ens_id","symbol"),c("id","gene")) %>% 
  .[,chr:=as.factor(sub("chr","",chr))]
```

### 2.8 数据解析

作者对数据都进行了简单的解析

```R
# Parse gene and feature metadata
feature_metadata_filt.met <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% met_dt[anno==y,id]] ) %>%
  rbindlist
feature_metadata_filt.acc <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% acc_dt[anno==y,id]] ) %>%
  rbindlist
gene_metadata_filt <- gene_metadata %>% .[,c("chr","start","end","gene")] %>% 
  .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
  setkey(chr,start,end)
  
## Parse RNA expression data ##

# Convert to data.table
rna_dt <- exprs(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))

## Parse accessibility data ##
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value

## Parse methylation data ##
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value

##############################
## Merge data with metadata ##
##############################

acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage","stage_lineage")], by="id_acc") %>% droplevels()
met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage","stage_lineage")], by="id_met") %>% droplevels()
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage","stage_lineage")], by="id_rna") %>% droplevels()

```

### 2.9 将基因组特征与重叠基因联系起来

```R
# Methylation
if (opts$overlapGenes) {
  met_list <- list()
  for (i in unique(met_dt$anno)){
    
    # Subset corresponding anno
    met_tmp <- met_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSMUSG", unique(met_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.met[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      met_list[[i]] <- merge(met_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSMUSG", unique(met_tmp$id)) == TRUE)) {
      met_list[[i]] <- merge(met_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  met_dt <- rbindlist(met_list)
  rm(met_list, met_tmp,feature_metadata_tmp,ov)
} else {
  met_dt[,gene:="NA"]
}

# Accessibility
if (opts$overlapGenes) {
  acc_list <- list()
  for (i in unique(acc_dt$anno)){
    
    # Subset corresponding anno
    acc_tmp <- acc_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.acc[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      acc_list[[i]] <- merge(acc_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == TRUE)) {
      acc_list[[i]] <- merge(acc_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  acc_dt <- rbindlist(acc_list)
  rm(acc_list, acc_tmp,feature_metadata_tmp,ov)
} else {
  acc_dt[,gene:="NA"]
}


```

### 2.10 数据过滤

```R
#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by  minimum number of cells
met_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$met_min.cells] %>% .[,N:=NULL]

# Filter features by variance
met_dt <- met_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
acc_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$acc_min.cells] %>% .[,N:=NULL]

# Filter features by variance
acc_dt <- acc_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

################################
## Filter RNA expression data ##
################################

# Remove lowly expressed genes
rna_dt <- rna_dt[,mean:=mean(expr),by="ens_id"] %>% .[mean>0.1] %>% .[,mean:=NULL]

# Remove genes with constant expression levels
rna_dt <- rna_dt[,var:=var(expr),by="ens_id"] %>% .[var>0.1] %>% .[,var:=NULL]

# Filter genes with low cellular detection rate and sites with low coverage across samples
rna_dt <- rna_dt[,cdr:=sum(expr>0)/length(opts$rna_cells), by="ens_id"] %>% .[cdr>=opts$rna_min.cdr] %>% .[,cdr:=NULL]

# Extract top N highly variable genes
rna_dt <- rna_dt[,var:=var(expr), by="ens_id"] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

```

### 2.11 协变量回归

作者使用了模型拟合的残差作为基因表达的标准化定量。

```R
# RNA: number of expressed genes
foo <- data.table(id_rna=colnames(sce), covariate=sce$total_features_by_counts/nrow(sce))
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~covariate)[["residuals"]], by=c("gene")] %>%
  .[,covariate:=NULL]
  
# RNA: strong batch effect between E7.5 plates
if (any(sample_metadata$stage=="E7.5")) {
  foo <- sample_metadata[stage=="E7.5",c("id_rna","plate")] %>%
    .[,plate:=as.factor(grepl("PS_VE",plate))]
  rna_dt <- rbind(
    rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
      .[,expr:=lm(formula=expr~plate)[["residuals"]], by=c("gene")] %>% .[,c("plate"):=NULL],
    rna_dt[!id_rna%in%foo$id_rna]
  )
}

# RNA: strong batch effect between E6.5 plates
if (any(sample_metadata$stage=="E6.5")) {
  foo <- sample_metadata[stage=="E6.5" & plate%in%c("E6.5_late_Plate1","E6.5_late_Plate2"),c("id_rna","plate")]
  rna_dt <- rbind(
    rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
      .[,expr:=lm(formula=expr~plate)[["residuals"]], by=c("gene")] %>% .[,plate:=NULL],
    rna_dt[!id_rna%in%foo$id_rna]
  )
}

# Methylation: differences in mean methylation rate
foo <- met_dt[,.(covariate=mean(m)),by=c("id_met")]
foo <- fread(io$met.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_met","mean")]
met_dt <- met_dt %>% merge(foo, by="id_met") %>%
  .[,m:=mean(m) + lm(formula=m~mean)[["residuals"]], by=c("id","anno")]


# Accessibility: differences in global accessibility rate
foo <- fread(io$acc.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_acc","mean")]
acc_dt <- acc_dt %>% merge(foo, by="id_acc") %>%
  .[,m:=mean(m) + lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

```

### 2.12 特征选择

在这里，作者选择了N个最大的样本方差最大的值，比如选择了前N个基因表达方差最大的基因

```R
# RNA: Extract top N highly variable genes
keep_hv_genes <- rna_dt[,.(var=var(expr)), by="ens_id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id
rna_dt <- rna_dt[ens_id%in%as.character(keep_hv_genes)] %>% droplevels()

# Accessibility: Extract top N most variable features
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

# Methylation: Extract top N most variable features
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

met_dt[,var(m),by=c("id","anno")]
```

### 2.13 合并数据

在这里，作者将三个表合并成了一个data.txt用于模型分析，我们在下一步的MOFA建模中会用到这个data.txt

```R
data1 <- rna_dt %>% .[,c("sample","gene","expr")] %>%  
  setnames(c("sample","feature","value")) %>% .[,c("feature_group","sample_group"):=list("RNA","gastrulation")]
data2 <- met_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("met_",feature), paste0("met_",feature_group), "gastrulation")]
data3 <- acc_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("acc_",feature), paste0("acc_",feature_group), "gastrulation")]

data <- rbind(data1,data2,data3)

outfile <- paste0(io$outdir,"/data.txt")
fwrite(data, file=outfile, col.names=T, quote=F, sep="\t")
```

作者后续的代码可能MOFA包跑得动，由于我是MOFA2包，所以后面的我就跑不动了，开始切环境训练模型

## 3、MOFA模型训练

该模型训练我用上了GPU，所以切到了colab上

### 3.1 依赖准备

mofapy2这个包有点神奇，我直接在jupyter里装会报错，所以我打开了终端，输入

```
pip install mofapy2
```

安装完后，我们导入包

```python
from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np

# initialise the entry point
ent = entry_point()
```

出现以下图案就是成功

```python
        #########################################################
        ###           __  __  ____  ______                    ### 
        ###          |  \/  |/ __ \|  ____/\    _             ### 
        ###          | \  / | |  | | |__ /  \ _| |_           ### 
        ###          | |\/| | |  | |  __/ /\ \_   _|          ###
        ###          | |  | | |__| | | / ____ \|_|            ###
        ###          |_|  |_|\____/|_|/_/    \_\              ###
        ###                                                   ### 
        ######################################################### 
```

### 3.2 数据清洗

我们前面得到了data.txt这个数据文件，于是我们需要稍微改一下columns，使其符合mofapy2的要求

```python
data=pd.read_csv('/content/drive/MyDrive/data/data.txt',sep='\t')
data.columns=['sample','feature','value','view','group']
data.head()
```

成功返回了如下内容

```python
	sample	feature	value	view	group
0	E6.75_Plate1_E1	Cdc45	-2.005848	RNA	gastrulation
1	E6.75_Plate1_E1	Narf	-1.188268	RNA	gastrulation
2	E6.75_Plate1_E1	Klf6	4.435413	RNA	gastrulation
3	E6.75_Plate1_E1	Wnt3	-0.727745	RNA	gastrulation
4	E6.75_Plate1_E1	Fer	-0.995191	RNA	gastrulation
```

### 3.3 模型数据设置

- **scale_groups**: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is False
- **scale_views**: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is False

```python
ent.set_data_options(
    scale_groups = False, 
    scale_views = False
)
```

### 3.4 导入数据到模型中

- **likelihoods**: a list of strings, either "gaussian", "poisson" or "bernoulli". If None (default), they are guessed internally

在这里，我们认为每一个组学数据都是gaussian分布

```python
ent.set_data_df(data, likelihoods = ["gaussian","gaussian","gaussian","gaussian","gaussian","gaussian","gaussian","gaussian","gaussian"])
```

程序输出结果：

```
Loaded group='gastrulation' view='RNA' with N=2147 samples and D=2500 features...
Loaded group='gastrulation' view='acc_H3K27ac_distal_E7.5_Ect_intersect12' with N=729 samples and D=1000 features...
Loaded group='gastrulation' view='acc_H3K27ac_distal_E7.5_End_intersect12' with N=729 samples and D=1000 features...
Loaded group='gastrulation' view='acc_H3K27ac_distal_E7.5_Mes_intersect12' with N=729 samples and D=1000 features...
Loaded group='gastrulation' view='acc_prom_2000_2000' with N=729 samples and D=1000 features...
Loaded group='gastrulation' view='met_H3K27ac_distal_E7.5_Ect_intersect12' with N=826 samples and D=1000 features...
Loaded group='gastrulation' view='met_H3K27ac_distal_E7.5_End_intersect12' with N=826 samples and D=1000 features...
Loaded group='gastrulation' view='met_H3K27ac_distal_E7.5_Mes_intersect12' with N=826 samples and D=1000 features...
Loaded group='gastrulation' view='met_prom_2000_2000' with N=826 samples and D=1000 features...
```

### 3.5 模型参数设置

- **factors**: number of factors
- **spikeslab_weights**: use spike-slab sparsity prior in the weights? default is TRUE
- **ard_factors**: use ARD prior in the factors? Default is TRUE if using multiple groups. This is guessed by default.
- **ard_weights**: use ARD prior in the weights? Default is TRUE if using multiple views. This is guessed by default.

```python
ent.set_model_options(
    factors = 10, 
    spikeslab_weights = True, 
    ard_factors = True,
    ard_weights = True
)
```

程序输出结果

```
Model options:
- Automatic Relevance Determination prior on the factors: True
- Automatic Relevance Determination prior on the weights: True
- Spike-and-slab prior on the factors: False
- Spike-and-slab prior on the weights: True
Likelihoods:
- View 0 (RNA): gaussian
- View 1 (acc_H3K27ac_distal_E7.5_Ect_intersect12): gaussian
- View 2 (acc_H3K27ac_distal_E7.5_End_intersect12): gaussian
- View 3 (acc_H3K27ac_distal_E7.5_Mes_intersect12): gaussian
- View 4 (acc_prom_2000_2000): gaussian
- View 5 (met_H3K27ac_distal_E7.5_Ect_intersect12): gaussian
- View 6 (met_H3K27ac_distal_E7.5_End_intersect12): gaussian
- View 7 (met_H3K27ac_distal_E7.5_Mes_intersect12): gaussian
- View 8 (met_prom_2000_2000): gaussian

```

### 3.6 模型训练设置

- **iter**: number of iterations. Default is 1000.
- **convergence_mode**: "fast", "medium", "slow". For exploration, the fast mode is good enough.
- **startELBO**: initial iteration to compute the ELBO (the objective function used to assess convergence)
- **freqELBO**: frequency of computations of the ELBO (the objective function used to assess convergence)
- **dropR2**: minimum variance explained criteria to drop factors while training gpu_mode: use GPU mode? (needs cupy installed and a functional GPU, see https://cupy.chainer.org/)
- **verbose**: verbose mode?
- **seed**: random seed

```python
ent.set_train_options(
    iter = 5000, 
    convergence_mode = "slow", 
    startELBO = 1, 
    freqELBO = 1, 
    dropR2 = 0.001, 
    gpu_mode = True, 
    verbose = False, 
    seed = 42
)
```

### 3.7 （可选）随机推理选项

如果样本数量非常大（以>1e4 的顺序），您可能想尝试随机推理方案。 然而，它需要一些额外的超参数，在某些数据集中可能需要优化（参见随机小插图），实际上都是机器学习的参数

- **batch_size**: float value indicating the batch size (as a fraction of the total data set: either 0.10, 0.25 or 0.50)
- **learning_rate**: learning rate (we recommend values from 0.25 to 0.50)
- **forgetting_rate**: forgetting rate (we recommend values from 0.75 to 1.0)

```python
# We do not want to use stochastic inference for this data set

# ent.set_stochastic_options(
#     batch_size = 0.5,
#     learning_rate = 0.5, 
#     forgetting_rate = 0.25
# )
```

### 3.8 训练MOFA模型

终于到这一步了，还是比较简单的

```python
ent.build()

ent.run()
```

```
######################################
## Training the model with seed 42 ##
######################################


ELBO before training: -80858378.90 

Iteration 1: time=8.77, ELBO=-17216582.96, deltaELBO=63641795.939 (78.70773173%), Factors=9
Iteration 2: time=3.07, ELBO=-17158861.17, deltaELBO=57721.785 (0.07138627%), Factors=9
Iteration 3: time=3.65, ELBO=-17147286.77, deltaELBO=11574.401 (0.01431441%), Factors=9
......
Iteration 156: time=3.25, ELBO=-17055635.59, deltaELBO=4.499 (0.00000556%), Factors=9
Iteration 157: time=3.22, ELBO=-17055631.25, deltaELBO=4.337 (0.00000536%), Factors=9
Iteration 158: time=3.23, ELBO=-17055627.07, deltaELBO=4.189 (0.00000518%), Factors=9
Iteration 159: time=3.34, ELBO=-17055623.01, deltaELBO=4.053 (0.00000501%), Factors=9
Iteration 160: time=3.34, ELBO=-17055619.08, deltaELBO=3.930 (0.00000486%), Factors=9
Iteration 161: time=3.36, ELBO=-17055615.26, deltaELBO=3.819 (0.00000472%), Factors=9

Converged!



#######################
## Training finished ##
#######################
```

### 3.9 模型的保存

```
ent.save(outfile='model_1.hdf5')
```

到此，我们成功通过多组学数据成功训练出了一个MOFA模型，不过进一步的分析要到R语言中进行了