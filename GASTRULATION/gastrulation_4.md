# Nature复现：增强子分析-Figure3c

前面的教程主要是单细胞分析，在本小节中，我们将复现作者的Figure3c，增强子部分，增强子是DNA上一小段可与蛋白质结合的区域，与蛋白质结合之后，基因的转录作用将会加强。增强子可能位于基因上游，也可能位于下游。作者的目的是看不同的胚胎细胞中，他们增强子的变化是怎样的。

## 1、数据预处理

### 1.1 加载依赖

```R
library(data.table)
library(purrr)
library(ggplot2)
```

### 1.2 路径设置

```R
io <- list()
io$basedir <- "/mnt/g/共享云端硬盘/794/scnmt_data"


# Sample metadata
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")

# Genomic contects
io$annos_dir <- paste0(io$basedir,"/features/genomic_contexts")

# DNA methylation and chromatin accessibility data
io$met.dir <- paste0(io$basedir,"/met/cpg_level")
io$acc.dir <- paste0(io$basedir,"/acc/gpc_level")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")

# Folders with the differential analysis results
io$diff.met <- paste0(io$basedir,"/met/results/differential")
io$diff.acc <- paste0(io$basedir,"/acc/results/differential")

# Output directory
io$pdfdir <- paste0(io$basedir,"/metacc/pseudobulked_profiles/lineage_enhancers")
```

### 1.3 操作设置

```R
opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Mesoderm",
  "E7.5_Endoderm"
  
)

# Define genomic contexts
opts$annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
)

# Define window positions and characteristics
opts$positions <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="center",
  "H3K27ac_distal_E7.5_End_intersect12"="center",
  "H3K27ac_distal_E7.5_Ect_intersect12"="center"
)
opts$window_size <- 2000
opts$met.tile <- 200
opts$acc.tile <- 150

# Thresholds in the differential analysis to subset lineage-defining enhancers
opts$diff.type <- 2
opts$min.fdr <- 0.10
opts$min.met.diff <- 5
opts$min.acc.diff <- 5

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% 
  .[stage_lineage%in%opts$stage_lineage] 
opts$met.cells <- tmp %>% .[pass_metQC==T,id_met]
opts$acc.cells <- tmp %>% .[pass_accQC==T,id_acc]

rm(tmp)
```

我们还需要配置一些RNA的设置

```R
# Define genomic contexts for methylation
opts$met.annos <- c(
  # "genebody"="Gene body",
  "prom_2000_2000"="Promoters",
  # "prom_2000_2000_cgi"="CGI promoters",
  # "prom_2000_2000_noncgi"="non-CGI promoters",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
  # "H3K4me3_E7.5_Mes"="Mes- H3K4me3",
  # "H3K4me3_E7.5_End"="End- H3K4me3",
  # "H3K4me3_E7.5_Ect"="Ect- H3K4me3"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
  # "genebody"="Gene body",
  "prom_2000_2000"="Promoters",
  # "prom_2000_2000_cgi"="CGI promoters",
  # "prom_2000_2000_noncgi"="non-CGI promoters",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
  # "H3K4me3_E7.5_Mes"="Mes- H3K4me3",
  # "H3K4me3_E7.5_End"="End- H3K4me3",
  # "H3K4me3_E7.5_Ect"="Ect- H3K4me3"
)

opts$annos <- opts$met.annos <- opts$acc.annos

# Overlap differential sites with nearby genes
opts$overlapGenes <- FALSE

# window length for the overlap between genes and features
opts$gene_window <- 25000

# How to select differential hits?
#   Option 1 (more liberal): (lineage_A) vs (lineage_B,lineage_C)
#   Option 2 (more conservative): (lineage_A vs lineage_B) AND (lineageA vs lineage_C)
opts$diff.type <- 2
opts$min.fdr <- 0.10
opts$min.acc.diff <- 5
opts$min.met.diff <- 5

# Lineage colors
opts$colors <- c(
  Ectoderm = "steelblue", 
  Endoderm = "#43CD80", 
  Mesoderm = "violetred"
)

######################
## Define functions ##
######################

gg_barplot <- function(tmp, title = "", ylim=NULL) {
  
  if (is.null(ylim)) {
    ylim <- c(min(tmp$value, na.rm=T), max(tmp$value, na.rm=T))
  }
  
  p <- ggplot(tmp, aes(x=anno, y=value)) +
    geom_bar(aes(fill=assay), color="black", stat="identity", position="dodge", size=0.25) +
    scale_fill_manual(values=c("met"="#F37A71", "acc"="#00BFC4")) +
    geom_hline(yintercept=0, color="black") +
    scale_y_continuous(limits=c(ylim[1],ylim[2])) +
    labs(title=i, x="", y="Number of hits") +
    theme_bw() +
    theme(
      plot.title = element_text(size=11, face='bold', hjust=0.5),
      axis.text = element_text(size=rel(1.0), color='black'),
      axis.text.x = element_text(size=rel(1.0), angle=60, hjust=1, vjust=1, color="black"),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size=rel(1.0), color='black'),
      axis.line = element_line(color="black"),
      legend.position="none"
    )
  
  return(p)
}

theme_pub <- function() {
  theme_bw() +
  theme(
    axis.text.x = element_text(size=rel(1.2), angle=60, hjust=1, vjust=1, color="black"),
    axis.text.y = element_text(size=rel(1.2), color="black"),
    axis.title.y = element_text(size=rel(1.2), color="black"),
    legend.position = "right"
    )
}
```

### 1.4 加载样本数据

```R
sample_metadata <- fread(io$sample.metadata) %>%
  .[,c("sample","id_acc","id_met","id_rna","stage","lineage10x_2")] %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep=" ")] %>%
  .[id_met%in%opts$met.cells | id_acc%in%opts$acc.cells]
  
sample_metadata %>% 
  .[,stage_lineage:=ifelse(stage_lineage=="E5.5 Epiblast","E6.5 Epiblast",stage_lineage)]
```

### 1.5 加载基因注释

```R
# Load genomic annotations
anno_list <- list()
for (anno in names(opts$annos)) {
  tmp <- fread(sprintf("%s/%s.bed",io$annos_dir,anno))[,c(1,2,3,4,5,6)]
  colnames(tmp) <- c("chr","start","end","strand","id","anno")
  
  # Define central position for the window approach
  if (opts$positions[anno] == "start") {
    tmp <- rbind(tmp[strand=="+",.(chr,start,strand,id,anno)] %>% .[,center:=start] %>% .[,c("start"):=NULL], 
                 tmp[strand=="-",.(chr,end,strand,id,anno)] %>% .[,center:=end] %>% .[,c("end"):=NULL]) 
  }
  if (opts$positions[anno] == "center") {
    stopifnot(all(tmp[,end] > tmp[,start]))
    tmp <- tmp[,.(chr,start,end,strand,id,anno)][,center:=round(end+start)/2][,c("start","end"):=NULL]
  }
  if (opts$positions[anno] == "end") {
    tmp <- rbind(tmp[strand=="+",.(chr,end,strand,id,anno)][,center:=end][,c("end"):=NULL], 
                 tmp[strand=="-",.(chr,start,strand,id,anno)][,center:=start][,c("start"):=NULL])
  }
  anno_list[[anno]] <- tmp %>% .[, c("start","end") := list(center-opts$window_size,center+opts$window_size)]
}

anno_df <- rbindlist(anno_list) %>% 
  .[,chr:=sub("chr","",chr)] %>%
  setkey(chr,start,end)

rm(anno_list)
```

### 1.6 加载RNA差异表达分析结果

```
if (!is.null(io$diff.rna)) {
  if (opts$diff.type==1) {
    
    # Mesoderm-specific
    diff.rna.mes <- fread(sprintf("%s/E7.5Mesoderm_vs_E7.5EndodermEctoderm.txt.gz",io$diff.rna)) %>%
      .[,lineage:="Mesoderm"] %>% .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff]
    
    # Ectoderm-specific
    diff.rna.ect <- fread(sprintf("%s/E7.5Ectoderm_vs_E7.5MesodermEndoderm.txt.gz",io$diff.rna)) %>%
      .[,lineage:="Ectoderm"] %>% .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff]
    
    # Endoderm-specific
    diff.rna.end <- fread(sprintf("%s/E7.5Endoderm_vs_E7.5MesodermEctoderm.txt.gz",io$diff.rna)) %>%
      .[,lineage:="Endoderm"] %>% .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff]
    
  } else if (opts$diff.type%in%c(2,3)) {
    
    # Mesoderm-specific
    diff.rna.mes <- rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Ectoderm"]
      ) %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff] %>% 
      data.table::dcast(id+symbol+lineage1~lineage2, value.var=c("logFC","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(logFC_Ectoderm)==sign(logFC_Endoderm)] %>%
      .[,logFC:=(logFC_Endoderm+logFC_Ectoderm)/2] %>%
      .[,c("id","symbol","lineage1","logFC","sig")] %>% setnames("lineage1","lineage")
    
    
    # Ectoderm-specific
    diff.rna.ect <- rbind(
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Endoderm"],
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Mesoderm"]
    ) %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff] %>% 
      data.table::dcast(id+symbol+lineage1~lineage2, value.var=c("logFC","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Endoderm==T & sign(logFC_Endoderm)==sign(logFC_Mesoderm)] %>%
      .[,logFC:=(logFC_Endoderm+logFC_Mesoderm)/2] %>%
      .[,c("id","symbol","lineage1","logFC","sig")] %>% setnames("lineage1","lineage")
    
    # Endoderm-specific
    diff.rna.end <- rbind(
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Ectoderm"],
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm.txt.gz",io$diff.rna)) %>% .[,lineage2:="Mesoderm"]
    ) %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(logFC)>opts$min.rna.diff] %>% 
      data.table::dcast(id+symbol+lineage1~lineage2, value.var=c("logFC","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(logFC_Ectoderm)==sign(logFC_Mesoderm)] %>%
      .[,logFC:=(logFC_Ectoderm+logFC_Mesoderm)/2] %>%
      .[,c("id","symbol","lineage1","logFC","sig")] %>% setnames("lineage1","lineage")
  }
  
  diff.rna <- do.call("rbind", list(diff.rna.mes,diff.rna.end,diff.rna.ect))
  rm(diff.rna.mes,diff.rna.ect,diff.rna.end)
}
```

### 1.7 加载Met差异表达分析结果

```R
if (!is.null(io$diff.met)) {
  
  # Mesoderm-specific
  if (opts$diff.type==1) {
    diff.met.mes <- lapply(names(opts$met.annos), function(x)
    # diff.met.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
      fread(sprintf("%s/E7.5Mesoderm_vs_E7.5EndodermEctoderm_%s.txt.gz",io$diff.met,x))
    ) %>% rbindlist() %>% .[,lineage:="Mesoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff]
  
    # Ectoderm-specific
    diff.met.ect <- lapply(names(opts$met.annos), function(x)
    # diff.met.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5MesodermEndoderm_%s.txt.gz",io$diff.met,x))
    ) %>% rbindlist() %>% .[,lineage:="Ectoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff]
  
    # Endoderm-specific
    diff.met.end <- lapply(names(opts$met.annos), function(x)
    # diff.met.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5MesodermEctoderm_%s.txt.gz",io$diff.met,x))
    ) %>% rbindlist() %>% .[,lineage:="Endoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff]
  }
  
  if (opts$diff.type==2) {
    
    # Mesoderm-specific
    diff.met.mes <- lapply(names(opts$met.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(diff_Ectoderm)==sign(diff_Endoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Ectoderm-specific
    diff.met.ect <- lapply(names(opts$met.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Mesoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Mesoderm==T & sign(diff_Endoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Mesoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Endoderm-specific
    diff.met.end <- lapply(names(opts$met.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Mesoderm"],
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(diff_Ectoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Mesoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
    }
    if (opts$diff.type==3) {
    
    # Mesoderm-specific
    diff.met.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(diff_Ectoderm)==sign(diff_Endoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Ectoderm-specific
    diff.met.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Mesoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Mesoderm==T & sign(diff_Endoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Mesoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Endoderm-specific
    diff.met.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Mesoderm"],
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.met,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.met.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(diff_Ectoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Mesoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
    }
  
  diff.met <- do.call("rbind", list(diff.met.mes,diff.met.end,diff.met.ect))
  rm(diff.met.mes,diff.met.ect,diff.met.end)
}

```

### 1.8 加载Acc差异表达分析结果

```R
if (!is.null(io$diff.acc)) {
  
  if (opts$diff.type==1) {
    
    # Mesoderm-specific
    diff.acc.mes <- lapply(names(opts$acc.annos), function(x)
      fread(sprintf("%s/E7.5Mesoderm_vs_E7.5EndodermEctoderm_%s.txt.gz",io$diff.acc,x))
    ) %>% rbindlist() %>% .[,lineage:="Mesoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff]
  
    # Ectoderm-specific
    diff.acc.ect <- lapply(names(opts$acc.annos), function(x)
      fread(sprintf("%s/E7.5Ectoderm_vs_E7.5MesodermEndoderm_%s.txt.gz",io$diff.acc,x))
    ) %>% rbindlist() %>% .[,lineage:="Ectoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff]
  
    # Endoderm-specific
    diff.acc.end <- lapply(names(opts$acc.annos), function(x)
      fread(sprintf("%s/E7.5Endoderm_vs_E7.5MesodermEctoderm_%s.txt.gz",io$diff.acc,x))
    ) %>% rbindlist() %>% .[,lineage:="Endoderm"] %>%
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff]
  
  }
  
  if (opts$diff.type==2) {
    # Mesoderm-specific
    diff.acc.mes <- lapply(names(opts$acc.annos), function(x)
    # diff.acc.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(diff_Endoderm)==sign(diff_Ectoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Ectoderm-specific
    # diff.acc.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
    diff.acc.ect <- lapply(names(opts$acc.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Mesoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Mesoderm==T & sign(diff_Endoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Mesoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Endoderm-specific
    # diff.acc.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
    diff.acc.end <- lapply(names(opts$acc.annos), function(x)
      rbind(
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Mesoderm"],
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(diff_Mesoderm)==sign(diff_Ectoderm)] %>%
      .[,diff:=(diff_Mesoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  }

  
  if (opts$diff.type==3) {
    
    # Mesoderm-specific
    diff.acc.mes <- lapply(list("H3K27ac_distal_E7.5_Mes_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Mesoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Mesoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Ectoderm==T & sign(diff_Endoderm)==sign(diff_Ectoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Ectoderm-specific
    diff.acc.ect <- lapply(list("H3K27ac_distal_E7.5_Ect_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Endoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Endoderm"],
        fread(sprintf("%s/E7.5Ectoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Mesoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Ectoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Endoderm==T & sig_Mesoderm==T & sign(diff_Endoderm)==sign(diff_Mesoderm)] %>%
      .[,diff:=(diff_Endoderm+diff_Mesoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  
    # Endoderm-specific
    diff.acc.end <- lapply(list("H3K27ac_distal_E7.5_End_intersect12"), function(x)
      rbind(
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Mesoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Mesoderm"],
        fread(sprintf("%s/E7.5Endoderm_vs_E7.5Ectoderm_%s.txt.gz",io$diff.acc,x)) %>% .[,lineage2:="Ectoderm"]
      )) %>% rbindlist() %>% .[,lineage1:="Endoderm"] %>% 
      .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.acc.diff] %>%
      data.table::dcast(id+anno+lineage1~lineage2, value.var=c("diff","padj_fdr","sig")) %>%
      .[,sig:=sig_Mesoderm==T & sig_Ectoderm==T & sign(diff_Mesoderm)==sign(diff_Ectoderm)] %>%
      .[,diff:=(diff_Mesoderm+diff_Ectoderm)/2] %>%
      .[,c("id","anno","lineage1","diff","sig")] %>% setnames("lineage1","lineage")
  }
  
  diff.acc <- do.call("rbind", list(diff.acc.mes,diff.acc.end,diff.acc.ect))
  rm(diff.acc.mes,diff.acc.ect,diff.acc.end)
}

```

### 1.9 加载注释

```R
# Methylation
anno_df.met <- anno_df %>% split(.$anno) %>%
  map2(.,names(.), function(x,y) x[id %in% diff.met[sig==T & anno==y,id]] ) %>%
  rbindlist %>% setkey(chr,start,end) %>%
  .[,anno:=stringr::str_replace_all(anno,opts$annos)]

# Accessibility
anno_df.acc <- anno_df %>% split(.$anno) %>%
  map2(.,names(.), function(x,y) x[id %in% diff.acc[sig==T & anno==y,id]] ) %>%
  rbindlist %>% setkey(chr,start,end) %>%
  .[,anno:=stringr::str_replace_all(anno,opts$annos)]
```

### 1.10 加载Met数据（较大）

```R
# Load methylation data
met_list <- list()
for (cell in opts$met.cells) {
  tmp <- fread(sprintf("%s/%s.tsv.gz",io$met.dir,cell), showProgress=F) %>%
    .[,c("chr","pos","rate")] %>% .[,id_met:=cell] %>% 
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    
    foverlaps(.,anno_df.met, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=opts$met.tile*round(dist/opts$met.tile)] %>%
    .[,list(rate=mean(rate), n=.N),by=.(id_met,id,dist,anno)]
  met_list[[cell]] <- tmp
}
met <- rbindlist(met_list) %>%
  .[,c("id_met","id","context"):=list(as.factor(id_met),as.factor(id),as.factor("CG"))]
  
rm(met_list)
```

### 1.11 加载Acc数据（较大）

```
# Load accessibility data
acc_list <- list()
for (cell in opts$acc.cells) {
  tmp <- fread(sprintf("%s/%s.tsv.gz",io$acc.dir,cell), showProgress=F) %>%
    .[,c("chr","pos","rate")] %>% .[,id_acc:=cell] %>% 
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    
    foverlaps(.,anno_df.acc, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=opts$acc.tile*round(dist/opts$acc.tile)] %>%
    .[,list(rate=mean(rate), n=.N),by=.(id_acc,id,dist,anno)]
  acc_list[[cell]] <- tmp
}
acc <- rbindlist(acc_list) %>%
  .[,c("id_acc","id","context"):=list(as.factor(id_acc),as.factor(id),as.factor("GC"))]
  
rm(acc_list)


```

### 1.12 合并数据并保存

```
# Merge data with sample metadata
met <- met %>% merge(sample_metadata, by="id_met") %>% droplevels()
acc <- acc %>% merge(sample_metadata, by="id_acc") %>% droplevels()

data <- rbind(
  met[,c("sample","stage","stage_lineage","id","anno","dist","rate","context")],
  acc[,c("sample","stage","stage_lineage","id","anno","dist","rate","context")]
)
data[,rate:=rate*100]
data[,anno:=stringr::str_replace_all(anno,opts$annos)]
saveRDS(data, "/Users/ricard/data/gastrulation/metacc/pseudobulked_profiles/lineage_enhancers/data.rds")
```

## 2、Met与Acc增强子分析

### 2.1 加载全局Met与Acc比例

```R
met.stats <- fread(io$met.stats) %>% .[,c("id_met","mean")] %>%
  merge(sample_metadata[,.(sample,id_met)], by="id_met") %>% .[,context:="CG"]

acc.stats <- fread(io$acc.stats) %>% .[,c("id_acc","mean")] %>%
  merge(sample_metadata[,.(sample,id_acc)], by="id_acc") %>% .[,context:="GC"]

stats <- rbind(
  met.stats[,c("sample","mean","context")],
  acc.stats[,c("sample","mean","context")]
) %>% merge(sample_metadata[,c("sample","stage","stage_lineage")],by="sample") %>%
  .[,.(mean=mean(mean)),by=c("stage_lineage","context")]
  
data[,stage_lineage:=stringr::str_replace_all(stage_lineage,"_"," ")]
stats[,stage_lineage:=stringr::str_replace_all(stage_lineage,"_"," ")]
```

### 2.2 绘制Met与Acc 增强子图谱

```R
p_list <- list()

for (i in unique(data$stage_lineage)) {
  print(i)
  
  tmp <- data[stage_lineage==i]
  
  p_list[[i]] <- ggplot(tmp, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
    facet_wrap(~anno, nrow=1, scales="fixed") +
    stat_summary(geom="ribbon", fun.data="mean_se", alpha=1) +
    stat_summary(geom="line", fun.data="mean_se") +
    geom_hline(yintercept=stats[context=="CG" & stage_lineage==i,median(mean,na.rm=T)], color="#F37A71", linetype="dashed", alpha=0.75, size=0.75) +
    geom_hline(yintercept=stats[context=="GC" & stage_lineage==i,median(mean,na.rm=T)], color="#00BFC4", linetype="dashed", alpha=0.75, size=0.75) +
    labs(x="Distance from center (bp)", y="Met/Acc levels (%)") +
    coord_cartesian(ylim=c(0,100)) +
    # scale_x_continuous(breaks=c(-1,0,1)) +
    xlim(-opts$window_size, opts$window_size) +
    guides(fill=FALSE, color=FALSE, linetype=FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(0.8), colour="black"),
      axis.text.y = element_text(size=rel(1.2), colour="black")
    )
  # print(p_list[[i]])

  #pdf(file=sprintf("%s/%s.pdf",io$pdfdir,i), width=8.5, height=5)
  print(p_list[[i]])
}
```

最终得到结果，这里只展示其中两幅图

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/image-20211008234859116.png" alt="image-20211008234859116" style="zoom:50%;" />

<img src="https://raw.githubusercontent.com/Starlitnightly/bioinformatic_galaxy/master/img/image-20211008234956453.png" alt="image-20211008234956453" style="zoom:50%;" />