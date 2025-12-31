if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 

BiocManager::install("DESeq2")
BiocManager::install("DESeq2")
library("DESeq2")

BiocManager::install("RColorBrewer")
BiocManager::install("pheatmap")
getwd
setwd("E:/Educational Courses/NGS course/RNAseq/10/files")
gene <- read.csv (file = "gene_count_matrix.csv", header = T, sep= ",", row.names = 1)
View(gene)
sample_key <- read.csv("sample_key.csv")
sample_key
View(sample_key)
view ("gene")
sample_key$condition <- factor ((c("1","2", "2","1", "2")))
sample_key
row.names(sample_key)<- sample_key$ sample
sample_key$sample_id <- NULL
library ("DESeq2")
dds<- DESeqDataSetFromMatrix (countData= gene, colData=sample_key, design = ~ condition)
dds <- DESeqDataSetFromMatrix(
  countData = gene,
  colData   = sample_key,
  design    = ~ condition
)
resultsNames (dds)
sample_key
dds <- DESeqDataSetFromMatrix(countData = gene, colData = sample_key, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
resgroup_gene <- results(dds, contrast= c("condition", "2", "1"))
write.csv(resgroup_gene, "resgroup_gene.csv")
resgroup_p05_gene <- subset (resgroup_gene, padj<0.05)
resgroup_p05_fold2_gene <- subset (resgroup_p05_gene, abs (log2FoldChange)> log2(2))
write.csv(resgroup_p05_fold2_gene, "resgrop_p05, fold2_gene.csv")
dim (resgroup_p05_fold2_gene)
rld <- rlogTransformation (dds)
library ("pheatmap")
if (!require("BiocManager", quietly = TRUE,))
  install.packages("BiocManager")
BiocManager:: install("pheatmap")
mat <- assay(rld) [head (order(resgroup_gene$padj), 20), ]
mat <- mat- rowMeans(mat)
df <- as.data.frame( colData(rld) [ , c("group")])
colnames(df)<- "sample"
row.names(df)<- colnames(mat)
pheatmap(mat, annotation_col=df)
library(pheatmap)
vsdata= vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup= "group")
plotMA(dds)
install.packages("pheatmap")
