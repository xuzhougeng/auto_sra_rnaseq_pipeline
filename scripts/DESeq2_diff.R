#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(!requireNamespace("DESeq2",quietly = TRUE))
    stop("DESeq2 is not installed")

if(!requireNamespace("data.table",quietly = TRUE) )
    stop("data.table is not installed")

if(!requireNamespace("ashr",quietly = TRUE))
    stop("ashr is not installed")

suppressMessages(library(DESeq2))
suppressMessages(library(data.table))

### kosa_diff

### 00.data path
countspath <- args[1]
metapath <- args[2]
kosadatapath <- args[3]

# countspath <- "03_merged_counts/dataset0_GSE113179_p53.tsv"
# metapath <- "metadata/GSE113179_p53_SRA_P53-1.txt"
# kosadatapath <- "kosaRdata/dataset0_GSE113179_p53.Rdata"

### 1.exprdata
exprSet <-  data.table::fread(countspath, data.table = F)
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
keep <- grepl("ENSG",rownames(exprSet))
exprSet <- exprSet[keep,]

### 2.metadata
metadata <- data.table::fread(metapath,data.table = F)

### 3.diff analysis

dds <- DESeqDataSetFromMatrix(countData=exprSet, colData=metadata, design=~group,tidy=F)
dds <- dds[rowSums(counts(dds))>1,]
dds <- DESeq(dds)
contrast <- c("group", "treat", "control")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
### MAplot
### plotMA(dd1, ylim=c(-5,5))
### lfcShrink
dd2 <- lfcShrink(dds, contrast=contrast, res=dd1,type="ashr")
### plotMA(dd2, ylim=c(-5,5))
### diffResults
diffResults <- as.data.frame(dd2)
diffResults <- cbind(gene_id=rownames(diffResults),diffResults)

### 4.save data
kosadata <- list(exprSet=exprSet,metadata=metadata,diffResults=diffResults)
save(kosadata,file = kosadatapath)
