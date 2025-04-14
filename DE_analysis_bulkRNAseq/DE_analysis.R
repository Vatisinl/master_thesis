library(dplyr)
library(ggplot2)
library(DOSE)
library(pheatmap)
library(tibble)
library(stringr)
library(DESeq2)
library(readr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(tidyr)
library(ggridges)
library(msigdbr)

setwd("C:/Users/Eva/Desktop/lab/DE")


## READING ##

gene_counts <- readRDS("salmon.merged.gene_counts.rds")
data <- assay(gene_counts, "counts")

meta <- data.frame(matrix(ncol = 3, nrow = 1))
colnames(meta) <- c('treatment', 'sample', 'individual')
for (i in colnames(data)) {
  f <- substring(i, 1, 2)
  id <- substring(i, 3)
  n <- substring(i, 3, 3)
  
  row_for_meta <- tibble(
    treatment = as.factor(f),
    sample = id,
    individual = as.factor(n)
  )
  meta <- bind_rows(meta, row_for_meta)
}
meta <- meta[-1, ]
rownames(meta) <- gene_counts$names

#To filter out counts <10
data <- data[rowSums(data[,]) > 10, ]







## VISUALIZATION ##
#To show - uncomment plots

dds <- DESeqDataSetFromMatrix(countData = round(data), 
                              colData = meta, 
                              design = ~ individual + treatment)

#is needed?
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

nts <- log2(assay(dds, normalized = TRUE)+1)
# ggplot(as.data.frame(nts), aes(x = ND1_S7, y = NS1_S8)) + 
# geom_point()+
# ggtitle('log2(x+1) transformation')+
# theme_dose(16)

#perform rlog-transformation
rld <- rlog(dds, blind = TRUE)

# ggplot(as.data.frame(assay(rld)), aes(x = ND1_S7, y = NS1_S8)) + 
# geom_point()+
# ggtitle('rlog transformation')+
# theme_dose(16)

plotPCA.mystyle <- function (object, ntop = 500, returnData = FALSE)
{
  font.size <- 18
  rv <- rowVars(assay(object), useNames = TRUE)
  r <- assay(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  d1 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], 
                   individual = meta$individual, 
                   treatment = meta$treatment,  
                   name = colnames(object))
  
  d1$treatment <- factor(d1$treatment)
  levels(d1$treatment) <- c("WT", "KO")
  
  ggplot(data = d1, aes_string(x = "PC1", y = "PC2")) +
    geom_point(aes_string(color = "treatment", shape = "individual"), size = 6) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_dose(font.size = font.size)+ #+ geom_label_repel(aes(label = colnames(data)), label.size = 0.1, box.padding = 0.2)
    theme(
      legend.key = element_rect(colour = NA, fill = NA), 
      legend.title= element_blank(), 
      legend.text=element_text(size=font.size-2)
    )
} 

plotPCA.mystyle(rld)
ggsave(
  'PCA.png',
  plot = last_plot(),
  dpi = 1600
)

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, annotation = meta[,c(1,3)])

meta_for_plot <- meta[,c(1,3)]
levels(meta_for_plot$treatment) <- c(levels(meta_for_plot$treatment), 'Control', 'Experiment')
meta_for_plot[meta_for_plot$treatment == 'ND',]$treatment <- 'Control'
meta_for_plot[meta_for_plot$treatment == 'NS',]$treatment <- 'Experiment'
pheatmap(rld_cor, annotation = meta_for_plot)


## DATA PRE-ANALYSIS ##
dds <- DESeqDataSetFromMatrix(countData = round(data), 
                              colData = meta, 
                              design = ~ individual+ treatment)

model.matrix(design(dds), data = colData(dds))
dds_analysis <- DESeq(dds)

#plotDispEsts(dds_analysis)

contrast <- c('treatment', 'NS', 'ND')

res_unshrunken <- results(dds_analysis, contrast=contrast, alpha = 0.05)

res <- lfcShrink(dds_analysis, 
                 contrast=contrast, 
                 res=res_unshrunken, 
                 type='normal')
#plotMA(res_unshrunken, ylim=c(-3,5))
#plotMA(res, ylim=c(-3,5))

#summary(res, alpha = 0.05)


## DATA ANALYSIS, WITH log2FC CUTOFF##
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

#convert the results table into a tibble
res_tb <- res %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

#select significant genes
significant <- res_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
nrow(significant)

map <- AnnotationDbi::select(org.Mm.eg.db, # database
                             keys = res_tb$gene,  # data to use for retrieval
                             columns = c("SYMBOL", "ENTREZID","GENENAME", 'ENSEMBL'), # information to retreive for given data
                             keytype = "ENSEMBL") # type of data given in 'keys' argument

#combine DE result tables with the obtained gene annotation
res_tb <- left_join(res_tb, map, by = c("gene" = "ENSEMBL"))

res_tb <- res_tb[!duplicated(res_tb$gene),]

#get the annotation only for the differentially expressed genes
significant <- res_tb %>% filter(gene %in% significant$gene)

p <- EnhancedVolcano(res_tb,
                     lab = res_tb$SYMBOL,
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'Exp vs Control',
                     subtitle = NULL,
                     pCutoff = 0.05,
                     FCcutoff = 0.58)

p

ggsave(
  'volcano_for_1.5LFC.png',
  plot = last_plot(),
  dpi = 1600
)
 
normalized_counts <- counts(dds_analysis, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var="gene")

norm_sig <- normalized_counts %>%
  filter(gene %in% significant$gene) %>% column_to_rownames('gene')

pheatmap(norm_sig, 
         cluster_rows = T,
         show_rownames = F,
         annotation = meta[,c(1,2)],
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         height = 10)

filtered <- filter(res_tb, ENTREZID != "NA")

foldchanges <- filtered$log2FoldChange
names(foldchanges) <- as.character(filtered$ENTREZID)
foldchanges <- sort(foldchanges, decreasing = TRUE)

sig <- msigdbr(species = "mouse")

msigdbr_t2g <- sig %>%
dplyr::distinct(gs_name, ensembl_gene) %>%
as.data.frame()

all_genes <- enricher(gene = significant$gene, TERM2GENE = msigdbr_t2g)

dotplot(all_genes)

sig_h <- msigdbr(species = "mouse", category = "H")

msigdbr_t2g <- sig_h %>%
  dplyr::distinct(gs_name, ensembl_gene) %>%
  as.data.frame()

hallmarks <- enricher(gene = significant$gene, TERM2GENE = msigdbr_t2g)

dotplot(hallmarks)


## DATA ANALYSIS, WITHOUT log2FC CUTOFF##
padj.cutoff <- 0.05
#lfc.cutoff <- 0.58

#convert the results table into a tibble
res_tb <- res %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

#select significant genes
significant <- res_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 0)
nrow(significant)

map <- AnnotationDbi::select(org.Mm.eg.db, # database
                             keys = res_tb$gene,  # data to use for retrieval
                             columns = c("SYMBOL", "ENTREZID","GENENAME", 'ENSEMBL'), # information to retreive for given data
                             keytype = "ENSEMBL") # type of data given in 'keys' argument

#combine DE result tables with the obtained gene annotation
res_tb <- left_join(res_tb, map, by = c("gene" = "ENSEMBL"))

res_tb <- res_tb[!duplicated(res_tb$gene),]

#get the annotation only for the differentially expressed genes
significant <- res_tb %>% filter(gene %in% significant$gene)

p <- EnhancedVolcano(res_tb,
                     lab = res_tb$SYMBOL, 
                     x = 'log2FoldChange',
                     y = 'padj', 
                     title = 'Exp vs Control',
                     subtitle = NULL,
                     pCutoff = 0.05, 
                     FCcutoff = 0.58)

p

normalized_counts <- counts(dds_analysis, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var="gene")

norm_sig <- normalized_counts %>%
  filter(gene %in% significant$gene) %>% column_to_rownames('gene')

# pheatmap(norm_sig, 
#          cluster_rows = T,
#          show_rownames = F,
#          annotation = meta[,c(1,2)],
#          border_color = NA,
#          fontsize = 10,
#          scale = "row",
#          fontsize_row = 10,
#          height = 10)
# 
# ego <- enrichGO(gene = significant$gene, 
#                 universe = res_tb$gene, 
#                 keyType = "ENSEMBL",
#                 OrgDb = org.Mm.eg.db,
#                 ont = "BP",
#                 pAdjustMethod = "BH", 
#                 pvalueCutoff = 0.05)
# 
# barplot(ego)

ego_down <- enrichGO(gene = significant[significant$log2FoldChange < 0, ]$gene,
                     universe = res_tb$gene,
                     keyType = "ENSEMBL",
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
barplot(ego_down)

ego_up <- enrichGO(gene = significant[significant$log2FoldChange > 0, ]$gene,
                   universe = res_tb$gene,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

barplot(ego_up)

ekegg <- enrichKEGG(gene = na.omit(significant$ENTREZID),
                    organism = 'mmu',
                    pvalueCutoff = 0.05,
                    universe = na.omit(res_tb$ENTREZID))



## DATA ANALYSIS FOR SAMPLES 4, 6, 7, WITH log2FC CUTOFF ##
data <- data[, c(3, 4, 5, 8, 9, 10)]
meta <- meta[c(3, 4, 5, 8, 9, 10), ]

dds <- DESeqDataSetFromMatrix(countData = round(data), 
                              colData = meta, 
                              design = ~ individual + treatment)

model.matrix(design(dds), data = colData(dds))
dds_analysis <- DESeq(dds)

#plotDispEsts(dds_analysis)

contrast <- c('treatment', 'NS', 'ND')

res_unshrunken <- results(dds_analysis, contrast=contrast, alpha = 0.05)

res <- lfcShrink(dds_analysis, 
                 contrast=contrast, 
                 res=res_unshrunken, 
                 type='normal')
plotMA(res_unshrunken, ylim=c(-3,5))
#plotMA(res, ylim=c(-3,5))

#summary(res, alpha = 0.05)

padj.cutoff <- 0.05
lfc.cutoff <- 0.58

#convert the results table into a tibble
res_tb <- res %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

#select significant genes
significant <- res_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
nrow(significant)

map <- AnnotationDbi::select(org.Mm.eg.db, # database
                             keys = res_tb$gene,  # data to use for retrieval
                             columns = c("SYMBOL", "ENTREZID","GENENAME", 'ENSEMBL'), # information to retreive for given data
                             keytype = "ENSEMBL") # type of data given in 'keys' argument

#combine DE result tables with the obtained gene annotation
res_tb <- left_join(res_tb, map, by = c("gene" = "ENSEMBL"))

res_tb <- res_tb[!duplicated(res_tb$gene),]

#get the annotation only for the differentially expressed genes
significant <- res_tb %>% filter(gene %in% significant$gene)

p <- EnhancedVolcano(res_tb,
                     lab = res_tb$SYMBOL,
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'Exp vs Control',
                     subtitle = NULL,
                     pCutoff = 0.05,
                     FCcutoff = 0.58)

#p
 
normalized_counts <- counts(dds_analysis, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var="gene")

norm_sig <- normalized_counts %>%
  filter(gene %in% significant$gene) %>% column_to_rownames('gene')

pheatmap(norm_sig,
         cluster_rows = T,
         show_rownames = F,
         annotation = meta[,c(1,2)],
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         height = 10)


filtered <- filter(res_tb, ENTREZID != "NA")

foldchanges <- filtered$log2FoldChange
names(foldchanges) <- as.character(filtered$ENTREZID)
foldchanges <- sort(foldchanges, decreasing = TRUE)

sig <- msigdbr(species = "mouse")     #, category = "H")

msigdbr_t2g <- sig %>%
  dplyr::distinct(gs_name, ensembl_gene) %>%
  as.data.frame()

all_genes <- enricher(gene = significant$gene, TERM2GENE = msigdbr_t2g)

dotplot(all_genes)

sig_h <- msigdbr(species = "mouse", category = "H")

msigdbr_t2g <- sig_h %>%
  dplyr::distinct(gs_name, ensembl_gene) %>%
  as.data.frame()

hallmarks <- enricher(gene = significant$gene, TERM2GENE = msigdbr_t2g)

dotplot(hallmarks)

#write.csv(significant, file = "significant_genes_samples_4_6_7.csv", row.names = FALSE)
write.csv(significant[significant$log2FoldChange > 0, ], file = "significant_genes_up_samples_467.csv", row.names = FALSE)
write.csv(significant[significant$log2FoldChange < 0, ], file = "significant_genes_down_samples_467.csv", row.names = FALSE)


## DATA ANALYSIS FOR SAMPLES 4, 6, 7, WITHOUT log2FC CUTOFF ##
padj.cutoff <- 0.05
#lfc.cutoff <- 0.58

#convert the results table into a tibble
res_tb <- res %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

#select significant genes
significant <- res_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 0)
nrow(significant)

map <- AnnotationDbi::select(org.Mm.eg.db, # database
                             keys = res_tb$gene,  # data to use for retrieval
                             columns = c("SYMBOL", "ENTREZID","GENENAME", 'ENSEMBL'), # information to retreive for given data
                             keytype = "ENSEMBL") # type of data given in 'keys' argument

#combine DE result tables with the obtained gene annotation
res_tb <- left_join(res_tb, map, by = c("gene" = "ENSEMBL"))

res_tb <- res_tb[!duplicated(res_tb$gene),]

#get the annotation only for the differentially expressed genes
significant <- res_tb %>% filter(gene %in% significant$gene)

p <- EnhancedVolcano(res_tb,
                     lab = res_tb$SYMBOL, 
                     x = 'log2FoldChange',
                     y = 'padj', 
                     title = 'Exp vs Control',
                     subtitle = NULL,
                     pCutoff = 0.05) 
                     #FCcutoff = 0.58)

# p

normalized_counts <- counts(dds_analysis, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var="gene")

norm_sig <- normalized_counts %>%
  filter(gene %in% significant$gene) %>% column_to_rownames('gene')


pheatmap(norm_sig, 
         cluster_rows = T,
         show_rownames = F,
         annotation = meta[,c(1,2)],
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         height = 10)

ego <- enrichGO(gene = significant$gene, 
                universe = res_tb$gene, 
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH", 
                pvalueCutoff = 0.05)

barplot(ego)

ego_down <- enrichGO(gene = significant[significant$log2FoldChange < 0, ]$gene, 
                     universe = res_tb$gene, 
                     keyType = "ENSEMBL",
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.05)
barplot(ego_down)

ego_down_results <- ego_down@result

# Create a data frame with GO category and associated genes
genes_by_category <- ego_down_results %>%
  select(Description, geneID) %>%
  mutate(geneID = strsplit(as.character(geneID), "/"))

for (i in 1:8) {
  g <- data.frame(genes_by_category[i,2])
  colnames(g) <- 'geneID'
  a <- left_join(g, map, by = c("geneID" = "ENSEMBL"))
  s <- paste(a$SYMBOL, collapse = ",")
  print(paste(genes_by_category[i,1], s, collapse = ":"))
}

ego_up <- enrichGO(gene = significant[significant$log2FoldChange > 0, ]$gene, 
                   universe = res_tb$gene, 
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05)

barplot(ego_up)

ego_up_results <- ego_up@result

# Create a data frame with GO category and associated genes
genes_by_category <- ego_up_results %>%
  select(Description, geneID) %>%
  mutate(geneID = strsplit(as.character(geneID), "/"))

for (i in 1:8) {
  g <- data.frame(genes_by_category[i,2])
  colnames(g) <- 'geneID'
  a <- left_join(g, map, by = c("geneID" = "ENSEMBL"))
  s <- paste(a$SYMBOL, collapse = ",")
  print(paste(genes_by_category[i,1], s, collapse = ":"))
}

ekegg <- enrichKEGG(gene = na.omit(significant$ENTREZID),
                    organism = 'mmu',
                    pvalueCutoff = 0.05,
                    universe = na.omit(res_tb$ENTREZID))
#write.csv(significant, file = "significant_genes_samples_without_log2FC_4_6_7.csv", row.names = FALSE)
write.csv(significant[significant$log2FoldChange > 0, ], file = "significant_genes_up_without_log2FC_samples_467.csv", row.names = FALSE)
write.csv(significant[significant$log2FoldChange < 0, ], file = "significant_genes_down_without_log2FC_samples_467.csv", row.names = FALSE)


map <- AnnotationDbi::select(org.Mm.eg.db, # database
                             keys = res_tb$gene,  # data to use for retrieval
                             columns = c("SYMBOL", "ENTREZID","GENENAME", 'ENSEMBL', 'GO'), # information to retreive for given data
                             keytype = "ENSEMBL")
#combine DE result tables with the obtained gene annotation
res_tb <- left_join(res_tb, map, by = c("gene" = "ENSEMBL"))

res_tb <- res_tb[!duplicated(res_tb$gene),]

#get the annotation only for the differentially expressed genes
significant <- res_tb %>% filter(gene %in% significant$gene)

write.csv(significant, file='significant_all_info.csv', row.names = FALSE)

write.csv(res_tb, file='universe_for_enrichment.csv', row.names = FALSE)



