library(Rsubread)
library(DESeq2)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(limma)
library(edgeR)


################################################################################
#####          LCM RNA-seq data analysis                                   #####
################################################################################

register(MulticoreParam(20))

# elegant color map 
tissue.col <- c("#C9A1CA", '#B46DA9', "#EDC66A", "#018B38", '#E87B1E',"#F9B3AD",
                "#D9DEE7","#9FDAF7", "#5CB0C3", "#8498AB", "#C74546" )

# read all bam files
all_files_lcm <- list.files(pattern = ".bam", path = "/home/go67saw/LCM_RNA/clean/bam_v4/")
all_files_lcm_path <- paste0("/home/go67saw/LCM_RNA/clean/bam_v4/", all_files_lcm)

# compute raw read counts
fc <- featureCounts(all_files_lcm_path, annot.ext = "/home/go67saw/LCM_RNA/Zea_mays_reference/v4/Zea_mays.B73_RefGen_v4.50.chr.gtf.gz", 
                    isGTFAnnotationFile = TRUE, isPairedEnd = TRUE, 
                    countMultiMappingReads = FALSE, useMetaFeatures = TRUE,
                    countChimericFragments = FALSE, nthreads = 20)

write_xlsx(data.frame(fc$counts), "/home/go67saw/LCM_RNA/intermediate_data/LCM_raw_counts_table_ref_v4.xlsx")
write.csv(fc$stat, file = "/home/go67saw/LCM_RNA/intermediate_data/featurecounts_stat.csv", row.names = F)

# check sequencing assignment
assigned.pct = data.frame(pct = t(fc$stat[1, -1]/colSums(fc$stat[, -1])))
assigned.pct = cbind.data.frame(rownames(assigned.pct), assigned.pct$X1)
colnames(assigned.pct) = c("sampleID", "assignedPercent")
write.csv(assigned.pct, "./results/rna-seq/featurecounts_stat/assigned_percent.csv", row.names = F)

# save raw gene counts table
raw.counts = fc$counts
saveRDS(raw.counts, '~/LCM_RNA/intermediate_data/LCM_raw_counts_table_ref_v4.rds')
# sample metadata creation
sample.md = data.frame(sampleID = colnames(raw.counts))
sample.md$Compartment = sub('\\..*', '', sample.md$sampleID)
sample.md$Tissue = c('cells_betw_Meta_Xylems', 'proto_xylem', ' proto_pholem',
                     'outer_Cortex', 'inner_Cortex', 'Endodermis', 'Pholem_pole_pericycle',
                     'Xylem_pole_pericycle', 'Epidermis', 'Meta_xylem', 'Pith')

write.csv(sample.md, "~/LCM_RNA/intermediate_data/sample_metadata.csv", row.names = F)

## Filtering to remove lowly expressed genes 
thresh <- raw.counts > 5
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
# Subset the rows of count data to keep the more highly expressed genes
counts.keep <- raw.counts[keep, ]
dim(counts.keep)   # 23333  x 11

write.csv(counts.keep, file = "~/LCM_RNA/intermediate_data/filtered_counts_table_ref_v4.csv", row.names = TRUE)

hist(rowSums(counts.keep), breaks = 100, freq = T)

table(sample.md$sampleID == colnames(counts.keep))
# since there is no replicates for each tissue, Deseq2 methods cannot be used!!
dge <- DGEList(counts = counts.keep, group = sample.md$Compartment)
dge <- calcNormFactors(dge, method = "TMM")
# get normalized counts
norm_counts <- edgeR::cpm(dge, log = T)
saveRDS(norm_counts, '~/LCM_RNA/intermediate_data/LCM_filtered_counts_table_ref_v4_normalized.rds')

# MDS plot
mds.pl = plotMDS(norm_counts, top=2000,  dim.plot = c(1,2), plot = FALSE)
points(mds.pl$x, mds.pl$y, col = "darkgreen", pch = 16, cex = 1.5)

pdf('~/LCM_RNA/results/LCM_RNA_MDS_plot_top2000_genes.pdf', width = 5, height = 5)
# Plot manually with custom limits
plot(mds.pl$x, mds.pl$y,
     xlim = c(-5, 4), ylim = c(-3, 4),         # << set limits here
     col = 'darkgreen',
     pch = 16,
     xlab = "MDS1 28%", ylab = "MDS2 19%",
     main = "")

# add labels
text(mds.pl$x, mds.pl$y, labels = sample.md$Tissue, pos = 3, cex = 0.8)
dev.off()

# plot gene expression
gene.exp = t(norm_counts)
gene.exp = data.frame(x=sample.md$Tissue, gene.exp)

ggplot(gene.exp, aes(x = x, y = Zm00001d027500)) +
  geom_point(size = 3)  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#-----------------------------------------------#
#      Differential expression analysis
#-----------------------------------------------#

# to set same contrast in model fitting, add LCM to let every tissue in base line
sample.md$Tissue = paste0('LCM_', sample.md$Tissue)

all.DE.genes = c()

for (t in sample.md$Tissue) {
  print(t)
  sample.md$comp2 = 'other'
  sample.md[sample.md$Tissue == t, ]$comp2 = t
  
  # Fit a quasi-likelihood model (QL) without replicates
  design <- model.matrix(~ comp2, data = sample.md)
  dge1 <- estimateDisp(dge, design, robust = TRUE)
  
  # Apply GLM fitting and testing
  fit <- glmQLFit(dge1, design)
  de_results <- glmQLFTest(fit, contrast = c(0, 1))  # 0 is control, 1 is other
  # so when genes express more in Control, logFC < 0
  results <- de_results$table
  results$FDR <- p.adjust(results$PValue, method = "BH")
  # Filter significant genes (FDR < 0.01 is more confident !!)
  DE_genes <- subset(results, PValue < 0.01 & logFC < -1)
  
  all.DE.genes <- rbind.data.frame(all.DE.genes, cbind(DE_genes, Tissue=t))
}

write_xlsx(cbind.data.frame(geneID = rownames(all.DE.genes), all.DE.genes), './LCM_glmQLfit_all_DEGs.xlsx')





sessionInfo()         # show the version of packages










