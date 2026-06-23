library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)


################################################################################
#####      Identification of non-additive genes                            #####
################################################################################

sp_merge_filtered_int <- readRDS('~/spatial_RNA/intermediate_data/all_sp_merge_filtered_integrated_remove_noise_bin50.RDS')

# t-test for identifying non-additive genes
non.additive.test = function(vars.names, het.input ){
  
  heterosis = c()
  mat = 'B73'
  pat = 'A554'
  hyb = 'B73xA554'
  
  for (v in vars.names) {

      # Subset sample data to only include the 3 genotypes 
      mat.data = het.input[het.input$Genotype == mat, v] 
      pat.data = het.input[het.input$Genotype == pat, v]
      hyb.data = het.input[het.input$Genotype == hyb, v]
      
      # Calculate mean for hybrid, maternal, and paternal
      hyb.Mean = mean(hyb.data, na.rm=TRUE) 
      mat.Mean = mean(mat.data, na.rm=TRUE) 
      pat.Mean = mean(pat.data, na.rm=TRUE)
      
      # Calculate mid-parent value
      midparent = sum(mat.Mean, pat.Mean)/2 
      
      # T-test for non-additive effects test
      tryCatch(
        {tt.MPH = t.test(hyb.data, mu = midparent, alternative = 'two.sided')
        if(tt.MPH$p.value < 0.05){
          # Save results
          tt.results = data.frame('Maternal' = mat, 'Paternal' = pat, 'Hybrid' = hyb,
                                  'MatMean'= mat.Mean, 'PatMean' = pat.Mean, 'HybMean' = hyb.Mean,
                                  'MidParentValue' = midparent, 'Trait'=v,
                                  'MPH.t.stat' = tt.MPH$statistic, 
                                  'MPH.t.df' = tt.MPH$parameter, 
                                  'MPH.p' = tt.MPH$p.value,
                                  'MPH.pct' = (hyb.Mean-midparent)/midparent)
          heterosis = rbind(heterosis, tt.results)
        }
        },
        error = function(e) NA)
    }
    
  return(heterosis)
}

het.res = c()
# perform non-additive gene test for each cluster
for (c in as.vector(unique(sp_merge_filtered_int$cluster_annotation))) {
  print(paste0('cluster...', c))
  sp.cluster <- subset(sp_merge_filtered_int, cluster_annotation == c)
  # subset cluster gene matrix
  gene_table <- GetAssayData(sp.cluster, layer = "data", assay = "SCT")
  gene.het.input = data.frame(t(gene_table), Genotype = sp.cluster@meta.data[, 'Genotype'])
  # choose genes that expressed at least in 10% cells of  in this cluster
  vars.names <- colnames(gene.het.input)[colSums(gene.het.input >0) > 0.1*nrow(gene.het.input)]
  print(paste0('number of expressed genes is...', length(vars.names)))
  # non-additive test
  res = non.additive.test(vars.names, gene.het.input) 
  het.res = rbind(het.res, cbind(res, cluster = c, num.gene = length(vars.names)))
  
}

# Adjust p-values to correct for multiple tests 
sp.heterosis.results <- het.res[het.res$HybMean > 0 & is.finite(het.res$MPH.pct), ]
sp.heterosis.results <- group_by(sp.heterosis.results, cluster) %>% 
  mutate(MPH.padj = p.adjust(MPH.p, method='bonferroni')) %>% 
  ungroup %>% as.data.frame 

sp.heterosis.results$pattern <- ifelse(sp.heterosis.results$MPH.padj <0.05, 'Non-additive', 'Additive')

write.csv(sp.heterosis.results, 
           '~/spatial_RNA/results/integrated_3_samples/heterotic_genes/new/sp_heterotic_genes_results_10percent_cells_filtered.csv', row.names = F)


# check which cluster has the most heterotic genes
gene.freq = sp.heterosis.results %>% 
  group_by(., cluster, pattern) %>% 
  summarise(., gene.fq = n())

g1 <- ggplot(gene.freq, aes(x = cluster, y = gene.fq, fill=pattern)) +
  geom_bar(stat="identity", position = 'fill') +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percentage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(g1, filename = './heterotic_genes/new/percent_het_genes_for_each_cluster.pdf', 
       width = 16, height = 8, units = 'cm')

## proportion test ##
sp_non_add_genes_fq <- sp.heterosis.results[sp.heterosis.results$pattern == 'Non-additive', ] %>%
  group_by(cluster) %>%
  summarise(n_genes = n())

sp_genes_fq <- sp.heterosis.results %>%
  group_by(cluster) %>%
  summarise(n_genes = n())

sig.prop <- c()
for (c in sp_non_add_genes_fq$cluster) {
    #One-vs-others test
    test.res <- prop.test(
      x = c(sp_non_add_genes_fq[sp_non_add_genes_fq$cluster == c, ]$n_genes, mean(sp_non_add_genes_fq[sp_non_add_genes_fq$cluster != c, ]$n_genes)),
      n = c(sp_genes_fq[sp_genes_fq$cluster == c, ]$n_genes, mean(sp_genes_fq[sp_genes_fq$cluster != c, ]$n_genes)),
      alternative = "greater")
    print(test.res)
    if(test.res$p.value <0.05) {
      sig.prop <- rbind(sig.prop, cbind(cluster = c, 
                                        curr.cluster.prop = test.res$estimate[1],
                                        other.cluster.prop = test.res$estimate[2],
                                        pval = test.res$p.value))
    }
}

sig.prop.df <- data.frame(sig.prop) %>% mutate(., padj = p.adjust(pval, method = 'fdr'))
write.csv(sig.prop.df, './heterotic_genes/new/sig_cluster_for_heterotic_genes_prop_test_in_each_cluster.csv', row.names = F)






