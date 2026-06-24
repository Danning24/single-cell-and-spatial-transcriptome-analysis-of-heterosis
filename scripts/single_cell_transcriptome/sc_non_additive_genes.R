library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(AnnotationHub)
library(clusterProfiler)
library(biomaRt)



################################################################################
#####          Non-additive genes in each cluster                          #####
################################################################################

setwd('~/single_cell_RNA_seq/')

scRNA_merge_filtered_int <- readRDS('./intermediate_data/new_all_samples/scRNA_all_annotated_filtered_CCA.RDS')

non.additive.test = function(hb.lev, vars.names, het.input){
  heterosis=c()
  for (h in hb.lev){ 
    for (v in vars.names) {
      mat = strsplit(h, 'x')[[1]][1] # identify maternal
      pat = strsplit(h, 'x')[[1]][2] # identify paternal
      
      # Subset sample data to only include the 3 genotypes from this cross
      mat.data = het.input[het.input$Genotype %in% c('B73_Rep1', 'B73_Rep2', 'B73_Rep3'), v] 
      pat.data = het.input[het.input$Genotype == pat, v]
      hyb.data = het.input[het.input$Genotype == h, v]
      
      # Calculate mean for hybrid, maternal, and paternal
      hyb.Mean = mean(hyb.data, na.rm=TRUE) 
      mat.Mean = mean(mat.data, na.rm=TRUE) 
      pat.Mean = mean(pat.data, na.rm=TRUE)
      # Calculate mid-parent value
      midparent = sum(mat.Mean, pat.Mean)/2 
      
      # T-test for non-additive effects test
      tryCatch(
        {tt.MPH = t.test(hyb.data, mu=midparent, alternative='two.sided')
        if(tt.MPH$p.value < 0.05){

          tt.results = data.frame('Maternal'=mat,'Paternal'=pat,'Hybrid'=h,
                                  'MatMean'=mat.Mean,'PatMean'=pat.Mean,'HybMean'=hyb.Mean,
                                  'MidParentValue'=midparent, 'Trait'=v,
                                  'MPH.t.stat'=tt.MPH$statistic,'MPH.t.df'=tt.MPH$parameter, 'MPH.p'=tt.MPH$p.value,
                                  'MPH.pct'=(hyb.Mean-midparent)/midparent)
          heterosis = rbind(heterosis, tt.results)
        }
        },
        error = function(e) NA)
    }
  }
  return(heterosis)
}

hybrids = unique(scRNA_merge_filtered_int$orig.ident)[(grepl('x', unique(scRNA_merge_filtered_int$orig.ident)))]
cluster.expr.genes = c()
het.res = c()

for (c in as.vector(unique(scRNA_merge_filtered_int$annotation))) {
  print(paste0('cluster...', c))
  sc.cluster <- subset(scRNA_merge_filtered_int, annotation==c)
  # get gene matrix
  gene_table <- GetAssayData(sc.cluster, layer = "data", assay = "SCT")
  gene.het.input = data.frame(t(gene_table), Genotype = sc.cluster@meta.data[, 'orig.ident'])
  # choose genes that expressed at least in 10% cells of  in this cluster
  vars.names <- colnames(gene.het.input)[colSums(gene.het.input >0) > 0.1*nrow(gene.het.input)]
  print(paste0('number of expressed genes is...', length(vars.names)))
  cluster.expr.genes <- rbind(cluster.expr.genes, cbind(cluster = c, geneID = vars.names))
  
  res = non.additive.test(hybrids, vars.names, gene.het.input) 
  het.res = rbind(het.res, cbind(res, cluster = c, num.gene = length(vars.names)))
}

# Adjust p-values to correct for multiple tests 
heterosis.results = group_by(het.res, Hybrid, cluster) %>% 
  mutate(MPH.padj = p.adjust(MPH.p, method='bonferroni')) %>% 
  ungroup %>% as.data.frame %>%
  mutate(HetType=case_when(
    MPH.padj < 0.05 ~ 'Non-additive',
    MPH.padj > 0.05 ~ 'Additive'))

write_xlsx(heterosis.results, 
           './results/heterotic_genes/scRNA_heterotic_genes_results_10percent_cells.xlsx')

# check which cluster has the most non-additive genes
hyb.freq = heterosis.results %>% 
  group_by(., cluster, HetType) %>% 
  summarise(., cluster.fq = n())

#####  GO enrichment analysis  #####

gene.GO <- read.csv('~/single_cell_RNA_seq/results/Zm-B73-REFERNECE-GRAMENE-4.0-GMs-GOTerms.csv')
zm.anno <- read.delim('~/spatial_RNA/intermediate_data/maize.B73.AGPv4.aggregate.gaf', header = FALSE, comment.char = "!", stringsAsFactors = FALSE)

hub <- AnnotationHub()
query(hub, "zea")
maize <- hub[['AH117408']]
columns(maize)
length(keys(maize))

selected_data <- select(maize, keys = head(keys(maize)), 
                        columns = c("ACCNUM", "ALIAS", "SYMBOL", "ENTREZID", "GID", "UNIGENE", "PMID","REFSEQ"))
head(selected_data)

require(GO.db)
all.gos <- as.list(GOTERM)

go.df <- data.frame()

for (go in names(all.gos)) {
  if(!is.na(go)){
    temp <- data.frame(
      GOID = all.gos[[go]]@GOID,
      Term = all.gos[[go]]@Term,
      Ontology = all.gos[[go]]@Ontology,
      Definition = all.gos[[go]]@Definition,
      stringsAsFactors = FALSE
    )
    go.df <- rbind.data.frame(go.df, temp)
  }
}

cluster.go <- enricher(heterosis.results[heterosis.results$cluster == 0 &
                                           heterosis.results$HetType == 'Non-additive', ]$Trait, 
                       TERM2GENE = zm.anno[, c(5,2)])                

enrich.go <- cluster.go@result
enrich.go <- enrich.go[enrich.go$p.adjust <0.05, ]
enrich.go.terms <- cbind(enrich.go, go.df[match(enrich.go$ID, go.df$GOID), ])

# save results
write_xlsx(enrich.go.terms, './results/heterotic_genes/cluster0_heterotic_genes_enriched_GO.xlsx')

# barplot for BP GO terms
g1 <- ggplot(enrich.go.terms[enrich.go.terms$Ontology == 'BP' & !is.na(enrich.go.terms$Ontology), ], 
             aes(x = Term, y = FoldEnrichment, fill=p.adjust)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  coord_flip()

ggsave(g1, filename = './results/heterotic_genes/cluster0_heterotic_genes_enriched_GO.pdf', 
       width = 25, height = 45, units = 'cm')

# find unique non-additive genes for each cluster
uniq.genes.all = NULL
for (c in unique(heterosis.results$cluster)) {
    curr.genes <- heterosis.results[heterosis.results$cluster == c, 'Trait']
    other.genes <- heterosis.results[heterosis.results$cluster != c, 'Trait']
    if(length(curr.genes) >0){
      uniq.genes <- setdiff(curr.genes, other.genes)
      if(length(uniq.genes) >0) {
        uniq.genes.all <- rbind(uniq.genes.all, cbind(cluster=c, uniq.genes))
        
      }
    }
}

uniq.genes.all <- data.frame(uniq.genes.all)

uniq.gene <- uniq.genes.all %>% 
  group_by(., cluster) %>% 
  summarize(., num.uniq.het = n()) %>% data.frame() 

g2 <- ggplot(uniq.gene, aes(x = cluster, y = num.uniq.het)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=num.uniq.het), vjust=-1, color='black', size=3)

ggsave(g2, filename = './heterotic_genes/number_unique_het_genes_in_each_cluster.pdf', 
       width = 40, height = 15, units = 'cm')






