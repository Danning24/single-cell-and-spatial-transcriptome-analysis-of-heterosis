# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

setwd(project_root)
library(topGO)
library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(AnnotationHub)
library(GO.db)
library(dplyr)

# Dependency only; install outside the analysis script: BiocManager::install("topGO")
# Configure proxies outside the analysis script if required.
# Configure proxies outside the analysis script if required.
# Optional proxy configuration removed from the public workflow.
# Optional proxy configuration removed from the public workflow.
# Optional proxy configuration removed from the public workflow.
hub = AnnotationHub()
# query(hub, "zea")
maize = hub[["AH114308"]]

nonaddgenes <- read_csv(file.path(project_root, "data", "sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv"))
nonadd_cortex <- filter(nonaddgenes,Hybrid == "B73xOh40B" & pattern == "Non-additive" &
                          cluster_anno %in% c("Cortex early",
                                              "Cortex I mature",
                                              "Cortex II mature",
                                              "Cortex outer"))
nonadd_cortex <-  nonadd_cortex[,c(4,5,8)]
keytypes(maize)
write.csv(nonadd_cortex,"nonadd_cortex_b73xoh40b.csv")

#在maizegdb上进行id转化
file="gene_translate.php.txt"
gene <- read.delim(file,header = T)
deTable = merge(nonadd_cortex, gene,by.x="Trait", by.y="input", all=F)

res <- deTable %>%
  group_by(GenBank.Gene) %>%
  summarise(
    status = case_when(
      all(MPH.pct > 0, na.rm = TRUE) ~ "up_cortex",
      all(MPH.pct < 0, na.rm = TRUE) ~ "down_cortex",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(status))

up_cortex <- res %>%
  filter(status == "up_cortex") %>%
  pull(GenBank.Gene)

down_cortex <- res %>%
  filter(status == "down_cortex") %>%
  pull(GenBank.Gene)

gene_lists <- list(
  up_cortex = up_cortex,
  down_cortex = down_cortex
)


for (i in names(gene_lists)) {
  current_genes <- gene_lists[[i]]
  cat("正在处理:", i, "- 基因数量:", length(current_genes), "\n")
  gsea_res<-enrichGO(gene = current_genes,OrgDb = maize,
                     keyType = "ENTREZID", ont = "ALL",
           qvalueCutoff = 0.05,readable = TRUE)
  file_name = paste0("GO_res_",i,".txt")
  write.table(gsea_res, file_name, quote=F, sep="\t")
}

# 如果 GeneRatio 是 "12/345" 这种格式，转成数值
parse_geneRatio <- function(x){
  if (is.numeric(x)) return(x)
  sapply(x, function(i){
    if (grepl("/", i)) {
      a <- as.numeric(strsplit(i, "/")[[1]][1])
      b <- as.numeric(strsplit(i, "/")[[1]][2])
      return(a / b)
    } else {
      return(as.numeric(i))
    }
  })
}

parse_geneRatio <- function(x){
  if (is.numeric(x)) return(x)
  sapply(x, function(i){
    if (grepl("/", i)) {
      sp <- strsplit(i, "/")[[1]]
      as.numeric(sp[1]) / as.numeric(sp[2])
    } else {
      as.numeric(i)
    }
  })
}

go_down <- read.delim("GO_res_down_cortex.txt",header = T)
go_up <- read.delim("GO_res_up_cortex.txt",header = T)


go_up$GeneRatio_num   <- parse_geneRatio(go_up$GeneRatio)
go_down$GeneRatio_num <- parse_geneRatio(go_down$GeneRatio)

go_up$logP   <- -log10(go_up$p.adjust)
go_down$logP <- -log10(go_down$p.adjust)

top_up <- go_up %>%
  arrange(p.adjust) %>%
  slice(1:10) %>%
  mutate(
    group = "up",
    row_id = 20:11
  )

top_down <- go_down %>%
  arrange(p.adjust) %>%
  slice(1:10) %>%
  mutate(
    group = "down",
    row_id = 10:1
  )

plot_df <- bind_rows(top_up, top_down)

# 不换行，直接用原始通路名
plot_df$label <- plot_df$Description

# up / down 三角位置
plot_df$x <- ifelse(plot_df$group == "up", 1.2, 1.8)

# 标签统一放左侧框外
label_df <- plot_df %>%
  mutate(label_x = 0.15)

p <- ggplot(plot_df, aes(x = x, y = row_id)) +
  geom_text(
    data = label_df,
    aes(x = label_x, y = row_id, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    size = 4
  ) +
  geom_point(
    data = subset(plot_df, group == "up"),
    aes(size = GeneRatio_num, fill = logP),
    shape = 24,
    color = "black",
    stroke = 0.4
  ) +
  geom_point(
    data = subset(plot_df, group == "down"),
    aes(size = GeneRatio_num, fill = logP),
    shape = 25,
    color = "black",
    stroke = 0.4
  ) +
  scale_x_continuous(
    breaks = c(1.2, 1.8),
    labels = c("up", "down"),
    limits = c(0.9, 2.4)
  ) +
  scale_y_continuous(breaks = NULL) +
  scale_size_continuous(name = "GeneRatio") +
  scale_fill_gradient(
    low = "#c7e9c0",
    high = "#006d2c",
    name = expression(-log[10](p.adjust))
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    text = element_text(size = 13),
    plot.margin = margin(10, 30, 10, 220)
  )

p

ggsave("GO_top10_up_down_separate_green_outside_labels.pdf",
       p, width = 12, height = 10)
ggsave("GO_top10_up_down_separate_green_outside_labels.png",
       p, width = 12, height = 10, dpi = 300)




##表皮
nonadd_Epidermis <- filter(nonaddgenes,Hybrid == "B73xA554" & pattern == "Non-additive" &
                          cluster_anno %in% c("Epidermis"))
nonadd_Epidermis <-  nonadd_Epidermis[,c(4,5)]
keytypes(maize)
write.csv(nonadd_Epidermis,"nonadd_Epidermis_b73xa554.csv")

#在maizegdb上进行id转化
file="gene_translate_epidermics_php.txt"
gene <- read.delim(file,header = T)
deTable = merge(nonadd_Epidermis, gene,by.x="Trait", by.y="input", all=F)

res <- deTable %>%
  group_by(GenBank.Gene) %>%
  summarise(
    status = case_when(
      all(MPH.pct > 0, na.rm = TRUE) ~ "up_epi",
      all(MPH.pct < 0, na.rm = TRUE) ~ "down_epi",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(status))

up_epi <- res %>%
  filter(status == "up_epi") %>%
  pull(GenBank.Gene)

down_epi <- res %>%
  filter(status == "down_epi") %>%
  pull(GenBank.Gene)

gene_lists <- list(
  up_epi = up_epi,
  down_epi = down_epi
)


for (i in names(gene_lists)) {
  current_genes <- gene_lists[[i]]
  cat("正在处理:", i, "- 基因数量:", length(current_genes), "\n")
  gsea_res<-enrichGO(gene = current_genes,OrgDb = maize,
                     keyType = "ENTREZID", ont = "ALL",
                     qvalueCutoff = 0.05,readable = TRUE)
  file_name = paste0("GO_res_",i,".txt")
  write.table(gsea_res, file_name, quote=F, sep="\t")
}


go_down <- read.delim("GO_res_down_epi.txt",header = T)
go_up <- read.delim("GO_res_up_epi.txt",header = T)

go_up$GeneRatio_num   <- parse_geneRatio(go_up$GeneRatio)
go_down$GeneRatio_num <- parse_geneRatio(go_down$GeneRatio)

go_up$logP   <- -log10(go_up$p.adjust)
go_down$logP <- -log10(go_down$p.adjust)

top_up <- go_up %>%
  arrange(p.adjust) %>%
  slice(1:10) %>%
  mutate(
    group = "up",
    row_id = 20:11
  )

top_down <- go_down %>%
  arrange(p.adjust) %>%
  slice(1:10) %>%
  mutate(
    group = "down",
    row_id = 10:1
  )

plot_df <- bind_rows(top_up, top_down)

# 不换行，直接用原始通路名
plot_df$label <- plot_df$Description

# up / down 三角位置
plot_df$x <- ifelse(plot_df$group == "up", 1.2, 1.8)

# 标签统一放左侧框外
label_df <- plot_df %>%
  mutate(label_x = 0.15)

p <- ggplot(plot_df, aes(x = x, y = row_id)) +
  geom_text(
    data = label_df,
    aes(x = label_x, y = row_id, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    size = 4
  ) +
  geom_point(
    data = subset(plot_df, group == "up"),
    aes(size = GeneRatio_num, fill = logP),
    shape = 24,
    color = "black",
    stroke = 0.4
  ) +
  geom_point(
    data = subset(plot_df, group == "down"),
    aes(size = GeneRatio_num, fill = logP),
    shape = 25,
    color = "black",
    stroke = 0.4
  ) +
  scale_x_continuous(
    breaks = c(1.2, 1.8),
    labels = c("up", "down"),
    limits = c(0.9, 2.4)
  ) +
  scale_y_continuous(breaks = NULL) +
  scale_size_continuous(name = "GeneRatio") +
  scale_fill_gradient(
    low = "#c7e9c0",
    high = "#006d2c",
    name = expression(-log[10](p.adjust))
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    text = element_text(size = 13),
    plot.margin = margin(10, 30, 10, 220)
  )

p

ggsave("GOepi_top10_up_down_separate_green_outside_labels.pdf",
       p, width = 12, height = 10)
ggsave("GOepi_top10_up_down_separate_green_outside_labels.png",
       p, width = 12, height = 10, dpi = 300)


plot_df$Description
