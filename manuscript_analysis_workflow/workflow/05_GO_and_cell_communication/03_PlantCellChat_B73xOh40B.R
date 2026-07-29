project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(ggplot2)
library(scales)
library(readxl)
library(PlantCellChat)
library(readr)
library(svglite)

message(">>> 加载 Seurat 对象...")
seob <- readRDS(file.path(project_root, "data", "scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS"))
seob <- subset(seob, orig.ident == "B73xOh40B" & Annotation %in% c("Cortex", "Epidermis"))

message(">>> 细胞数: ", ncol(seob))
message(">>> 细胞类型分布:")
print(table(seob$Annotation, seob$seurat_clusters))

seob$celltype <- paste(seob$Annotation, seob$seurat_clusters, sep = "_")
message(">>> celltype 分类:")
print(table(seob$celltype))

message(">>> 创建 PlantCellChat 对象...")
pccob <- CreatePlantCellChat(
  object   = seob,
  meta     = seob@meta.data,
  input.ident = "celltype",
  assay    = "SCT"
)

zma_v4 <- read_csv("zma_v4.csv", show_col_types = FALSE)
zma_v4 <- zma_v4 %>% distinct(Interaction_name, .keep_all = TRUE)
message(">>> 加载 LR 数据库: ", nrow(zma_v4), " 行 (去重后)")
pccob@database$interaction <- zma_v4

message(">>> 重建 genelist...")
all_ids <- unique(c(zma_v4$Ligand, zma_v4$Receptor))
all_ids <- all_ids[!is.na(all_ids)]
gene_only <- all_ids[grepl("^Zm00001d", all_ids)]
pccob@database$genelist <- data.frame(Gene_ID = gene_only, stringsAsFactors = FALSE)
message(">>> genelist 基因数: ", nrow(pccob@database$genelist),
        " (SCT中: ", sum(gene_only %in% rownames(seob@assays$SCT)), ")")

message(">>> ExtractSignalingData...")
pccob <- ExtractSignalingData(pccob)

message(">>> IdentifyOverExpressedGenes...")
pccob <- IdentifyOverExpressedGenes(pccob)

message(">>> ExtractOverExpressedInteractions...")
pccob <- ExtractOverExpressedInteractions(pccob)

source("patch_PlantCellChat.R")

message(">>> CalculateAvgExp...")
pccob <- CalculateAvgExp(pccob, methods = "average")

message(">>> CalculateCommunStrength (100 permutations)...")
pccob <- CalculateCommunStrength(pccob)

message(">>> 跳过 CalculateSignalingStrength (数据库仅含 Interacting protein 类型)")

saveRDS(pccob, "B73xOh40B_PlantCellChat.rds")
message(">>> PlantCellChat 对象已保存至 B73xOh40B_PlantCellChat.rds")

message(">>> 提取通讯结果...")

net_prob <- pccob@net$prob
net_pval <- pccob@net$pvalue
lr_signaling <- pccob@lrpairs$lr.signaling

targets  <- dimnames(net_prob)[[1]]
sources  <- dimnames(net_prob)[[2]]
lr_indices <- as.integer(dimnames(net_prob)[[3]])

message(">>> net$prob 维度: ", paste(dim(net_prob), collapse=" × "))
message(">>> 细胞类型: ", paste(targets, collapse=", "))

comm_list <- list()
for (idx in seq_along(lr_indices)) {
  lr_row <- lr_indices[idx]
  lr_info <- lr_signaling[lr_row, ]
  prob_mat <- net_prob[, , idx]
  pval_mat <- net_pval[, , idx]

  for (ti in seq_along(targets)) {
    for (si in seq_along(sources)) {
      comm_list[[length(comm_list) + 1]] <- data.frame(
        Interaction_name = lr_info$Interaction_name,
        Signal           = lr_info$Signal,
        Ligand           = lr_info$Ligand,
        Receptor         = lr_info$Receptor,
        Source           = sources[si],
        Target           = targets[ti],
        Prob             = prob_mat[ti, si],
        Pvalue           = pval_mat[ti, si],
        Interaction_type = lr_info$Interaction_type,
        stringsAsFactors = FALSE
      )
    }
  }
}
comm_df <- bind_rows(comm_list)
message(">>> 通讯结果总行数: ", nrow(comm_df))

comm_sig <- comm_df %>% filter(Pvalue < 0.05 & Prob > 0.05)
message(">>> 显著通讯 (P<0.05): ", nrow(comm_sig), " 行")

message(">>> 筛选 Cortex-Epidermis 通讯...")

comm_sig <- comm_sig %>%
  mutate(
    source_tissue = gsub("_[0-9]+$", "", Source),
    target_tissue = gsub("_[0-9]+$", "", Target)
  )

comm_ce <- comm_sig %>%
  filter(source_tissue != target_tissue) %>%
  filter(source_tissue %in% c("Cortex", "Epidermis"),
         target_tissue %in% c("Cortex", "Epidermis")) %>%
  mutate(Direction = paste0(source_tissue, "→", target_tissue))

message(">>> Cortex-Epidermis 显著通讯: ", nrow(comm_ce), " 行")
message(">>> 各方向通讯数:")
print(table(comm_ce$Direction))

message(">>> 读取 ROS 基因列表...")
ros_genes_df <- read_excel(file.path(project_root, "data", "spatial_fish_36_genes.xlsx"))
ros_genes <- ros_genes_df %>%
  filter(Mod == "ROS") %>%
  pull(geneID)

message(">>> ROS 基因 (", length(ros_genes), " 个):")
print(ros_genes)

comm_ce_ros <- comm_ce %>%
  filter(sapply(seq_len(n()), function(i) {
    any(ros_genes %in% c(Ligand[i], Receptor[i]))
  }))

message(">>> Cortex-Epidermis ROS 相关通讯: ", nrow(comm_ce_ros), " 行")

if (nrow(comm_ce_ros) > 0) {
  message(">>> ROS 通讯详情:")
  print(comm_ce_ros %>%
          select(Interaction_name, Signal, Ligand, Receptor, Direction, Prob, Pvalue) %>%
          as.data.frame())
}

all_spatial_genes <- ros_genes_df$geneID

comm_ce_spatial <- comm_ce %>%
  filter(sapply(seq_len(n()), function(i) {
    any(all_spatial_genes %in% c(Ligand[i], Receptor[i]))
  }))

message(">>> Cortex-Epidermis spatial_fish (36基因) 相关通讯: ", nrow(comm_ce_spatial), " 行")

message(">>> 保存结果...")

if (nrow(comm_ce_ros) > 0) {
  write_csv(comm_ce_ros, "B73xOh40B_ROS_communications.csv")
  message("  已保存: B73xOh40B_ROS_communications.csv")
}

if (nrow(comm_ce_spatial) > 0) {
  write_csv(comm_ce_spatial, "B73xOh40B_spatial_fish_communications.csv")
  message("  已保存: B73xOh40B_spatial_fish_communications.csv")
}

write_csv(comm_ce, "B73xOh40B_all_CE_communications.csv")
message("  已保存: B73xOh40B_all_CE_communications.csv")

message(">>> 绘图...")
dir.create("figures", showWarnings = FALSE)

nclusters <- length(unique(pccob@idents))
color_list <- rep(c('#1f7ab4', '#fd8d3c', '#66c1a7', '#e4191c',
                    '#a6cee5', '#bc80bd', '#ffed6f', 'grey60'),
                  length.out = nclusters)

setwd(project_root)
svglite("paracrine_network.svg", width = 8, height = 8)
PlottingCommunNetwork(pccob,
                      ligand.type   = "lrpairs",
                      comm.pattern  = "paracrine",
                      input.color   = color_list)
dev.off()

svglite("autocrine_chord.svg", width = 8, height = 8)
PlottingCommunNetwork(pccob,
                      ligand.type   = "lrpairs",
                      comm.pattern  = "autocrine",
                      plot.type     = "chord",
                      input.color   = color_list)
dev.off()

svglite("paracrine_chord.svg", width = 8, height = 8)
PlottingCommunNetwork(pccob,
                      ligand.type   = "lrpairs",
                      comm.pattern  = "paracrine",
                      plot.type     = "chord",
                      input.color   = color_list)
dev.off()

unique(comm_ce$Interaction_name)

library(ggplot2)
library(dplyr)

df <- comm_ce[order(-comm_ce$Prob), ][1:40, ]

extra_pairs <- c("Zm00001d030557_Zm00001d022160",
                 "Zm00001d051804_Zm00001d022160")

extra_rows <- comm_ce[comm_ce$Interaction_name %in% extra_pairs, ]
df <- unique(rbind(df, extra_rows))

df_plot <- df %>%
  mutate(Source_Target = paste(Source, "->", Target)) %>%
  mutate(neg_log10_Pvalue = -log10(Pvalue+0.000001))

p <- ggplot(df_plot, aes(x = Source_Target, y = Interaction_name)) +

  geom_point(aes(size = Prob, color = neg_log10_Pvalue)) +

  scale_color_gradient(low = "#4575b4", high = "#d73027", name = "-log10(P-value)") +

  scale_size_continuous(name = "Probability") +

  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +

  labs(
    x = "Cell-Cell Communication (Source -> Target)",
    y = "Interaction Name",
    title = "Receptor-Ligand Interactions"
  )

svglite("interaction_dotplot.svg", width = 6, height = 8)
print(p)
dev.off()
