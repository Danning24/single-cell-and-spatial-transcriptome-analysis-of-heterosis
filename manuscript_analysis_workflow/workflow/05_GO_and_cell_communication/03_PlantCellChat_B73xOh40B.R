# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

# =========================================================================
# PlantCellChat: B73xOh43 玉米根 Cortex-Epidermis 细胞通讯 × ROS 基因筛选
# =========================================================================

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(ggplot2)
library(scales)
library(readxl)
library(PlantCellChat)
library(readr)
library(svglite)

# =========================================================================
# 1. 加载数据 & 筛选 B73xOh43 的 Cortex 和 Epidermis
# =========================================================================
message(">>> 加载 Seurat 对象...")
seob <- readRDS(file.path(project_root, "data", "scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS"))
seob <- subset(seob, orig.ident == "B73xOh40B" & Annotation %in% c("Cortex", "Epidermis"))

message(">>> 细胞数: ", ncol(seob))
message(">>> 细胞类型分布:")
print(table(seob$Annotation, seob$seurat_clusters))

# 创建精细的 celltype 标识: Annotation_seurat_cluster
seob$celltype <- paste(seob$Annotation, seob$seurat_clusters, sep = "_")
message(">>> celltype 分类:")
print(table(seob$celltype))

# =========================================================================
# 2. 创建 PlantCellChat 对象 & 加载 v4 基因数据库
# =========================================================================
message(">>> 创建 PlantCellChat 对象...")
pccob <- CreatePlantCellChat(
  object   = seob,
  meta     = seob@meta.data,
  input.ident = "celltype",
  assay    = "SCT"
)

# 加载已转换的 v4 基因ID 配体-受体数据库
zma_v4 <- read_csv("zma_v4.csv", show_col_types = FALSE)
# 去除重复的 Interaction_name（可能出现于 v5→v4 转换中的 1:多映射）
zma_v4 <- zma_v4 %>% distinct(Interaction_name, .keep_all = TRUE)
message(">>> 加载 LR 数据库: ", nrow(zma_v4), " 行 (去重后)")
pccob@database$interaction <- zma_v4



# 重建 genelist（替换数据库后必须手动重建，因为 PlantCellChat 不会自动更新）
# 只保留 Zm00001d 格式的基因ID，排除小分子信号（如 ABA, BR, Ca2+）
message(">>> 重建 genelist...")
all_ids <- unique(c(zma_v4$Ligand, zma_v4$Receptor))
all_ids <- all_ids[!is.na(all_ids)]
gene_only <- all_ids[grepl("^Zm00001d", all_ids)]
pccob@database$genelist <- data.frame(Gene_ID = gene_only, stringsAsFactors = FALSE)
message(">>> genelist 基因数: ", nrow(pccob@database$genelist),
        " (SCT中: ", sum(gene_only %in% rownames(seob@assays$SCT)), ")")

# =========================================================================
# 3. 运行 PlantCellChat 管线
# =========================================================================
message(">>> ExtractSignalingData...")
pccob <- ExtractSignalingData(pccob)

message(">>> IdentifyOverExpressedGenes...")
pccob <- IdentifyOverExpressedGenes(pccob)

message(">>> ExtractOverExpressedInteractions...")
pccob <- ExtractOverExpressedInteractions(pccob)

# ★ 加载 PlantCellChat bug 热修复补丁
# （修复 ExtractExpMatrix 的 tibble [,] → [[ ]] 问题）
source("patch_PlantCellChat.R")

message(">>> CalculateAvgExp...")
pccob <- CalculateAvgExp(pccob, methods = "average")

message(">>> CalculateCommunStrength (100 permutations)...")
pccob <- CalculateCommunStrength(pccob)

# CalculateSignalingStrength 仅适用于 "Signaling molecule" 类型，
# 方案1过滤后只有 Interacting protein 类型，跳过此步骤。
# 通讯结果已在 CalculateCommunStrength 的 net slot 中。
message(">>> 跳过 CalculateSignalingStrength (数据库仅含 Interacting protein 类型)")


# 保存 PlantCellChat 对象（方便后续直接加载分析）
saveRDS(pccob, "B73xOh40B_PlantCellChat.rds")
message(">>> PlantCellChat 对象已保存至 B73xOh40B_PlantCellChat.rds")

# =========================================================================
# 4. 提取所有通讯结果 (从 net slot 手动构建)
# =========================================================================
message(">>> 提取通讯结果...")

# net$prob 维度: [Target, Source, Interaction_index]
# dimnames: target celltypes × source celltypes × interaction row numbers
net_prob <- pccob@net$prob
net_pval <- pccob@net$pvalue
lr_signaling <- pccob@lrpairs$lr.signaling

targets  <- dimnames(net_prob)[[1]]
sources  <- dimnames(net_prob)[[2]]
lr_indices <- as.integer(dimnames(net_prob)[[3]])

message(">>> net$prob 维度: ", paste(dim(net_prob), collapse=" × "))
message(">>> 细胞类型: ", paste(targets, collapse=", "))

# 构建完整的通讯结果表 (Target×Source×Interaction → long format)
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

# 筛选 Pvalue < 0.05 的显著通讯
comm_sig <- comm_df %>% filter(Pvalue < 0.05 & Prob > 0.05)
message(">>> 显著通讯 (P<0.05): ", nrow(comm_sig), " 行")

# =========================================================================
# 5. 筛选 Cortex - Epidermis 通讯 (双向)
# =========================================================================
message(">>> 筛选 Cortex-Epidermis 通讯...")

# 标记 Source 和 Target 的组织类型
comm_sig <- comm_sig %>%
  mutate(
    source_tissue = gsub("_[0-9]+$", "", Source),
    target_tissue = gsub("_[0-9]+$", "", Target)
  )

# 筛选 Cortex↔Epidermis 通讯 (跨组织)
comm_ce <- comm_sig %>%
  filter(source_tissue != target_tissue) %>%
  filter(source_tissue %in% c("Cortex", "Epidermis"),
         target_tissue %in% c("Cortex", "Epidermis")) %>%
  mutate(Direction = paste0(source_tissue, "→", target_tissue))

message(">>> Cortex-Epidermis 显著通讯: ", nrow(comm_ce), " 行")
message(">>> 各方向通讯数:")
print(table(comm_ce$Direction))


# =========================================================================
# 6. 筛选 ROS 相关通讯
# =========================================================================
message(">>> 读取 ROS 基因列表...")
ros_genes_df <- read_excel(file.path(project_root, "data", "spatial_fish_36_genes.xlsx"))
ros_genes <- ros_genes_df %>%
  filter(Mod == "ROS") %>%
  pull(geneID)

message(">>> ROS 基因 (", length(ros_genes), " 个):")
print(ros_genes)

# 筛选 Ligand 或 Receptor 含 ROS 基因的通讯
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

# =========================================================================
# 7. 筛选全部 36 个 spatial_fish 基因相关通讯 (不仅 ROS)
# =========================================================================
all_spatial_genes <- ros_genes_df$geneID

comm_ce_spatial <- comm_ce %>%
  filter(sapply(seq_len(n()), function(i) {
    any(all_spatial_genes %in% c(Ligand[i], Receptor[i]))
  }))

message(">>> Cortex-Epidermis spatial_fish (36基因) 相关通讯: ", nrow(comm_ce_spatial), " 行")

# =========================================================================
# 8. 输出结果
# =========================================================================
message(">>> 保存结果...")

# ROS 通讯结果
if (nrow(comm_ce_ros) > 0) {
  write_csv(comm_ce_ros, "B73xOh40B_ROS_communications.csv")
  message("  已保存: B73xOh40B_ROS_communications.csv")
}

# 全部 spatial_fish 基因通讯结果
if (nrow(comm_ce_spatial) > 0) {
  write_csv(comm_ce_spatial, "B73xOh40B_spatial_fish_communications.csv")
  message("  已保存: B73xOh40B_spatial_fish_communications.csv")
}

# 全部 Cortex-Epidermis 显著通讯 (含非 ROS)
write_csv(comm_ce, "B73xOh40B_all_CE_communications.csv")
message("  已保存: B73xOh40B_all_CE_communications.csv")

# =========================================================================
# 9. 可视化 (简化版，仅 Interacting protein 类型)
# =========================================================================
message(">>> 绘图...")
dir.create("figures", showWarnings = FALSE)

nclusters <- length(unique(pccob@idents))
color_list <- rep(c('#1f7ab4', '#fd8d3c', '#66c1a7', '#e4191c',
                    '#a6cee5', '#bc80bd', '#ffed6f', 'grey60'),
                  length.out = nclusters)

# --- 图1: 总体通讯网络 ---
setwd(project_root)
# 保存 paracrine 网络为 svg
svglite("paracrine_network.svg", width = 8, height = 8)
PlottingCommunNetwork(pccob,
                      ligand.type   = "lrpairs",
                      comm.pattern  = "paracrine",
                      input.color   = color_list)
dev.off()

# 保存 autocrine 网络为 svg (chord 图)
svglite("autocrine_chord.svg", width = 8, height = 8)
PlottingCommunNetwork(pccob,
                      ligand.type   = "lrpairs",
                      comm.pattern  = "autocrine",
                      plot.type     = "chord",
                      input.color   = color_list)
dev.off()

# 保存 paracrine 网络为 svg (chord 图)
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

# 指定要加入的交互
extra_pairs <- c("Zm00001d030557_Zm00001d022160",
                 "Zm00001d051804_Zm00001d022160")

# 筛选并合并
extra_rows <- comm_ce[comm_ce$Interaction_name %in% extra_pairs, ]
df <- unique(rbind(df, extra_rows))



# ==========================================

# 3. 数据预处理
df_plot <- df %>%
  # 拼接 Source 和 Target，生成横坐标标签，例如 "Cortex_1 -> Epidermis_1"
  mutate(Source_Target = paste(Source, "->", Target)) %>%
  # 计算 -log10(Pvalue)，用于颜色映射
  mutate(neg_log10_Pvalue = -log10(Pvalue+0.000001))

# 4. 使用 ggplot2 绘制气泡图
p <- ggplot(df_plot, aes(x = Source_Target, y = Interaction_name)) +

  # 映射点的大小为 Prob，颜色为 -log10(Pvalue)
  geom_point(aes(size = Prob, color = neg_log10_Pvalue)) +

  # 设置颜色渐变 (你可以根据喜好修改颜色，这里用经典的 蓝-红 渐变)
  scale_color_gradient(low = "#4575b4", high = "#d73027", name = "-log10(P-value)") +

  # 设置点大小的图例名称
  scale_size_continuous(name = "Probability") +

  # 调整主题，使其看起来更干净专业
  theme_bw() +
  theme(
    # 横坐标标签倾斜 45 度，防止文字重叠
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    # 去除多余的网格线，突出数据点
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +

  # 设置坐标轴标题
  labs(
    x = "Cell-Cell Communication (Source -> Target)",
    y = "Interaction Name",
    title = "Receptor-Ligand Interactions"
  )

# 显示图表
svglite("interaction_dotplot.svg", width = 6, height = 8)
print(p)
dev.off()
