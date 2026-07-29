project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

rm(list=ls())
gc()
setwd(project_root)
mmd_summary <- read_csv("data_output/mmd_result/mmd_summary.csv")
library(tidyverse)
library(rstatix)

df_long <- mmd_summary %>%
  pivot_longer(cols = everything(),
               names_to = "comparison",
               values_to = "mmd")

df_long$comparison <- factor(df_long$comparison,
                             levels = colnames(mmd_summary))

df_long$pair <- sub("^[^_]*_([^_]*)_.*$", "\\1", df_long$comparison)

df_long %>%
  group_by(pair) %>%
  t_test(mmd ~ comparison)

dev.off()

png("./Fig_output/scmmd.png",height = 1500,width = 1200,res=400)
ggplot(df_long, aes(x = comparison, y = mmd,fill = pair)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~pair, scales = "free_x") +
  scale_fill_manual(values = c(
    "#A6CEE3",
    "#B2DF8A",
    "#FDBF6F"
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, colour = "black"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.title.y = element_text(colour = "black",size = 14),
    axis.text.y = element_text(colour = "black",size = 9)
  ) +
  labs(
    x = "",
    y = "MMD",
    title = " "
  )+NoLegend()
dev.off()
