project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

Cortex$logp<- -log10(as.numeric(Cortex$`p...8`))
library(ggplot2)
library(ggrepel)

ggplot(Cortex, aes(x = `r...4`,
                   y = `r...7`,
                   size = n_traits)) +
  geom_point(shape = 21,
             fill = "red",
             color = "black") +
  scale_size_continuous(range = c(1, 6),
                        breaks = pretty(Cortex$n_traits, n = 3)) +
  theme_bw() +
  geom_text_repel(aes(label = hyb),
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "grey50")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

labs(
  x = "Pearson’s r(Nonaddgene expression level and E_Score)",
  y = "Pearson’s r(Nonaddgene expression level and M_Score)",
  size = "n_traits",
  title = "Association between Nonaddgene expression and development"
)
