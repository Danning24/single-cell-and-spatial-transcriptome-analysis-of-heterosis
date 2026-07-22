# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

# =========================================================================
# PlantCellChat v0.0.9 Bug 热修复补丁
# Bug: ExtractExpMatrix() 中 interaction[, comps] 对 tibble 返回单列 tibble
#      导致 match() 只返回 1 个 NA，表达值全被替换为 1，Prob 统一为 0.667
# Fix: [, comps] → [[comps]]（总是返回向量，data.frame 和 tibble 通用）
#
# 用法：在 library(PlantCellChat) 之后 source 本文件即可
# =========================================================================

patch_ExtractExpMatrix <- function() {
  ns <- asNamespace("PlantCellChat")
  func_text <- deparse(get("ExtractExpMatrix", envir = ns))

  # 检查是否已经修复
  if (any(grepl("interaction\\[\\[comps\\]\\]", func_text))) {
    message("PlantCellChat 已经打过补丁，跳过。")
    return(invisible())
  }

  # 替换: interaction[, comps] → interaction[[comps]]
  func_text <- gsub(
    "index <- interaction[, comps]",
    "index <- interaction[[comps]]",
    func_text,
    fixed = TRUE
  )

  new_func <- eval(parse(text = func_text))
  environment(new_func) <- ns

  unlockBinding("ExtractExpMatrix", ns)
  assign("ExtractExpMatrix", new_func, envir = ns)
  lockBinding("ExtractExpMatrix", ns)

  message("✓ PlantCellChat ExtractExpMatrix bug 已热修复")
}

patch_ExtractExpMatrix()
