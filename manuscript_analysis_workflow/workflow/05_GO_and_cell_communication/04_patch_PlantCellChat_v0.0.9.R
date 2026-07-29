project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

patch_ExtractExpMatrix <- function() {
  ns <- asNamespace("PlantCellChat")
  func_text <- deparse(get("ExtractExpMatrix", envir = ns))

  if (any(grepl("interaction\\[\\[comps\\]\\]", func_text))) {
    message("PlantCellChat 已经打过补丁，跳过。")
    return(invisible())
  }

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
