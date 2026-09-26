options(tinytex.install_packages = FALSE)

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
manuscript_directory <- dirname(script_path)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The existing R installation must include rmarkdown; nothing will be installed.")
}
if (!rmarkdown::pandoc_available()) {
  bundled_pandoc <- file.path(Sys.getenv("ProgramFiles"),
                             "RStudio/resources/app/bin/quarto/bin/tools")
  if (file.exists(file.path(bundled_pandoc, "pandoc.exe"))) {
    Sys.setenv(RSTUDIO_PANDOC = bundled_pandoc)
  }
}

rmarkdown::render(
  file.path(manuscript_directory, "manuscript.Rmd"),
  output_dir = file.path(manuscript_directory, "../output/pdf"),
  envir = new.env(parent = globalenv())
)
