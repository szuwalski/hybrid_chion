# Manuscript

Edit `manuscript.Rmd` and click **Knit to PDF** in RStudio, or run this from the
repository root:

```powershell
./manuscript/render_pdf.ps1
```

The command writes `output/pdf/manuscript.pdf`. RStudio's Knit button writes
`manuscript/manuscript.pdf` by default. Both use the same R Markdown source.

The PDF setup follows `full_manuscript_working_copy.Rmd`: R Markdown's
`pdf_document`, a table of contents, left-side line numbers, and the same float
and landscape packages. The example is unchanged. Its unrelated analysis code
and missing CSL reference are not copied into the new manuscript.

The renderer uses the existing R, rmarkdown, RStudio-bundled Pandoc, and LaTeX
installation. Automatic TinyTeX package installation is disabled, including
when knitting in RStudio. Missing dependencies produce an error rather than
triggering an installation.

The original Markdown draft is now `manuscript.Rmd` so there is only one editable
manuscript source. Equations and section placeholders are retained, and HTML
drafting comments do not appear in the PDF. The existing Word file is unchanged.
