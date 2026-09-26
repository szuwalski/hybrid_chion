param([string]$RscriptPath = "")

$ErrorActionPreference = "Stop"
if (-not $RscriptPath) {
    $rscriptCommand = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($rscriptCommand) {
        $RscriptPath = $rscriptCommand.Source
    } else {
        $RscriptPath = Join-Path $env:ProgramFiles "R/R-4.4.2/bin/Rscript.exe"
    }
}
if (-not (Test-Path -LiteralPath $RscriptPath)) {
    throw "Provide -RscriptPath pointing to your existing Rscript.exe."
}

& $RscriptPath (Join-Path $PSScriptRoot "render_pdf.R")
if ($LASTEXITCODE -ne 0) {
    throw "PDF rendering failed with exit code $LASTEXITCODE."
}
