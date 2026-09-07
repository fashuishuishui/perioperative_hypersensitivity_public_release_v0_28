[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$CanadaDir,

    [Parameter(Mandatory = $true)]
    [string]$JaderDir,

    [string]$PythonExe = "python",

    [string]$RscriptExe = "Rscript",

    [string]$FaersZipDir,

    [string]$MedDraAscii,

    [string]$MedDraReferenceAscii,

    [switch]$RebuildFaersRaw,

    [switch]$ReuseValidatedFaersRawIntegrity,

    [switch]$SkipLicensedSmq
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$RawBuild = Join-Path $Root "raw_rebuild"
$Work = Join-Path $Root "work"
$FaersBuild = Join-Path $Work "faers_raw_2004q1_2025q4_schema_fixed"
$Output = Join-Path $Root "aggregate_outputs\schema_fixed"
$DefaultFaersZipDir = Join-Path $Root "raw_inputs\faers_ascii_2004q1_2025q4_reretrieved_20260730"
if (-not $FaersZipDir) {
    $FaersZipDir = $DefaultFaersZipDir
}

foreach ($Path in @($CanadaDir, $JaderDir)) {
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "Required source directory was not found: $Path"
    }
}
if (-not (Test-Path -LiteralPath $FaersZipDir)) {
    throw "FAERS ZIP directory was not found: $FaersZipDir"
}

foreach ($Command in @($PythonExe, $RscriptExe)) {
    if (-not (Get-Command $Command -ErrorAction SilentlyContinue) -and -not (Test-Path -LiteralPath $Command)) {
        throw "Executable was not found: $Command"
    }
}

$env:PV_REPO_ROOT = $Root
$env:PV_FAERS_RAW_DB = Join-Path $FaersBuild "faers_raw_2004q1_2025q4.sqlite"
$env:PV_FAERS_CACHE_DB = Join-Path $FaersBuild "faers_latest_2004q1_2025q4.sqlite"
$env:PV_FAERS_FLAGS_DB = Join-Path $Work "faers_flags_schema_fixed.sqlite"
$env:PV_REPORT_FLAGS_DB = Join-Path $Work "report_level_flags_canada_jader_schema_fixed.sqlite"
$env:PV_CANADA_JADER_FLAGS_DB = $env:PV_REPORT_FLAGS_DB
$env:PV_CANADA_DIR = (Resolve-Path -LiteralPath $CanadaDir).Path
$env:PV_JADER_DIR = (Resolve-Path -LiteralPath $JaderDir).Path
$env:PV_ANALYSIS_OUTPUT_DIR = $Output
$env:PV_AGGREGATE_OUTPUT_DIR = $Output
$env:PV_ESTIMAND_OUTPUT_DIR = $Output
$env:PV_R_LIBS_USER = Join-Path $Root ".r_libs"

if ($RebuildFaersRaw) {
    & $PythonExe (Join-Path $RawBuild "build_faers_raw_branch.py") `
        --manifest (Join-Path $RawBuild "faers_ascii_manifest_2004q1_2025q4.csv") `
        --zip-dir $FaersZipDir `
        --output-dir $FaersBuild `
        --force
} else {
    $RawReuseArgs = @(
        (Join-Path $RawBuild "build_faers_raw_branch.py"),
        "--manifest", (Join-Path $RawBuild "faers_ascii_manifest_2004q1_2025q4.csv"),
        "--zip-dir", $FaersZipDir,
        "--output-dir", $FaersBuild,
        "--reuse-raw",
        "--force"
    )
    if ($ReuseValidatedFaersRawIntegrity) {
        $RawReuseArgs += "--skip-raw-integrity-check"
    }
    & $PythonExe @RawReuseArgs
}

& $PythonExe (Join-Path $Root "code\build_faers_flags.py")
& $PythonExe (Join-Path $Root "code\build_report_level_flags.py")
& $PythonExe (Join-Path $Root "code\run_estimand_context_sensitivities.py")
& $PythonExe (Join-Path $Root "code\run_reference_setting_audits.py")
& $PythonExe (Join-Path $Root "code\run_country_role_probes.py")
& $PythonExe (Join-Path $Root "code\run_broad_pt_sensitivity.py")
& $PythonExe (Join-Path $Root "code\run_teicoplanin_time_audit.py")
& $PythonExe (Join-Path $Root "code\run_jader_structure_audit.py")

if (-not $SkipLicensedSmq) {
    if (-not $MedDraAscii) {
        throw "Set -MedDraAscii or use -SkipLicensedSmq. Licensed SMQ files are required for this optional sensitivity."
    }
    $env:PV_MEDDRA_ASCII = (Resolve-Path -LiteralPath $MedDraAscii).Path
    $env:PV_MEDDRA_REFERENCE_ASCII = if ($MedDraReferenceAscii) {
        (Resolve-Path -LiteralPath $MedDraReferenceAscii).Path
    } else {
        $env:PV_MEDDRA_ASCII
    }
    $CanadaReactionFile = Get-ChildItem -LiteralPath $env:PV_CANADA_DIR -File -Filter "reactions.txt*" | Select-Object -First 1
    if (-not $CanadaReactionFile) {
        throw "No Canada reactions file matching reactions.txt* was found in $($env:PV_CANADA_DIR)"
    }
    $env:PV_CANADA_REACTIONS = $CanadaReactionFile.FullName
    & $PythonExe (Join-Path $Root "code\run_official_smq_sensitivity.py")
}

& $PythonExe (Join-Path $Root "code\build_v0_26_canonical_artifacts.py")
& $RscriptExe (Join-Path $Root "code\validate_breslow_day.R")
& $RscriptExe (Join-Path $Root "code\generate_figures_v0_26.R")
& $PythonExe (Join-Path $Root "code\build_supplementary_v0_26.py")

Write-Host "v0.26 rebuild completed. Check aggregate_outputs/schema_fixed/canonical_artifact_lock_v0_26.json for the validation result."
