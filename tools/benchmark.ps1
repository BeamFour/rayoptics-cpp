# Time this build of LensTool2 against the Java it was ported from, as described in
# PERFORMANCE.md. Runs one process at a time, C++ and Java alternating, so that drift
# in the machine affects both equally. Reports ratios; absolute times on a laptop are
# not comparable between sittings.
#
# usage:
#   tools\benchmark.ps1                       # plain run, 3 pairs
#   tools\benchmark.ps1 -Workload flags       # adds the optional report outputs
#   tools\benchmark.ps1 -Pairs 5 -Profile     # more pairs, plus the path-cache profile
#
# Close editors and browsers first, and run nothing else while it works.
#
# This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
# Copyright 2025-2026 by Dibyendu Majumdar
# License GPL v3
# See LICENSE-GPL-3.0.txt
param(
    [ValidateSet('plain', 'flags')][string] $Workload = 'plain',
    [int]    $Pairs   = 3,
    [switch] $Profile,
    [string] $Exe     = "$PSScriptRoot\..\build\Release\lenstool2.exe",
    [string] $Beam43  = "$PSScriptRoot\..\..\Beam43",
    [string] $Spec    = "$PSScriptRoot\..\Examples\canon-rf70-200mm-f2.8LZ\US20250155694_Example01P.txt",
    # Kept out of %TEMP%: Defender's scanning of newly written files there inflates
    # absolute times. cmake-build-* is already gitignored.
    [string] $WorkDir = "$PSScriptRoot\..\cmake-build-perf"
)

$ErrorActionPreference = 'Stop'
if (-not (Test-Path $Exe))  { throw "lenstool2 not found at $Exe -- build it first, or pass -Exe." }
if (-not (Test-Path $Spec)) { throw "prescription not found at $Spec" }

$extra = @()
if ($Workload -eq 'flags') {
    $extra = @('--output-ray-aberration-plots', '--output-wavelength-mtfs', '--do-wideangle-layout')
}

# The Java side is optional: without a built Beam43 only the C++ is timed.
$cp = "$Beam43\rayoptics\target\classes;$Beam43\mathlib\target\classes"
$withJava = (Test-Path "$Beam43\rayoptics\target\classes\org\redukti\tools\LensTool2.class")
if (-not $withJava) {
    Write-Host "No Beam43 classes found; timing the C++ only. Build them with: mvn -o compile"
}

# Every artifact is written next to the prescription unless --outdir says otherwise, so
# each side works on its own copy and the repository's Examples are left alone.
$name = Split-Path $Spec -Leaf
$dirs = @{ cpp = Join-Path $WorkDir 'cpp'; java = Join-Path $WorkDir 'java' }
foreach ($d in $dirs.Values) {
    if (Test-Path $d) { Remove-Item -Recurse -Force $d }
    New-Item -ItemType Directory -Force $d | Out-Null
    Copy-Item $Spec $d
}
$specC = Join-Path $dirs.cpp $name
$specJ = Join-Path $dirs.java $name

# Nothing left over from an interrupted run may compete with this one.
Get-Process lenstool2 -ErrorAction SilentlyContinue | Stop-Process -Force
Get-CimInstance Win32_Process -Filter "Name='java.exe'" |
    Where-Object { $_.CommandLine -like '*LensTool2*' } |
    ForEach-Object { Stop-Process -Id $_.ProcessId -Force }

Write-Host "workload $Workload, $Pairs pair(s), work dir $WorkDir"
$ratios = @()
for ($i = 1; $i -le $Pairs; $i++) {
    $c = Measure-Command { & $Exe --specfile $specC @extra 2>&1 | Out-Null }
    if ($withJava) {
        $j = Measure-Command {
            & java -cp $cp org.redukti.tools.LensTool2 --specfile $specJ @extra 2>&1 | Out-Null
        }
        $ratio = $c.TotalSeconds / $j.TotalSeconds
        $ratios += $ratio
        "pair $i  C++/Java {0:N3}" -f $ratio | Write-Host
    }
    else {
        "pair $i  cpp {0:N2} s" -f $c.TotalSeconds | Write-Host
    }
}
if ($ratios.Count) {
    "mean C++/Java {0:N3}  (below 1 means the C++ is faster)" -f ($ratios | Measure-Object -Average).Average |
        Write-Host
}

if ($Profile) {
    # The path-cache counters go to stderr at normal exit.
    $env:RAYOPTICS_PROFILE_PATH = '1'
    $prof = Join-Path $WorkDir 'profile.txt'
    & $Exe --specfile $specC @extra 2> $prof | Out-Null
    Remove-Item Env:RAYOPTICS_PROFILE_PATH
    Select-String -Path $prof -Pattern 'path profile|hit_rate' | ForEach-Object { $_.Line } | Write-Host
}

$counts = $dirs.GetEnumerator() | ForEach-Object {
    $n = (Get-ChildItem $_.Value -File | Where-Object { $_.Name -ne $name }).Count
    "$($_.Key)=$n"
}
Write-Host ("artifacts written: " + ($counts -join ' '))

# Which artifacts the two versions disagree on: the documented hexapolar sin/cos
# divergence reaches the MTF plots and the spot reports.
if ($withJava) {
    $differ = @()
    foreach ($f in Get-ChildItem $dirs.cpp -File | Where-Object { $_.Name -ne $name }) {
        $other = Join-Path $dirs.java $f.Name
        if (-not (Test-Path $other)) { $differ += "$($f.Name) (missing in java)"; continue }
        $a = (Get-Content $f.FullName -Raw) -replace "`r`n", "`n"
        $b = (Get-Content $other -Raw) -replace "`r`n", "`n"
        if ($a -ne $b) { $differ += $f.Name }
    }
    Write-Host "artifacts differing: $($differ.Count)"
    if ($differ.Count) { $differ | ForEach-Object { Write-Host "    $_" } }
}
