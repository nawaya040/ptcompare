$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
$runtime = 'C:/Program Files/R/R-4.5.2/bin/x64/Rterm.exe'
$env:LC_ALL = 'C'
$cases = @(
  @(1000,2,64,'range'), @(1000,10,64,'range'), @(1000,50,64,'range'), @(1000,100,64,'range'),
  @(10000,10,64,'range'), @(10000,20,64,'range'), @(50000,20,32,'range'), @(200000,20,32,'range'),
  @(2000,10,64,'range'), @(2000,10,256,'range'), @(2000,10,1024,'range'),
  @(1000,2,64,'rank'), @(1000,10,64,'rank'), @(1000,50,64,'rank'), @(1000,100,64,'rank'),
  @(10000,20,64,'rank'), @(50000,20,32,'rank'), @(200000,20,32,'rank')
)
$results = @()
for ($caseId=0; $caseId -lt $cases.Count; $caseId++) {
  $case = $cases[$caseId]
  $relativeOutput = "evaluation/results/benchmark-$caseId.csv"
  $arguments = @('--vanilla','--slave','--file=evaluation/benchmark-one.R','--args', $case[0], $case[1], $case[2], $relativeOutput, $case[3])
  $process = Start-Process -FilePath $runtime -ArgumentList $arguments -WorkingDirectory $projectRoot -WindowStyle Hidden -PassThru
  $peakBytes = 0L
  while (-not $process.HasExited) {
    $process.Refresh()
    try { $peakBytes = [Math]::Max($peakBytes, $process.PeakWorkingSet64) } catch {}
    Start-Sleep -Milliseconds 100
  }
  $process.WaitForExit()
  if ($process.ExitCode -ne 0) { throw "Benchmark $caseId failed" }
  $row = Import-Csv -LiteralPath (Join-Path $projectRoot $relativeOutput)
  $row | Add-Member -NotePropertyName peak_working_set_mb -NotePropertyValue ($peakBytes/1MB)
  $results += $row
  $results | Export-Csv -LiteralPath (Join-Path $projectRoot 'evaluation/results/benchmarks.csv') -NoTypeInformation
  Write-Output "Benchmark $caseId completed: n=$($case[0]), d=$($case[1]), M=$($case[2])"
}
