# PubChem bulk downloader
# Downloads compounds in batches of 100 CIDs

$outputDir = "c:\Users\wmhor\dev\GuessWho\data\pubchem"
New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

$baseUrl = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
$props = "Title,MolecularFormula,CanonicalSMILES,MolecularWeight"
$batchSize = 100
$totalBatches = 500  # 500 batches * 100 = 50,000 compounds
$startCid = 1

Write-Host "Downloading $($totalBatches * $batchSize) compounds from PubChem..."
Write-Host "Output: $outputDir"
Write-Host ""

for ($batch = 0; $batch -lt $totalBatches; $batch++) {
    $startId = $startCid + ($batch * $batchSize)
    $endId = $startId + $batchSize - 1
    $cids = ($startId..$endId) -join ","
    $outFile = "$outputDir\batch_$($batch.ToString('D4')).csv"
    
    if (Test-Path $outFile) {
        Write-Host "[$($batch+1)/$totalBatches] Skipping (exists): $outFile"
        continue
    }
    
    try {
        $url = "$baseUrl/$cids/property/$props/CSV"
        Invoke-WebRequest -Uri $url -OutFile $outFile -TimeoutSec 60 -ErrorAction Stop
        $count = (Get-Content $outFile).Count - 1
        Write-Host "[$($batch+1)/$totalBatches] Downloaded $count compounds (CIDs $startId-$endId)"
        
        # Rate limiting - be nice to PubChem
        Start-Sleep -Milliseconds 200
    }
    catch {
        Write-Host "[$($batch+1)/$totalBatches] ERROR: $($_.Exception.Message)"
        # Continue to next batch
    }
    
    # Progress every 50 batches
    if (($batch + 1) % 50 -eq 0) {
        $files = Get-ChildItem "$outputDir\*.csv"
        $totalLines = ($files | ForEach-Object { (Get-Content $_).Count - 1 } | Measure-Object -Sum).Sum
        Write-Host "=== Progress: $totalLines compounds downloaded ==="
    }
}

# Final count
$files = Get-ChildItem "$outputDir\*.csv"
$totalLines = ($files | ForEach-Object { (Get-Content $_).Count - 1 } | Measure-Object -Sum).Sum
Write-Host ""
Write-Host "=== COMPLETE: $totalLines compounds in $($files.Count) files ==="
