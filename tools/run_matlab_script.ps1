param(
    [Parameter(Mandatory = $true)]
    [string]$ScriptFile,

    [Parameter(Mandatory = $true)]
    [string]$LogFile,

    [Parameter(Mandatory = $true)]
    [string]$ExitFile,

    [string]$MatlabExe = (Get-Command matlab.exe -ErrorAction Stop).Source
)

# Build the batch command for the supplied MATLAB script.
$matlabCommand = "run('$ScriptFile');"

# Run MATLAB and record its process exit code.
& $MatlabExe -batch $matlabCommand *> $LogFile
$exitCode = $LASTEXITCODE

Set-Content -LiteralPath $ExitFile -Value $exitCode -Encoding ASCII
exit $exitCode
