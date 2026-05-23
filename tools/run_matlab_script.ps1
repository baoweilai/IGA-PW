param(
    [Parameter(Mandatory = $true)]
    [string]$ScriptFile,

    [Parameter(Mandatory = $true)]
    [string]$LogFile,

    [Parameter(Mandatory = $true)]
    [string]$ExitFile,

    [string]$MatlabExe = 'C:\Matlab2025a\bin\matlab.exe'
)

$matlabCommand = "run('$ScriptFile');"

& $MatlabExe -batch $matlabCommand *> $LogFile
$exitCode = $LASTEXITCODE

Set-Content -LiteralPath $ExitFile -Value $exitCode -Encoding ASCII
exit $exitCode
