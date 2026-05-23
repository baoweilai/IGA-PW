param(
    [Parameter(Mandatory = $true)]
    [string]$RepoRoot,

    [Parameter(Mandatory = $true)]
    [string]$Example,

    [Parameter(Mandatory = $true)]
    [string]$Runner,

    [Parameter(Mandatory = $true)]
    [string]$Workflow,

    [Parameter(Mandatory = $true)]
    [string]$LogFile,

    [Parameter(Mandatory = $true)]
    [string]$ExitFile,

    [string]$MatlabExe = 'C:\Matlab2025a\bin\matlab.exe'
)

$matlabCommand = "cd('$RepoRoot'); addpath(fullfile(pwd,'examples','$Example')); $Runner('$Workflow');"

& $MatlabExe -batch $matlabCommand *> $LogFile
$exitCode = $LASTEXITCODE

Set-Content -LiteralPath $ExitFile -Value $exitCode -Encoding ASCII
exit $exitCode
