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

    [string]$MatlabExe = (Get-Command matlab.exe -ErrorAction Stop).Source
)

# Build the batch command for one example workflow.
$matlabCommand = "cd('$RepoRoot'); addpath(fullfile(pwd,'examples','$Example')); $Runner('$Workflow');"

# Run MATLAB and record its process exit code.
& $MatlabExe -batch $matlabCommand *> $LogFile
$exitCode = $LASTEXITCODE

Set-Content -LiteralPath $ExitFile -Value $exitCode -Encoding ASCII
exit $exitCode
