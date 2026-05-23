function E = cached_dg_error(cacheFile, refFile, runFile, opt)
%Load or compute cached DG error data.
refStamp = file_stamp_local(refFile);
runStamp = file_stamp_local(runFile);
optKey = opt_key_local(opt);

if exist(cacheFile, 'file')
    S = load(cacheFile, 'E');
    if isfield(S, 'E') && is_valid_cache_local(S.E, refFile, runFile, refStamp, runStamp, optKey)
        E = S.E;
        return;
    end
end

Sref = load(refFile, 'run');
Scase = load(runFile, 'run');
E = struct();
E.referenceRunFile = refFile;
E.caseRunFile = runFile;
E.referenceStamp = refStamp;
E.caseStamp = runStamp;
E.optKey = optKey;
E.lambdaRef = real(Sref.run.lambda(1));
E.lambdaCase = real(Scase.run.lambda(1));
E.errLambda = max(abs(E.lambdaRef - E.lambdaCase), eps);
E.errDG = max(dg_error(refFile, runFile, opt), eps);

cacheDir = fileparts(cacheFile);
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
save(cacheFile, 'E');
end

function tf = is_valid_cache_local(E, refFile, runFile, refStamp, runStamp, optKey)
%Check whether cached error data is current.
tf = isfield(E, 'referenceRunFile') && isfield(E, 'caseRunFile') && ...
    isfield(E, 'referenceStamp') && isfield(E, 'caseStamp') && isfield(E, 'optKey') && ...
    strcmp(E.referenceRunFile, refFile) && strcmp(E.caseRunFile, runFile) && ...
    E.referenceStamp == refStamp && E.caseStamp == runStamp && strcmp(E.optKey, optKey) && ...
    isfield(E, 'errLambda') && isfield(E, 'errDG');
end

function stamp = file_stamp_local(fileName)
%Return a file timestamp key.
d = dir(fileName);
if isempty(d), error('Missing run file: %s', fileName); end
stamp = d.datenum + d.bytes * eps;
end

function key = opt_key_local(opt)
%Build a cache key from options.
fields = {'innerGridN', 'outerGridN', 'faceGridN', 'chunkSize', 'Csigma'};
parts = cell(1, numel(fields));
for i = 1:numel(fields)
    if isfield(opt, fields{i})
        parts{i} = sprintf('%s=%.15g', fields{i}, opt.(fields{i}));
    else
        parts{i} = sprintf('%s=[]', fields{i});
    end
end
key = strjoin(parts, ';');
end
