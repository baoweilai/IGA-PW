function compare_results()
%Compare baseline and refactored MAT data.

rootDir = fileparts(fileparts(mfilename('fullpath')));
baselineDir = fullfile(rootDir, 'paper_outputs', 'data', 'baseline');
refactoredDir = fullfile(rootDir, 'paper_outputs', 'data', 'refactored');
reportFile = fullfile(rootDir, 'paper_outputs', 'data', 'baseline_vs_refactored_report.txt');

assert(isfolder(baselineDir), 'Missing baseline directory: %s', baselineDir);
assert(isfolder(refactoredDir), 'Missing refactored directory: %s', refactoredDir);

baselineFiles = dir(fullfile(baselineDir, '*.mat'));
assert(~isempty(baselineFiles), 'No baseline MAT files found in %s.', baselineDir);
fid = fopen(reportFile, 'w');
assert(fid > 0, 'Cannot open comparison report: %s', reportFile);
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, 'Baseline versus refactored comparison\n');
fprintf(fid, 'Baseline:   %s\n', baselineDir);
fprintf(fid, 'Refactored: %s\n\n', refactoredDir);

tol = 1e-12;
for k = 1:numel(baselineFiles)
    oldFile = fullfile(baselineDir, baselineFiles(k).name);
    newName = strrep(baselineFiles(k).name, '_baseline.mat', '_refactored.mat');
    newFile = fullfile(refactoredDir, newName);
    assert(isfile(newFile), 'Missing refactored result: %s', newFile);

    oldData = load(oldFile);
    newData = load(newFile);
    compare_structs(oldData, newData, baselineFiles(k).name, tol, fid);
end

fprintf(fid, '\nAll compared numerical arrays satisfy relative tolerance %.1e.\n', tol);
end

function compare_structs(oldData, newData, label, tol, fid)
%Compare matching numeric fields in two structs.
oldNames = fieldnames(oldData);
for i = 1:numel(oldNames)
    name = oldNames{i};
    assert(isfield(newData, name), 'Missing field %s in refactored data for %s.', name, label);
    oldValue = oldData.(name);
    newValue = newData.(name);

    if isnumeric(oldValue) || islogical(oldValue)
        assert(isequal(size(oldValue), size(newValue)), ...
            'Size mismatch for %s in %s.', name, label);
        relErr = norm(double(newValue(:)) - double(oldValue(:))) / max(1, norm(double(oldValue(:))));
        assert(relErr < tol, 'Refactored result differs for %s in %s.', name, label);
        fprintf(fid, '%s :: %s :: relErr = %.3e\n', label, name, relErr);
    end
end
end
