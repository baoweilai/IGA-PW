function snapshot = snapshot_example_outputs(exampleDir)
%Record generated files for one example.

assert(isfolder(exampleDir), 'Missing example directory: %s', exampleDir);

scanRoots = {
    fullfile(exampleDir, 'data')
    fullfile(exampleDir, 'results')
    fullfile(exampleDir, 'figures')
    };

snapshot = struct('path', {}, 'datenum', {}, 'bytes', {});
for ir = 1:numel(scanRoots)
    if ~isfolder(scanRoots{ir})
        continue;
    end
    files = dir(fullfile(scanRoots{ir}, '**', '*'));
    for k = 1:numel(files)
        if files(k).isdir
            continue;
        end
        filePath = fullfile(files(k).folder, files(k).name);
        if ~is_candidate_output_file_local(filePath)
            continue;
        end
        snapshot(end + 1).path = filePath; %#ok<AGROW>
        snapshot(end).datenum = files(k).datenum;
        snapshot(end).bytes = files(k).bytes;
    end
end
end

function tf = is_candidate_output_file_local(filePath)
%Select generated output files.

[~, fileName, ext] = fileparts(filePath);
ext = lower(ext);
filePathLower = lower(filePath);

allowedExt = {'.png', '.pdf', '.fig', '.eps', '.svg', ...
    '.csv', '.mat', '.tex', '.txt', '.md', '.xlsx'};
tf = ismember(ext, allowedExt);
if ~tf
    return;
end

if strcmpi(fileName, 'run')
    tf = false;
    return;
end

tf = ~contains(filePathLower, [filesep 'cache']);
end
