function copiedFiles = sync_paper_outputs(exampleDir, exampleName, workflowName, beforeSnapshot)
%Copy new outputs into the shared output tree.

assert(isfolder(exampleDir), 'Missing example directory: %s', exampleDir);
assert(isstruct(beforeSnapshot), 'sync_paper_outputs requires a snapshot.');

workflowName = char(workflowName);
dirs = paper_output_dirs(exampleName);
afterSnapshot = snapshot_example_outputs(exampleDir);
changedFiles = changed_outputs_local(afterSnapshot, beforeSnapshot);

copiedFiles = strings(0, 1);
for k = 1:numel(changedFiles)
    kind = output_kind_local(changedFiles(k).path);
    if kind == ""
        continue;
    end
    relPath = relative_output_path_local(exampleDir, changedFiles(k).path, workflowName);
    destFile = fullfile(dirs.(char(kind)), workflowName, relPath);
    ensure_directory(fileparts(destFile));
    copyfile(changedFiles(k).path, destFile, 'f');
    remove_local_copy_local(changedFiles(k).path, kind);
    copiedFiles(end + 1, 1) = string(destFile); %#ok<AGROW>
end
end

function changedFiles = changed_outputs_local(afterSnapshot, beforeSnapshot)
%Keep files created or updated by this task.

changedFiles = struct('path', {}, 'datenum', {}, 'bytes', {});
beforePaths = {beforeSnapshot.path};
for k = 1:numel(afterSnapshot)
    idx = find(strcmp(beforePaths, afterSnapshot(k).path), 1);
    if isempty(idx) || afterSnapshot(k).datenum ~= beforeSnapshot(idx).datenum ...
            || afterSnapshot(k).bytes ~= beforeSnapshot(idx).bytes
        changedFiles(end + 1) = afterSnapshot(k); %#ok<AGROW>
    end
end
end

function kind = output_kind_local(filePath)
%Classify a generated output file.

[~, ~, ext] = fileparts(filePath);
switch lower(ext)
    case {'.png', '.pdf', '.fig', '.eps', '.svg'}
        kind = "figures";
    case {'.tex', '.xlsx'}
        kind = "tables";
    case {'.csv', '.mat', '.txt', '.md'}
        kind = "data";
    otherwise
        kind = "";
end
end

function relPath = relative_output_path_local(exampleDir, filePath, workflowName)
%Keep only the result-specific branch in paper_outputs.

roots = {
    fullfile(exampleDir, 'data', 'result')
    fullfile(exampleDir, 'data')
    fullfile(exampleDir, 'results')
    fullfile(exampleDir, 'figures')
    };

for k = 1:numel(roots)
    rootPath = roots{k};
    prefix = [rootPath filesep];
    if strncmpi(filePath, prefix, numel(prefix))
        relPath = filePath(numel(prefix) + 1:end);
        if k == 1
            relPath = remove_example_label_local(relPath);
        end
        relPath = remove_workflow_label_local(relPath, workflowName);
        return;
    end
end

[~, name, ext] = fileparts(filePath);
relPath = [name ext];
end

function relPath = remove_example_label_local(relPath)
%Drop the repeated Example_i folder name.

parts = regexp(relPath, regexptranslate('escape', filesep), 'split');
if numel(parts) > 1 && ~isempty(regexp(parts{1}, '^Example_\d+$', 'once'))
    relPath = fullfile(parts{2:end});
end
end

function relPath = remove_workflow_label_local(relPath, workflowName)
%Drop a repeated workflow folder name.

parts = regexp(relPath, regexptranslate('escape', filesep), 'split');
if numel(parts) > 1 && strcmp(parts{1}, char(workflowName))
    relPath = fullfile(parts{2:end});
end
end

function remove_local_copy_local(filePath, kind)
%Keep final figures and tables only under paper_outputs.

if kind == "figures" || kind == "tables"
    delete(filePath);
end
end
