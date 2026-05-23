function addedPaths = add_workflow_paths(workflowDir, expectedSubdirs)
%Add model subfolders to the MATLAB path.

assert(isfolder(workflowDir), 'Missing workflow directory: %s', workflowDir);

addedPaths = {workflowDir};
for k = 1:numel(expectedSubdirs)
    folderPath = fullfile(workflowDir, expectedSubdirs{k});
    assert(isfolder(folderPath), 'Missing workflow dependency: %s', folderPath);
    addedPaths{end + 1} = folderPath;
end

for k = 1:numel(addedPaths)
    addpath(addedPaths{k}, '-begin');
end
end
