function workflowDir = activate_example_workflow(workflowName, expectedSubdirs)
% Activate paths for one example task.

stack = dbstack('-completenames');
assert(numel(stack) >= 2, 'activate_example_workflow must be called from an example script.');

callerDir = fileparts(stack(2).file);
workflowDir = '';

[parentDir, leafName] = fileparts(callerDir);
if ismember(leafName, {'figures', 'tables', 'data'})
    exampleDir = parentDir;
else
    [modelDir, workflowLeaf] = fileparts(callerDir);
    [exampleDirCandidate, modelLeaf] = fileparts(modelDir);
    if strcmp(modelLeaf, 'model') && strcmp(workflowLeaf, workflowName)
        workflowDir = callerDir;
        exampleDir = exampleDirCandidate;
    else
        exampleDir = callerDir;
    end
end

if isempty(workflowDir)
    workflowDir = fullfile(exampleDir, 'model', workflowName);
end
add_workflow_paths(workflowDir, expectedSubdirs);

commonDir = fullfile(exampleDir, 'model', 'common');
if isfolder(commonDir)
    add_workflow_paths(commonDir, {'operators'});
end
end
