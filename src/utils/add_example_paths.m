function addedPaths = add_example_paths(exampleDir)
%Add one example folder to the MATLAB path.

requiredDirs = {
    exampleDir
    fullfile(exampleDir, 'figures')
    fullfile(exampleDir, 'tables')
    fullfile(exampleDir, 'data')
    };

modelDir = fullfile(exampleDir, 'model');
if isfolder(modelDir)
    requiredDirs{end + 1} = modelDir;
end

for k = 1:numel(requiredDirs)
    assert(isfolder(requiredDirs{k}), 'Missing example directory: %s', requiredDirs{k});
end

addedPaths = requiredDirs(:);
for k = 1:numel(addedPaths)
    addpath(addedPaths{k}, '-begin');
end
end
