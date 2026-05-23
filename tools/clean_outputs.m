function clean_outputs()
%Remove generated figures, tables, and refactored data.

rootDir = fileparts(fileparts(mfilename('fullpath')));
outputDirs = {
    fullfile(rootDir, 'paper_outputs', 'figures')
    fullfile(rootDir, 'paper_outputs', 'tables')
    fullfile(rootDir, 'paper_outputs', 'data', 'refactored')
    };

for k = 1:numel(outputDirs)
    assert(isfolder(fileparts(outputDirs{k})), 'Missing output parent: %s', fileparts(outputDirs{k}));
    if isfolder(outputDirs{k})
        delete(fullfile(outputDirs{k}, '*'));
    else
        mkdir(outputDirs{k});
    end
end
end
