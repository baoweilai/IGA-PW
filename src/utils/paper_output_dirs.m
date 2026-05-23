function dirs = paper_output_dirs(exampleName)
%Create output folders for one example.

assert(ischar(exampleName) || isstring(exampleName), ...
    'paper_output_dirs requires an example name.');

exampleName = char(exampleName);
rootDir = project_root();

dirs.root = fullfile(rootDir, 'paper_outputs');
dirs.figures = fullfile(dirs.root, 'figures', exampleName);
dirs.tables = fullfile(dirs.root, 'tables', exampleName);
dirs.data = fullfile(dirs.root, 'data', exampleName);

ensure_directory(dirs.figures);
ensure_directory(dirs.tables);
ensure_directory(dirs.data);
end
