function rootDir = project_root()
% Return the repository root directory.

thisFile = mfilename('fullpath');
rootDir = fileparts(fileparts(fileparts(thisFile)));
end
