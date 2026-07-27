function ensure_directory(folderPath)
% Create a directory when it is missing.

if ~isfolder(folderPath)
    mkdir(folderPath);
end
end
