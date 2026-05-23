function audit_codebase()
%Write an inventory of the MATLAB codebase.

rootDir = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(rootDir, 'paper_outputs', 'data');
assert(isfolder(outputDir), 'Missing audit output directory.');

outFile = fullfile(outputDir, 'codebase_audit_before.txt');
folders = list_relative_folders(rootDir);
mFiles = list_relative_files(rootDir, '.m');

records = build_records(rootDir, mFiles);
duplicateFiles = find_duplicates({records.fileName}.', {records.path}.');
duplicateFunctions = find_duplicates({records.functionName}.', {records.path}.');
duplicateFunctions = duplicateFunctions(strlength(string({duplicateFunctions.name}.')) > 0);

fid = fopen(outFile, 'w');
assert(fid > 0, 'Unable to open audit output file.');

write_header(fid, rootDir, folders, mFiles);
write_section(fid, 'All folders', folders);
write_section(fid, 'All .m files', mFiles);
write_record_section(fid, 'Entry-point run scripts', records, [records.isEntryPoint]);
write_record_section(fid, 'Plotting scripts', records, [records.isPlotScript]);
write_record_section(fid, 'Table-generation scripts', records, [records.isTableScript]);
write_record_section(fid, 'Files that save figures', records, [records.savesFigure]);
write_record_section(fid, 'Files that save .mat data', records, [records.savesMatData]);
write_record_section(fid, 'Files containing nargin', records, [records.hasNargin]);
write_record_section(fid, 'Files containing try', records, [records.hasTry]);
write_record_section(fid, 'Files containing catch', records, [records.hasCatch]);
write_record_section(fid, 'Files containing fprintf', records, [records.hasFprintf]);
write_record_section(fid, 'Files containing disp', records, [records.hasDisp]);
write_record_section(fid, 'Files containing warning', records, [records.hasWarning]);
write_record_section(fid, 'Files containing Chinese comments', records, [records.hasChineseComment]);
write_duplicate_section(fid, 'Duplicated file names', duplicateFiles);
write_duplicate_section(fid, 'Duplicated function names', duplicateFunctions);
write_record_section(fid, 'Scripts that appear to generate paper figures', records, [records.appearsPaperFigure]);
write_record_section(fid, 'Scripts that appear to generate paper tables', records, [records.appearsPaperTable]);

fclose(fid);
end

function folders = list_relative_folders(rootDir)
%List folders below the project root.
items = dir(fullfile(rootDir, '**', '*'));
items = items([items.isdir]);
names = fullfile({items.folder}.', {items.name}.');
skip = ismember({items.name}.', {'.', '..'});
names = names(~skip);
folders = relative_paths(rootDir, names);
folders = sort(folders);
end

function files = list_relative_files(rootDir, extension)
%List files with the requested extension.
items = dir(fullfile(rootDir, '**', ['*' extension]));
names = fullfile({items.folder}.', {items.name}.');
files = relative_paths(rootDir, names);
files = sort(files);
end

function rel = relative_paths(rootDir, paths)
%Convert absolute paths to project-relative paths.
rootPrefix = [rootDir filesep];
rel = strings(numel(paths), 1);
for k = 1:numel(paths)
    rel(k) = string(erase(paths{k}, rootPrefix));
end
end

function records = build_records(rootDir, mFiles)
%Collect metadata for each MATLAB file.
records = repmat(empty_record(), numel(mFiles), 1);
for k = 1:numel(mFiles)
    relPath = char(mFiles(k));
    absPath = fullfile(rootDir, relPath);
    text = fileread(absPath);
    [~, fileStem, fileExt] = fileparts(relPath);
    lowerStem = lower(fileStem);
    lowerPath = lower(relPath);

    records(k).path = relPath;
    records(k).fileName = [fileStem fileExt];
    records(k).functionName = first_function_name(text);
    records(k).isFunctionFile = strlength(string(records(k).functionName)) > 0;
    records(k).isEntryPoint = startsWith(lowerStem, 'run') || startsWith(lowerStem, 'startup') ...
        || contains(lowerStem, 'main') || contains(lowerStem, 'driver');
    records(k).isPlotScript = startsWith(lowerStem, 'plot') || contains(lowerStem, '_plot') ...
        || contains(lowerPath, [filesep 'plot' filesep]) || contains(lowerPath, [filesep 'plots' filesep]);
    records(k).isTableScript = contains(lowerStem, 'table') || contains(lowerStem, 'tab') ...
        || contains(text, 'writetable') || contains(text, 'latex') || contains(text, '.tex');
    records(k).savesFigure = contains(text, 'exportgraphics') || contains(text, 'saveas') ...
        || contains(text, 'savefig') || contains(text, 'print(') || contains(text, '.pdf') ...
        || contains(text, '.eps') || contains(text, '.png') || contains(text, '.fig');
    records(k).savesMatData = ~isempty(regexp(text, '\<save\s*\(', 'once')) || contains(text, '.mat');
    records(k).hasNargin = contains(text, 'nargin');
    records(k).hasTry = contains_keyword(text, 'try');
    records(k).hasCatch = contains_keyword(text, 'catch');
    records(k).hasFprintf = contains(text, 'fprintf');
    records(k).hasDisp = contains(text, 'disp(') || contains(text, 'disp ');
    records(k).hasWarning = contains(text, 'warning(') || contains(text, 'warning ');
    records(k).hasChineseComment = has_chinese_comment(text);
    records(k).appearsPaperFigure = records(k).isPlotScript || records(k).savesFigure;
    records(k).appearsPaperTable = records(k).isTableScript;
end
end

function record = empty_record()
%Create one empty audit record.
record = struct( ...
    'path', '', ...
    'fileName', '', ...
    'functionName', '', ...
    'isFunctionFile', false, ...
    'isEntryPoint', false, ...
    'isPlotScript', false, ...
    'isTableScript', false, ...
    'savesFigure', false, ...
    'savesMatData', false, ...
    'hasNargin', false, ...
    'hasTry', false, ...
    'hasCatch', false, ...
    'hasFprintf', false, ...
    'hasDisp', false, ...
    'hasWarning', false, ...
    'hasChineseComment', false, ...
    'appearsPaperFigure', false, ...
    'appearsPaperTable', false);
end

function name = first_function_name(text)
%Read the first function name in a file.
tokens = regexp(text, '^\s*function\s+(?:\[[^\]]*\]\s*=\s*|[A-Za-z]\w*\s*=\s*)?([A-Za-z]\w*)', 'tokens', 'once', 'lineanchors');
name = '';
if ~isempty(tokens)
    name = tokens{1};
end
end

function result = contains_keyword(text, keyword)
%Check for a standalone MATLAB keyword.
result = ~isempty(regexp(text, ['(^|\n)\s*' keyword '(\s|$)'], 'once'));
end

function result = has_chinese_comment(text)
%Detect Chinese characters in comments.
lines = splitlines(string(text));
result = false;
for k = 1:numel(lines)
    line = char(lines(k));
    commentIndex = strfind(line, '%');
    if ~isempty(commentIndex)
        commentText = line(commentIndex(1):end);
        codePoints = double(commentText);
        if any(codePoints >= 19968 & codePoints <= 40959)
            result = true;
            return;
        end
    end
end
end

function duplicates = find_duplicates(values, paths)
%Locate an index or object used by the computation.
values = string(values);
paths = string(paths);
uniqueValues = unique(values);
duplicates = struct('name', {}, 'paths', {});
for k = 1:numel(uniqueValues)
    value = uniqueValues(k);
    if strlength(value) == 0
        continue;
    end
    mask = values == value;
    if nnz(mask) > 1
        duplicates(end + 1).name = char(value);
        duplicates(end).paths = cellstr(paths(mask));
    end
end
end

function write_header(fid, rootDir, folders, mFiles)
%Write the audit summary header.
fprintf(fid, 'IGA--PW--DG codebase audit before refactoring\n');
fprintf(fid, 'Root: %s\n', rootDir);
fprintf(fid, 'Generated: %s\n', datestr(now, 31));
fprintf(fid, 'Folder count: %d\n', numel(folders));
fprintf(fid, 'MATLAB file count: %d\n\n', numel(mFiles));
end

function write_section(fid, titleText, values)
%Write a text list section.
fprintf(fid, '============================================================\n');
fprintf(fid, '%s (%d)\n', titleText, numel(values));
fprintf(fid, '============================================================\n');
for k = 1:numel(values)
    fprintf(fid, '%s\n', values(k));
end
fprintf(fid, '\n');
end

function write_record_section(fid, titleText, records, mask)
%Write selected audit records.
selected = records(mask);
fprintf(fid, '============================================================\n');
fprintf(fid, '%s (%d)\n', titleText, numel(selected));
fprintf(fid, '============================================================\n');
for k = 1:numel(selected)
    fprintf(fid, '%s\n', selected(k).path);
end
fprintf(fid, '\n');
end

function write_duplicate_section(fid, titleText, duplicates)
%Write duplicate-name records.
fprintf(fid, '============================================================\n');
fprintf(fid, '%s (%d)\n', titleText, numel(duplicates));
fprintf(fid, '============================================================\n');
for k = 1:numel(duplicates)
    fprintf(fid, '%s\n', duplicates(k).name);
    for j = 1:numel(duplicates(k).paths)
        fprintf(fid, '  %s\n', duplicates(k).paths{j});
    end
end
fprintf(fid, '\n');
end
