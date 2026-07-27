function [summary, csvPath] = build_comparison_table()
% BUILD_COMPARISON_TABLE Build the eight-row comparison table.

exampleDir = fileparts(fileparts(mfilename('fullpath')));
resultDir = fullfile(exampleDir, 'data', 'result', 'hartree_comparison');

comparison = readtable(fullfile(resultDir, 'comparison.csv'), ...
    'TextType', 'string');
assert(height(comparison) == 8, ...
    'The comparison table must contain eight rows.');
summary = comparison(:, ...
    {'method', 'parameters', 'nDof', 'energyError', 'totalTime'});

csvPath = fullfile(resultDir, 'table.csv');
writetable(summary, csvPath);
end
