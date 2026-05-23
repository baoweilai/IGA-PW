function runFile = dg_case_file(rootDir, studyName, Nc, pdeg, nElem)
%Return the DG error data filename.
arguments
    rootDir
    studyName
    Nc
    pdeg
    nElem
end
runFile = fullfile(rootDir, 'result', upper(studyName), ...
    sprintf('K_%d', Nc), sprintf('p_%d', pdeg), ...
    sprintf('nelem_%02d', nElem), 'run.mat');
end
