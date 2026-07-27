function mexFile = build_rfp_mex(cfg)

% Compile the RFP sparse-factorization MEX function.
arguments
    cfg struct
end

sourceFile = fullfile(cfg.solverDir, 'rfp_chol_mex.c');
assert(isfile(sourceFile), 'Missing RFP MEX source: %s', sourceFile);

mex('-R2018a', ...
    '-outdir', cfg.solverDir, ...
    '-output', 'rfp_chol_mex', ...
    sourceFile, '-lmwlapack', '-lmwblas');

mexFile = fullfile(cfg.solverDir, ['rfp_chol_mex.', mexext]);
assert(isfile(mexFile), 'RFP MEX build did not create: %s', mexFile);
end
