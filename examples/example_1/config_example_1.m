function cfg = config_example_1()
% Return fixed Example 1 parameters.

cfg.example = 1;
cfg.name = 'Example 1';
cfg.outputPrefix = 'ex1';
cfg.figureScripts = {
    'plot_h_convergence'
    'plot_cutoff_convergence'
    'plot_eigen'
    'plot_scaled_errors'
    'plot_penalty'
    'plot_method_comparison'
    };
cfg.dataScripts = {
    'model/h_convergence/run_h_data.m'
    'model/cutoff_convergence/run_cutoff_data.m'
    'model/scaled_errors/run_scaled_data.m'
    'model/penalty_condition/run_penalty_data.m'
    };
cfg.hCases = [ ...
    struct('Nc', 30, 't', 0, 'refines', 1:6), ...
    struct('Nc', 30, 't', 1, 'refines', 1:6)];
end
