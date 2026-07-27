function cfg = config_example_3()
% Return fixed Example 3 parameters.

cfg.example = 3;
cfg.name = 'Example 3';
cfg.outputPrefix = 'ex3';
cfg.figureScripts = {
    'plot_h_convergence'
    'plot_cutoff_dg'
    'plot_dg_error'
    'plot_preconditioner'
    'plot_preconditioner_residuals'
    'plot_method_fields'
    };
cfg.dataScripts = {
    'model/h_convergence/run_h_data.m'
    'model/h_convergence/run_reference.m'
    'model/cutoff_convergence/run_cutoff_data.m'
    'model/cutoff_convergence/run_reference.m'
    'model/method_comparison/run_method_data.m'
    };
end
