function cfg = config_example_2()
% Return fixed Example 2 parameters.

cfg.example = 2;
cfg.name = 'Example 2';
cfg.outputPrefix = 'ex2';
cfg.figureScripts = {
    'plot_h_convergence'
    'plot_cutoff_convergence'
    'plot_eigen_error'
    'plot_method_fields'
    };
cfg.dataScripts = {
    'model/h_convergence/run_h_data.m'
    'model/cutoff_convergence/run_cutoff_data.m'
    };
cfg.hRefines = 2:6;
cfg.reference = struct('Nc', 45, 'pdeg', 3, 'refines', 8);
end
