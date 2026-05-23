function cfg = config_example_4()
%Return fixed Example 4 parameters.

cfg.example = 4;
cfg.name = 'Example 4';
cfg.outputPrefix = 'ex4';
cfg.figureScripts = {
    'plot_h_convergence'
    'plot_cutoff_convergence'
    'plot_energy'
    'plot_helium_fields'
    };
cfg.dataScripts = {
    'model/reference/run_reference_solution.m'
    'model/h_convergence/run_h_data.m'
    'model/cutoff_convergence/run_cutoff_data.m'
    'model/energy/run_energy_data.m'
    'model/vout_fourier/run_vout_2048_simple.m'
    };
end
