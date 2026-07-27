function run_example_2(workflowName)
% Run one selected Example 2 task.

assert(exist('workflowName', 'var') == 1, 'run_example_2 requires a workflow name.');

% Set paths and enter the Example 2 data directory.
exampleDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(exampleDir));
cfg = config_example_2();
add_project_paths_local(projectDir);
add_example_paths(exampleDir);
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(exampleDir, 'data'));

% Run the selected data and figure workflow.
switch string(workflowName)
    case "h_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'h_convergence'), ...
            {'nurbs', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core'});
        run_h_data(cfg.hRefines);
        run_reference(cfg.reference.Nc, cfg.reference.pdeg, ...
            cfg.reference.refines);
        plot_h_convergence();
        plot_eigen_error();

    case "cutoff_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'cutoff_convergence'), ...
            {'nurbs', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core'});
        run_cutoff_data();
        plot_cutoff_convergence();

    case "method_comparison"
        add_workflow_paths(fullfile(exampleDir, 'model', 'method_comparison'), ...
            {'nurbs', 'iga', 'operators', 'core'});
        plot_method_fields();

    otherwise
        error('Unknown Example 2 workflow: %s', workflowName);
end

end

function add_project_paths_local(projectDir)
% Add source and external solver paths.

requiredDirs = {
    fullfile(projectDir, 'src', 'assembly')
    fullfile(projectDir, 'src', 'dg')
    fullfile(projectDir, 'src', 'iga')
    fullfile(projectDir, 'src', 'pw')
    fullfile(projectDir, 'src', 'nurbs')
    fullfile(projectDir, 'src', 'error_norms')
    fullfile(projectDir, 'src', 'postprocess')
    fullfile(projectDir, 'src', 'plotting')
    fullfile(projectDir, 'src', 'utils')
    fullfile(projectDir, 'external', 'primme', 'Matlab')
    };

for k = 1:numel(requiredDirs)
    assert(isfolder(requiredDirs{k}), 'Missing project directory: %s', requiredDirs{k});
    addpath(requiredDirs{k});
end
end
