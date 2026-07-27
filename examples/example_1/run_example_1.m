function run_example_1(workflowName)
% Run one selected Example 1 task.

assert(exist('workflowName', 'var') == 1, 'run_example_1 requires a workflow name.');

% Set paths and enter the Example 1 data directory.
exampleDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(exampleDir));
cfg = config_example_1();
add_project_paths_local(projectDir);
add_example_paths(exampleDir);
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(exampleDir, 'data'));

% Run the selected data and figure workflow.
switch string(workflowName)
    case "h_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'h_convergence'), ...
            {'nurbs', 'dg', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core'});
        run_h_data(cfg.hCases);
        plot_h_convergence();
        plot_eigen();

    case "cutoff_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'cutoff_convergence'), ...
            {'nurbs', 'dg', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core'});
        run_cutoff_data();
        plot_cutoff_convergence();

    case "scaled_errors"
        add_workflow_paths(fullfile(exampleDir, 'model', 'scaled_errors'), ...
            {'nurbs', 'dg', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core'});
        run_scaled_data();
        plot_scaled_errors();

    case "penalty_condition"
        add_workflow_paths(fullfile(exampleDir, 'model', 'penalty_condition'), ...
            {'nurbs', 'iga', 'assembly', 'operators', 'core'});
        run_penalty_data(struct());
        plot_penalty(struct());

    case "method_comparison"
        add_workflow_paths(fullfile(exampleDir, 'model', 'method_comparison'), ...
            {'nurbs', 'iga', 'assembly', 'operators', 'core'});
        plot_method_comparison(struct());

    otherwise
        error('Unknown Example 1 workflow: %s', workflowName);
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
