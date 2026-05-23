function run_example_1(workflowName)
%Run one selected Example 1 task.

assert(exist('workflowName', 'var') == 1, 'run_example_1 requires a workflow name.');

exampleDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(exampleDir));
add_project_paths_local(projectDir);
add_example_paths(exampleDir);
outputSnapshot = snapshot_example_outputs(exampleDir);
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(exampleDir, 'data'));

switch string(workflowName)
    case "h_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'h_convergence'), ...
            {'nurbs', 'dg', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core', 'solver'});
        run_h_data();
        plot_h_convergence();
        plot_eigen_error();

    case "cutoff_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'cutoff_convergence'), ...
            {'nurbs', 'dg', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core', 'solver'});
        run_cutoff_data();
        plot_cutoff_convergence();

    case "scaled_errors"
        add_workflow_paths(fullfile(exampleDir, 'model', 'scaled_errors'), ...
            {'nurbs', 'dg', 'iga', 'assembly', ...
            'operators', 'error_norms', 'core', 'solver'});
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

sync_paper_outputs(exampleDir, 'example_1', workflowName, outputSnapshot);
end

function add_project_paths_local(projectDir)
%Add source and external solver paths.

requiredDirs = {
    fullfile(projectDir, 'src', 'assembly')
    fullfile(projectDir, 'src', 'dg')
    fullfile(projectDir, 'src', 'iga')
    fullfile(projectDir, 'src', 'pw')
    fullfile(projectDir, 'src', 'nurbs')
    fullfile(projectDir, 'src', 'solvers')
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
