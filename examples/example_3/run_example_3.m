function run_example_3(workflowName)
%Run one selected Example 3 task.

assert(exist('workflowName', 'var') == 1, 'run_example_3 requires a workflow name.');

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
            {'nurbs', 'iga', 'assembly', 'operators', 'core', 'solver'});
        run_h_data();
        plot_h_convergence();
        run_reference();
        plot_dg_error();

    case "cutoff_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'cutoff_convergence'), ...
            {'nurbs', 'iga', 'assembly', 'operators', 'core', 'solver'});
        run_cutoff_data();
        run_reference();
        plot_cutoff_dg();

    case "preconditioner"
        add_workflow_paths(fullfile(exampleDir, 'model', 'preconditioner'), ...
            {'nurbs', 'iga', 'assembly', 'operators', 'core', 'solver'});
        plot_preconditioner(struct());
        plot_preconditioner_residuals();

    case "method_comparison"
        plot_method_fields(fullfile(exampleDir, 'data', 'method_comparison', ...
            'method_fields'));

    otherwise
        error('Unknown Example 3 workflow: %s', workflowName);
end

sync_paper_outputs(exampleDir, 'example_3', workflowName, outputSnapshot);
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
