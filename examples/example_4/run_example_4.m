function run_example_4(workflowName)
% Run one selected Example 4 task.

assert(exist('workflowName', 'var') == 1, 'run_example_4 requires a workflow name.');

% Set paths and enter the Example 4 data directory.
exampleDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(exampleDir));
add_project_paths_local(projectDir);
add_example_paths(exampleDir);
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(exampleDir, 'data'));

% Run the selected data and figure workflow.
switch string(workflowName)
    case "reference"
        add_workflow_paths(fullfile(exampleDir, 'model', 'reference'), ...
            {'config', 'core', 'operators', 'solver'});
        run_reference_solution();

    case "h_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'h_convergence'), ...
            {'config', 'core', 'operators', 'solver'});
        run_h_data([1 2]);
        plot_h_convergence();

    case "cutoff_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'cutoff_convergence'), ...
            {'config', 'core', 'operators', 'solver'});
        run_cutoff_data();
        plot_cutoff_convergence(struct());

    case "energy"
        add_workflow_paths(fullfile(exampleDir, 'model', 'energy'), ...
            {'config', 'core', 'operators', 'solver'});
        run_energy_data();
        plot_energy(struct());

    case "hartree_comparison"
        add_workflow_paths(fullfile(exampleDir, 'model', 'hartree_comparison'), ...
            {'config', 'core', 'operators', 'solver'});
        run_hartree_comparison("run");
        build_comparison_table();
        plot_comparison();

    case "hartree_comparison_postprocess"
        add_workflow_paths(fullfile(exampleDir, 'model', 'hartree_comparison'), ...
            {'config', 'core', 'operators', 'solver'});
        run_hartree_comparison("check");
        build_comparison_table();
        plot_comparison();

    case "helium_fields"
        add_workflow_paths(fullfile(exampleDir, 'model', 'helium_fields'), ...
            {'config', 'core', 'operators', 'solver'});
        plot_helium_fields();

    otherwise
        error('Unknown Example 4 workflow: %s', workflowName);
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
