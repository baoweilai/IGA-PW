function run_example_4(workflowName)
%Run one selected Example 4 task.

assert(exist('workflowName', 'var') == 1, 'run_example_4 requires a workflow name.');

exampleDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(exampleDir));
add_project_paths_local(projectDir);
add_example_paths(exampleDir);
outputSnapshot = snapshot_example_outputs(exampleDir);
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(exampleDir, 'data'));

switch string(workflowName)
    case "reference"
        add_workflow_paths(fullfile(exampleDir, 'model', 'reference'), ...
            {'config', 'core', 'operators', 'solver'});
        run_reference_solution(true, struct());

    case "h_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'h_convergence'), ...
            {'config', 'core', 'operators', 'solver'});
        run_h_data([1 2]);
        plot_h_convergence();

    case "cutoff_convergence"
        add_workflow_paths(fullfile(exampleDir, 'model', 'cutoff_convergence'), ...
            {'config', 'core', 'operators', 'solver'});
        run_cutoff_data(true, struct());
        plot_cutoff_convergence(struct());

    case "energy"
        add_workflow_paths(fullfile(exampleDir, 'model', 'energy'), ...
            {'config', 'core', 'operators', 'solver'});
        run_energy_data(true, struct());
        plot_energy(struct());

    case "helium_fields"
        add_workflow_paths(fullfile(exampleDir, 'model', 'helium_fields'), ...
            {'config', 'core', 'operators', 'solver'});
        plot_helium_fields();

    case "vout_fourier"
        add_workflow_paths(fullfile(exampleDir, 'model', 'vout_fourier'), ...
            {'config', 'core', 'operators', 'solver'});
        run_vout_2048_simple();

    otherwise
        error('Unknown Example 4 workflow: %s', workflowName);
end

sync_paper_outputs(exampleDir, 'example_4', workflowName, outputSnapshot);
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
