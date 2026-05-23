# IGA-PW-DG Paper Codes

This directory is the refactored paper-reproducibility entry point for the IGA--PW--DG examples. The original folders are kept beside this directory until baseline data and refactored data are compared.

## Folder Structure

- `src/`: shared MATLAB implementation code.
- `external/primme/`: third-party PRIMME MATLAB interface.
- `examples/example_1/` to `examples/example_4/`: paper examples with matching `run_example_i.m`, `config_example_i.m`, `figures/`, `tables/`, `data/`, and workflow-specific `model/` folders.
- `paper_outputs/figures/example_i/`: final manuscript figures grouped by example and workflow.
- `paper_outputs/tables/example_i/`: final manuscript tables grouped by example and workflow.
- `paper_outputs/data/example_i/`: manuscript data grouped by example and workflow.
- `tools/`: audit, comparison, and cleanup utilities.

## Running the Examples

There is no root-level setup script and no root-level all-example runner. Each example runner initializes the shared paths internally and runs one named workflow at a time:

```matlab
cd examples/example_1
run_example_1('h_convergence')
run_example_1('cutoff_convergence')
run_example_1('scaled_errors')
run_example_1('penalty_condition')
run_example_1('method_comparison')

cd ../example_2
run_example_2('h_convergence')
run_example_2('cutoff_convergence')
run_example_2('method_comparison')

cd ../example_3
run_example_3('h_convergence')
run_example_3('cutoff_convergence')
run_example_3('preconditioner')
run_example_3('method_comparison')

cd ../example_4
run_example_4('reference')
run_example_4('h_convergence')
run_example_4('cutoff_convergence')
run_example_4('energy')
run_example_4('helium_fields')
run_example_4('vout_fourier')
```

## Paper Output Map

| Paper output | Script |
|---|---|
| Example 1 h-convergence figure | `examples/example_1/figures/plot_h_convergence.m` |
| Example 1 cutoff-convergence figure | `examples/example_1/figures/plot_cutoff_convergence.m` |
| Example 1 eigenfunction-error figure | `examples/example_1/figures/plot_eigen_error.m` |
| Example 1 scaled-error figure | `examples/example_1/figures/plot_scaled_errors.m` |
| Example 1 penalty/condition figure | `examples/example_1/figures/plot_penalty.m` |
| Example 1 method-comparison figure | `examples/example_1/figures/plot_method_comparison.m` |
| Example 2 h-convergence figure | `examples/example_2/figures/plot_h_convergence.m` |
| Example 2 cutoff-convergence figure | `examples/example_2/figures/plot_cutoff_convergence.m` |
| Example 2 eigenfunction-error figure | `examples/example_2/figures/plot_eigen_error.m` |
| Example 2 method-comparison fields | `examples/example_2/figures/plot_method_fields.m` |
| Example 3 h-convergence figure | `examples/example_3/figures/plot_h_convergence.m` |
| Example 3 cutoff/DG-convergence figure | `examples/example_3/figures/plot_cutoff_dg.m` |
| Example 3 eigenfunction DG-error figure | `examples/example_3/figures/plot_dg_error.m` |
| Example 3 preconditioner comparison | `examples/example_3/figures/plot_preconditioner.m` |
| Example 3 preconditioner residual histories | `examples/example_3/figures/plot_preconditioner_residuals.m` |
| Example 3 method-error fields | `examples/example_3/figures/plot_method_fields.m` |
| Example 4 h-convergence figure | `examples/example_4/figures/plot_h_convergence.m` |
| Example 4 cutoff-convergence figure | `examples/example_4/figures/plot_cutoff_convergence.m` |
| Example 4 energy-convergence figure | `examples/example_4/figures/plot_energy.m` |
| Example 4 helium field figures | `examples/example_4/figures/plot_helium_fields.m` |
| Example 4 outer-potential Fourier summary | `examples/example_4/model/vout_fourier/run_vout_2048_simple.m` |

The current repository does not contain table-generation scripts under `examples/example_i/tables/`. Generated table artifacts are copied into `paper_outputs/tables/` when a workflow produces them.

## Baseline and Validation

- Before-refactoring audit: `paper_outputs/data/codebase_audit_before.txt`.
- Planned file-move map: `paper_outputs/data/planned_file_move_map_before.txt`.
- Comparison script: `tools/compare_results.m`.
- Comparison report target: `paper_outputs/data/baseline_vs_refactored_report.txt`.

Baseline MAT files should be stored in `paper_outputs/data/baseline/`. Refactored MAT files should be stored in `paper_outputs/data/refactored/` with matching names ending in `_refactored.mat`. `tools/compare_results.m` compares numerical arrays with relative tolerance `1e-12`.

## External Code

PRIMME is treated as third-party code and lives in `external/primme/`. Its internal `nargin`, printing, and defensive checks are not rewritten in this refactor.
