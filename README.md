# IGA-PW-DG Paper Codes

This repository contains the MATLAB code and retained numerical data for the revised manuscript *A Hybrid Discontinuous Galerkin Method with Isogeometric-Planewaves Coupling for Periodic Full-Potential Electronic Structure Calculations*.

## Folder Structure

- `src/`: shared assembly, discretization, solver, error, plotting, and utility code.
- `external/primme/`: third-party PRIMME MATLAB interface.
- `examples/example_1/` to `examples/example_4/`: the four numerical examples. Each example contains a runner, configuration, model code, figure or table code, and retained data.
- `examples/vout_table_test/`: the 2-D and 3-D FFT-Chebyshev performance tests used in Table 4.2.
- `tools/`: workflow and reference-data utilities.

All project paths are resolved from the repository location. Numerical-example outputs remain inside the corresponding example's `data/` folder, while the Table 4.2 CSV files are written inside `examples/vout_table_test/`.

## Running the Examples

There is no root-level all-example runner. Start MATLAB in the repository and run one workflow at a time:

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
run_example_4('hartree_comparison')
run_example_4('hartree_comparison_postprocess')
run_example_4('helium_fields')

cd ../vout_table_test
run_vout_table_2d()
run_vout_table_3d()
```

Run the Example 4 reference workflow before the other Example 4 workflows. The `hartree_comparison_postprocess` workflow reads retained comparison results without solving the eight cases again.

## Workflow Steps

1. The example runner selects a named workflow.
2. The path utilities add the shared source code and workflow-specific dependencies.
3. The configuration defines the physical model, discretization, solver, and retained cases.
4. The model code assembles the operators and solves or loads each case.
5. The figure or table code reads the saved data and creates the manuscript output.
6. Multi-panel manuscript figures are assembled from the separate PDF panels produced by the mapped scripts.

If a required file or folder is missing, MATLAB stops with an error.

## Paper Output Map

The numbering below follows the revised manuscript.

### Method and Implementation Results

| Paper item | Content | Workflow | Corresponding code |
|---|---|---|---|
| Table 4.1 | First-eigenvalue sensitivity to Chebyshev degree | From `examples/example_1`: `addpath('model'); run_cheb_sensitivity()` | `examples/example_1/model/run_cheb_sensitivity.m` |
| Table 4.2 | 2-D and 3-D masked FFT versus FFT-Chebyshev performance | `run_vout_table_2d()` and `run_vout_table_3d()` | `examples/vout_table_test/run_vout_table_2d.m` and `examples/vout_table_test/run_vout_table_3d.m` |
| Figure 4.2 | Condition number, total solution time, and residual history for the three preconditioners | `run_example_3('preconditioner')` | `examples/example_3/figures/plot_preconditioner.m` and `examples/example_3/figures/plot_preconditioner_residuals.m` |

### Numerical Examples

| Paper item | Content | Workflow | Data and output code |
|---|---|---|---|
| Figure 5.1 | Example 1 DOF and time comparison for PW, IGA, and IGA-PW | `run_example_1('method_comparison')` | Retained data: `examples/example_1/data/method_comparison/summary_method_compare.mat`; plot: `examples/example_1/figures/plot_method_comparison.m` |
| Figure 5.2 | Example 1 eigenvalue convergence in h, K = 30 | `run_example_1('h_convergence')` | Data: `examples/example_1/model/h_convergence/run_h_data.m`; plot: `examples/example_1/figures/plot_h_convergence.m` |
| Figure 5.3 | Example 1 eigenfunction convergence in DG norm, K = 30 | `run_example_1('h_convergence')` | Data: `examples/example_1/model/h_convergence/run_h_data.m`; error and plot: `examples/example_1/figures/plot_eigen.m` |
| Figure 5.5 | Example 1 eigenvalue and eigenfunction convergence in K | `run_example_1('cutoff_convergence')` | Data: `examples/example_1/model/cutoff_convergence/run_cutoff_data.m`; plot: `examples/example_1/figures/plot_cutoff_convergence.m` |
| Figure 5.6 | Example 1 scaled L2 and DG errors for alpha = 1, 2 | `run_example_1('scaled_errors')` | Data: `examples/example_1/model/scaled_errors/run_scaled_data.m`; plot: `examples/example_1/figures/plot_scaled_errors.m` |
| Figure 5.7 | Example 1 penalty sensitivity and conditioning | `run_example_1('penalty_condition')` | Error data and panel: `examples/example_1/model/penalty_condition/run_penalty_data.m`; condition-number panel: `examples/example_1/figures/plot_penalty.m` |
| Figure 5.8 | Example 2 eigenvalue and eigenfunction convergence in h | `run_example_2('h_convergence')` | Data: `examples/example_2/model/h_convergence/run_h_data.m` and `examples/example_2/model/h_convergence/run_reference.m`; plots: `examples/example_2/figures/plot_h_convergence.m` and `examples/example_2/figures/plot_eigen_error.m` |
| Figure 5.9 | Example 2 eigenvalue and eigenfunction convergence in K | `run_example_2('cutoff_convergence')` | Data: `examples/example_2/model/cutoff_convergence/run_cutoff_data.m`; plot: `examples/example_2/figures/plot_cutoff_convergence.m` |
| Figure 5.10 | Example 2 first eigenfunction and absolute errors for PW, IGA, and IGA-PW | `run_example_2('method_comparison')` | Solve, evaluate, and plot: `examples/example_2/figures/plot_method_fields.m` |
| Figure 5.11 | Example 3 h- and K-convergence of the first eigenpair | `run_example_3('h_convergence')` and `run_example_3('cutoff_convergence')` | Data: `examples/example_3/model/h_convergence/run_h_data.m`, `examples/example_3/model/h_convergence/run_reference.m`, and `examples/example_3/model/cutoff_convergence/run_cutoff_data.m`; plots: `examples/example_3/figures/plot_h_convergence.m`, `examples/example_3/figures/plot_dg_error.m`, and `examples/example_3/figures/plot_cutoff_dg.m` |
| Figure 5.12 | Example 3 absolute eigenfunction errors for PW, IGA, and IGA-PW | `run_example_3('method_comparison')` | Data: `examples/example_3/model/method_comparison/run_method_data.m`; plot: `examples/example_3/figures/plot_method_fields.m` |
| Figure 5.13 | Example 4 h- and K-convergence of the first eigenpair | `run_example_4('reference')`, `run_example_4('h_convergence')`, and `run_example_4('cutoff_convergence')` | Data: `examples/example_4/model/reference/run_reference_solution.m`, `examples/example_4/model/h_convergence/run_h_data.m`, and `examples/example_4/model/cutoff_convergence/run_cutoff_data.m`; plots: `examples/example_4/figures/plot_h_convergence.m` and `examples/example_4/figures/plot_cutoff_convergence.m` |
| Figure 5.14 | Example 4 helium orbital, density, and Hartree potential | `run_example_4('reference')`, `run_example_4('h_convergence')`, and `run_example_4('helium_fields')` | Data: `examples/example_4/model/reference/run_reference_solution.m` and `examples/example_4/model/h_convergence/run_h_data.m`; plot: `examples/example_4/figures/plot_helium_fields.m` |
| Figure 5.15 | Example 4 shifted periodic-energy convergence | `run_example_4('reference')` and `run_example_4('energy')` | Data: `examples/example_4/model/energy/run_energy_data.m`; plot: `examples/example_4/figures/plot_energy.m` |
| Table 5.1 | Eight-case IGA-PW and DG-PW comparison | `run_example_4('hartree_comparison')` | Data: `examples/example_4/model/hartree_comparison/run_hartree_comparison.m`; table: `examples/example_4/tables/build_comparison_table.m` |

## External Code

PRIMME is treated as third-party code and lives in `external/primme/`. Before running the Example 4 IGA-PW/DG-PW comparison, build its platform-specific RFP MEX solver as described in `examples/example_4/model/hartree_comparison/README.md`.
