# Periodic helium comparison

This Example 4 workflow compares IGA-PW and DG-PW for eight fixed cases.

| Method | Cases `(K,p,r)` or `(K,n_r,L_m)` |
|---|---|
| IGA-PW | `(10,2,2)`, `(10,2,4)`, `(15,2,2)`, `(15,2,4)` |
| DG-PW | `(10,2,2)`, `(10,3,3)`, `(15,2,2)`, `(15,3,3)` |

Both methods use the Example 4 reference at `K=30`, `p=2`, and
`nelem=32`. Exchange-correlation is disabled.

## Steps

1. `default_config` defines the eight cases, numerical settings, and
   repository-relative paths.
2. `solve_iga` computes the four IGA-PW cases.
3. `solve_dg` and `dg_scf` compute the four DG-PW cases.
4. `run_hartree_comparison` stores each case and creates the combined result.
5. `build_comparison_table` and `plot_comparison` create the table and plots.

## Run

From `examples/example_4`, run the eight cases and postprocess them:

```matlab
run_example_4('hartree_comparison')
```

Postprocess the retained results without recomputing:

```matlab
run_example_4('hartree_comparison_postprocess')
```

The numerical results are stored under
`data/result/hartree_comparison`. Summary files use short names:

- `iga.csv` and `iga.mat`;
- `dg.csv` and `dg.mat`;
- `comparison.csv` and `comparison.mat`;
- `table.csv`;
- `reference.csv` and `reference.mat`;
- `run.log`.

Individual cases use `iga_1.mat` through `iga_4.mat` and `dg_1.mat`
through `dg_4.mat`. DG-PW SCF histories and stage times use the same case
indices.

The IGA-PW path uses PRIMME and the packed/RFP preconditioner. Build the
platform-specific MEX file from the workflow directory with:

```matlab
cd solver
build_rfp_mex
```
