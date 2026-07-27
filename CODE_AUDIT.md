# IGA-PW-DG 完整代码核查报告

核查日期：2026-07-27

## 1. 核查基准与范围

- 核查目录：仓库根目录
- 上游仓库：<https://github.com/baoweilai/IGA-PW>
- 上游树：`origin/main`，提交 `03a955be72ff8705a8cf2c7d6b5a4943d4ff7cc1`
- 本地 HEAD：`97295ace220c1e194b1b6192ec2701e1c85d6434`
- 代码范围：所有 `.m`、`.ps1`、`.c` 和 `.cpp` 文件，包括被 Git 忽略但实际存在的运行脚本。
- 数据、图片、PDF、日志和 MAT 文件不计入代码数量。
- `external/primme` 是第三方代码，单独核查，不改写其内部实现。

本地与上游历史已经分叉，因此“多出来”和“缺少”均按两个提交的文件树直接比较，不按提交先后推断。

| 项目 | 数量 |
|---|---:|
| 当前全部代码文件 | 782 |
| 当前第一方代码文件 | 777 |
| 当前第一方 MATLAB 文件 | 774 |
| 上游代码文件 | 796 |
| 当前相对上游新增路径 | 41 |
| 当前相对上游缺少路径 | 55 |

## 2. 总体结论

| 核查项 | 结论 | 证据 |
|---|---|---|
| 本机绝对路径依赖 | 通过 | 第一方代码中的 Windows 绝对路径为 0 |
| 指定的防御性结构 | 通过 | 第一方代码中对应的五类关键字均为 0 |
| 后备分支或后备标记 | 通过 | 第一方代码中相关标记为 0，缺失输入直接停止 |
| 提醒类注释 | 通过 | TODO、FIXME、注意、提醒、务必、later、must、should、avoid 等提醒式注释为 0 |
| 图片输出格式 | 通过 | 第一方代码只生成 PDF，仓库内非 PDF 图片文件为 0 |
| MATLAB 静态分析 | 通过但有性能建议 | 774 个文件仅剩 136 条 `SPRIX`，均为稀疏索引可能较慢 |
| 不可达或未使用代码 | 通过 | `UNRCH`、`DEFNU`、`INUSD`、`NASGU`、`PREALL`、`AGROW` 均为 0 |
| 函数名与文件名 | 通过 | 781 个函数文件无主函数名不匹配 |
| 文件名长度 | 部分通过 | 超过 36 个字符的文件名为 0；超过 30 个字符的技术性文件名仍有 42 个 |
| 代码简洁性 | 部分通过 | 无效代码已清理，但仍有 99 组完全重复文件，共涉及 634 个文件 |
| 完整数值复算 | 未执行 | 已做轻量数值和编译验证；未运行全部正式高成本参数扫掠 |

第三方 `external/primme` 中仍存在本报告所检查的防御性写法，主要位于 `make.m`、`primme_eigs.m` 和 `primme_svds.m`。这些文件属于外部依赖，保持上游内容不变。

## 3. 已完成的主要修改

上游的组织逻辑是共享实现放在 `src`，PRIMME 放在 `external`，四个算例分别从 `run_example_i.m` 进入，再依次完成路径初始化、案例配置、数据计算和绘图。当前目录保留了这条主调用链，没有增加根目录全量运行器，也没有把不同工作流同时加入 MATLAB 路径。新增的 Hartree 和验证流程都作为独立具名工作流接入。

1. 所有第一方 MATLAB 和 PowerShell 启动代码均通过文件位置推导项目根目录；MATLAB 可执行文件由系统命令发现。
2. Example 1 和 Example 2 的正式参数列表已放入配置并由入口显式传入。
3. Example 4 的公共算子目录已纳入工作流路径初始化。
4. Vout 性能测试集中到 `examples/vout_table_test`，只保留论文表格中的 4 个 2-D 案例和 2 个 3-D 案例。
5. 删除 9 份未被调用且内容完全相同的 `eigifp.m`。
6. Example 3 的 12 个方向性界面装配文件由 3 个统一快速装配文件替代。
7. 3 个长名称的非线性装配文件改为 `assemble_nonlinear_nurbs.m`。
8. 历史运行脚本改为短文件名，并删除其中的本机路径。
9. 删除未使用变量、无效预分配、失效的分析器抑制、未调用局部函数和重复计算。
10. 注释只说明当前代码执行的计算、装配、读取、绘图或验证操作。
11. 删除 Example 1 的 Smooth 专属工作流 28 个代码文件和 54 个结果文件，并清理入口、配置和文档引用。
12. Example 4 Hartree 比较只保留 4 个 IGA-PW 和 4 个 DG-PW 案例，删除 15 个独立重复测试代码文件和 9 个重复测试结果文件。
13. 删除 Example 1 的门限测试和特征函数灵敏度测试及其结果；保留的 Chebyshev 灵敏度脚本改为 `run_cheb_sensitivity.m`，结果集中到 `cheb_sensitivity`。
14. `run_cheb_sensitivity.m` 只计算第一特征值，删除第 3/4 特征值相关计算；结果只保存为 `.mat`，结果目录中的日志、CSV 和状态记录均已删除；未启用的 `save_* = false`、空输出目录和强制重算开关也已移除。
15. 全目录删除固定为 `false` 的保存、绘图、强制重算、Smoke、交换关联、PW 磁盘缓存和 Legacy tensor 配置及其不可达分支；Example 4 的能量流程只保存所需的 `run.mat` 内容，六份 3-D SCF 求解器删除不可达的 PRIMME 历史输出和零次迭代回退路径。
16. 所有第一方绘图和 SCF 诊断代码统一只导出 PDF，删除 PNG 输出开关、路径、返回字段和旧格式同步列表，并删除仓库中 143 个 PNG 结果文件及 3 个相关历史日志。
17. Example 3 三方法对比统一为 `run_method_data.m` 和 `plot_method_fields.m`；正式参数为 PW `Nc=30`、IGA `p=1,nElem=64`、IGA-PW `Nc=10,p=1,nElem=40`，参考解为 `Nc=40,p=2,refine=7`；删除旧图并将数据和 PDF 放入唯一正式目录。
18. Vout 性能测试只保留 4 个 2-D 案例和 2 个 3-D 案例，每项计算 3 次并保存平均值；删除 24 个旧流程代码文件、2 份逐次结果 CSV、内存占比门槛和额外诊断输出。
19. 删除共享论文输出的快照、复制、目录创建和清理比较代码；四个 Example 的数据、表格和 PDF 均保留在各自的 `data` 目录。
20. 删除 Example 1 中未被代码读取的 `cheb_sensitivity/history` 旧参数归档；删除 Example 4 代码中重复保存 `scf_history` 的逻辑、Hartree 比较的四份 `_scf.csv` 及其生成代码；Example 3 的三份求解器迭代历史由残差图直接读取，因此保留。

## 4. 工作流分步说明

1. `run_example_i.m` 选择一个具名工作流。
2. 路径工具加入共享源码、当前工作流目录及其依赖目录。
3. 配置文件给出物理量、离散参数、求解器参数和案例列表。
4. 模型代码组装算子并求解指定案例。
5. 绘图或表格代码读取结果并生成论文输出。
6. 数值结果、表格和 PDF 保存在对应 `examples/example_i/data` 目录。

## 5. 相对上游多出来的全部代码

这里的“多出来”仅表示上游提交中没有同一路径，不等于代码无用。41 个路径分布为：Example 1 有 6 个、Example 2 有 1 个、Example 3 有 7 个、Example 4 有 25 个、Vout 性能测试有 2 个。

### Example 1：6 个

```text
examples/example_1/figures/cutoff_eigenfunction_errors.m
examples/example_1/figures/plot_eigen.m
examples/example_1/model/cutoff_convergence/assembly/assemble_DG_square_interface_fast.m
examples/example_1/model/h_convergence/assembly/assemble_DG_square_interface_fast.m
examples/example_1/model/run_cheb_sensitivity.m
examples/example_1/model/scaled_errors/assembly/assemble_DG_square_interface_fast.m
```

### Example 2：1 个

```text
examples/example_2/figures/cutoff_eigerr.m
```

### Example 3：7 个

```text
examples/example_3/model/method_comparison/run_method_data.m
examples/example_3/model/cutoff_convergence/assembly/assemble_DG_square_interface_fast.m
examples/example_3/model/cutoff_convergence/operators/assemble_nonlinear_nurbs.m
examples/example_3/model/h_convergence/assembly/assemble_DG_square_interface_fast.m
examples/example_3/model/h_convergence/operators/assemble_nonlinear_nurbs.m
examples/example_3/model/preconditioner/assembly/assemble_DG_square_interface_fast.m
examples/example_3/model/preconditioner/operators/assemble_nonlinear_nurbs.m
```

### Example 4：25 个

```text
examples/example_4/figures/plot_comparison.m
examples/example_4/model/common/operators/assemble_vext_direct_3D.m
examples/example_4/model/hartree_comparison/config/default_config.m
examples/example_4/model/hartree_comparison/core/dg_scf.m
examples/example_4/model/hartree_comparison/core/solve_dg.m
examples/example_4/model/hartree_comparison/core/solve_iga.m
examples/example_4/model/hartree_comparison/operators/build_ewald_potential.m
examples/example_4/model/hartree_comparison/operators/fftcheb3.m
examples/example_4/model/hartree_comparison/operators/fftcheb3_hartree.m
examples/example_4/model/hartree_comparison/operators/int_overlap.m
examples/example_4/model/hartree_comparison/operators/int_overlap_radius.m
examples/example_4/model/hartree_comparison/operators/int_overlap_radius_3.m
examples/example_4/model/hartree_comparison/operators/int_overlap_radius_d.m
examples/example_4/model/hartree_comparison/operators/r_basis.m
examples/example_4/model/hartree_comparison/operators/r_basis_d.m
examples/example_4/model/hartree_comparison/operators/spherical_bessel.m
examples/example_4/model/hartree_comparison/operators/spherical_harmonic.m
examples/example_4/model/hartree_comparison/operators/spherical_harmonic_xyz.m
examples/example_4/model/hartree_comparison/operators/Wigner3j.m
examples/example_4/model/hartree_comparison/run_hartree_comparison.m
examples/example_4/model/hartree_comparison/solver/build_rfp_mex.m
examples/example_4/model/hartree_comparison/solver/build_trace_rfp.m
examples/example_4/model/hartree_comparison/solver/rfp_chol_mex.c
examples/example_4/model/hartree_comparison/solver/tbprec_rfp.m
examples/example_4/tables/build_comparison_table.m
```

### Vout 性能测试：2 个

```text
examples/vout_table_test/run_vout_table_2d.m
examples/vout_table_test/run_vout_table_3d.m
```

## 6. 相对上游缺少路径的主要处理

这些路径均有明确处理结果，不是遗漏复制。

- `examples/example_1/figures/plot_eigen_error.m` 由 `plot_eigen.m` 替代。
- 9 份 `eigifp.m` 是未调用的完全重复文件，已经删除：

```text
src/solvers/eigifp.m
examples/example_1/model/cutoff_convergence/solver/eigifp.m
examples/example_1/model/h_convergence/solver/eigifp.m
examples/example_1/model/scaled_errors/solver/eigifp.m
examples/example_2/model/cutoff_convergence/solver/eigifp.m
examples/example_2/model/h_convergence/solver/eigifp.m
examples/example_3/model/cutoff_convergence/solver/eigifp.m
examples/example_3/model/h_convergence/solver/eigifp.m
examples/example_3/model/preconditioner/solver/eigifp.m
```

- Example 3 的 12 个方向性边界装配文件由各工作流的 `assemble_DG_square_interface_fast.m` 替代：

```text
examples/example_3/model/cutoff_convergence/assembly/IGA_DG_Bottom_Edge_Assemble.m
examples/example_3/model/cutoff_convergence/assembly/IGA_DG_Left_Edge_Assemble.m
examples/example_3/model/cutoff_convergence/assembly/IGA_DG_Right_Edge_Assemble.m
examples/example_3/model/cutoff_convergence/assembly/IGA_DG_Top_Edge_Assemble.m
examples/example_3/model/h_convergence/assembly/IGA_DG_Bottom_Edge_Assemble.m
examples/example_3/model/h_convergence/assembly/IGA_DG_Left_Edge_Assemble.m
examples/example_3/model/h_convergence/assembly/IGA_DG_Right_Edge_Assemble.m
examples/example_3/model/h_convergence/assembly/IGA_DG_Top_Edge_Assemble.m
examples/example_3/model/preconditioner/assembly/IGA_DG_Bottom_Edge_Assemble.m
examples/example_3/model/preconditioner/assembly/IGA_DG_Left_Edge_Assemble.m
examples/example_3/model/preconditioner/assembly/IGA_DG_Right_Edge_Assemble.m
examples/example_3/model/preconditioner/assembly/IGA_DG_Top_Edge_Assemble.m
```

- 3 个 `assemble_nonlinear_nurbs_from_samples.m` 改名为 `assemble_nonlinear_nurbs.m`。
- `paper_output_dirs.m`、`snapshot_example_outputs.m`、`sync_paper_outputs.m`、`audit_codebase.m`、`clean_outputs.m` 和 `compare_results.m` 随共享输出流程一起删除。

## 7. 保留的多参数测试和案例代码

结论：仓库其他工作流仍有多参数测试。Example 4 Hartree 比较已按指定范围收敛为固定 8 个案例。

| 文件或模式 | 保留的案例 |
|---|---|
| `run_cheb_sensitivity.m` | 只计算第一特征值；`p=2`；`(refine,K)=(4,25),(6,25),(7,35)` 与 `n=10,50,100` 的 9 个组合；FFT 网格 500；只写出 `.mat` |
| `run_vout_table_2d.m` | 20000² masked FFT，加上 `(grid,n)=(300²,50),(600²,50),(300²,100)`；每项计算 3 次并保存平均值 |
| `run_vout_table_3d.m` | 1000³ masked FFT 与 300³、`n=100` 的 FFT-Chebyshev；每项计算 3 次并保存平均值 |
| Hartree `run` | IGA-PW：`(10,2,2)`、`(10,2,4)`、`(15,2,2)`、`(15,2,4)`；DG-PW：`(10,2,2)`、`(10,3,3)`、`(15,2,2)`、`(15,3,3)` |

Vout 结果只保存两份平均值 CSV；逐次运行 CSV、2048³/`n=80` 旧流程、内存占比门槛和额外诊断报告均已删除。

## 8. 简洁性与命名核查

- 已删除所有静态分析识别出的未使用赋值、无效预分配、未使用输入和不可达分支。
- 19 个 MATLAB 文件仍超过 500 行，主要是 3-D SCF、复杂绘图和比较求解器。
- 99 组完全重复文件说明工作流本地仍复制了大量 NURBS 和 IGA 基础函数。
- 这些重复文件暂未集中到 `src`，因为当前上游逻辑依靠工作流局部路径隔离；一次性合并会扩大数值回归范围。
- 42 个超过 30 字符的文件名均不超过 36 字符，主要描述装配对象、势函数插值或 Chebyshev 修正。
- 所有函数文件的主函数名均与文件名一致。

超过 30 字符的 42 个文件对应 9 个不同的文件基名：

| 文件基名 | 长度 | 重复数 |
|---|---:|---:|
| `build_inner_chebyshev_correction_3D` | 35 | 5 |
| `assemble_NURBS_potential_interp_3D` | 34 | 5 |
| `generate_example1_scaled_reference` | 34 | 1 |
| `build_affine_cube_fast_support_3D` | 33 | 5 |
| `assemble_DG_square_interface_fast` | 33 | 8 |
| `update_pw_potential_from_grid_3D` | 32 | 5 |
| `assemble_nonlinear_pw_from_grid` | 31 | 3 |
| `IGA_DG_Left_Right_Edge_Assemble` | 31 | 5 |
| `IGA_DG_Top_Bottom_Edge_Assemble` | 31 | 5 |

因此，当前代码已经去掉明确的冗余和历史残留，但不能认定为“全局最简”。进一步压缩需要把工作流局部公共函数统一到 `src`，并对四个 Example 做完整数值回归。

## 9. 验证结果

| 验证 | 结果 |
|---|---|
| 第一方关键字、后备标记、本机路径和提醒注释扫描 | 全部 0 命中 |
| 固定关闭配置项扫描 | `opts/cfg/common/runOpts.* = false` 为 0 命中 |
| 非 PDF 图片输出和文件 | 代码引用为 0，仓库文件为 0 |
| Example 3 三方法对比 | 参数与三个 `run.mat` 一致；旧代码路径和旧图目录为 0 |
| MATLAB R2025a `checkcode` | 787 个文件；仅 136 条 `SPRIX` |
| Example 4 Hartree 案例集合 | 4 个 IGA-PW 和 4 个 DG-PW 参数、DOF、误差与时间全部匹配 |
| 18 个工作流目录及 Example 4 公共算子路径 | 通过 |
| 3 个 PowerShell 文件语法解析 | 通过 |
| 函数名与文件名一致性 | 781 个函数文件全部通过 |
| `git diff --check` | 通过 |
| 球谐加法定理 | 最大误差 `1.5543122344752192e-15` |
| Wigner 3-j 已知值 | 误差 `1.1102230246251565e-16` |
| Wigner 3-j 选择规则 | 通过 |
| NURBS 曲线两级细化 | 通过 |
| RFP C/MEX 编译 | MinGW64 编译通过 |
| RFP 2×2 线性系统 | 残差 `0` |

`SPRIX` 表示稀疏矩阵按索引反复更新可能较慢，不是语法错误或未使用代码。其位置集中在 48 个装配文件。改变这些语句需要重写稀疏装配策略，因此本轮保持上游数值逻辑。

本轮没有执行全部正式 SCF、20000² FFT、1000³ FFT 和所有收敛扫掠。现有验证证明入口、路径、语法、轻量数学恒等式和 MEX 核心可用，但不能替代完整论文结果复算。
