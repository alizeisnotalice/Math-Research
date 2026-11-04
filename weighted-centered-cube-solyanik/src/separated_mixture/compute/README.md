# Round 20：外层混合权数值核查

## 口径

本轮研究的是外层模型
\[
H_\alpha(E)=\{x:\sup_{r>0}|E\cap[x-r,x+r]|/(2r)>\alpha\},
\qquad w(H_\alpha(E))/w(E).
\]
因此 halo 仅由 Lebesgue 长度决定，混合权只在最后测量同一个 \(H_\alpha(E)\) 和 \(E\)。前一版把 \(w\) 放入平均分母，实际是内禀模型；该错误试跑、脚本和缓存都保存在 `rejected/`，不作为本报告任何结论的证据。

## 已完成的实际计算

- `outer_mixture_batch.py --workers 15 --audit-full-pool`：15 个配置、5 个固定阈值，共 75 个外层 halo/加权比值检查；全部完成、0 失败。使用 Manager barrier，记录了 15 个不同的实际 worker PID：54347–54361（15个连续编号）。每配置还有两个奇点的局部缩放，因安全条件 \(5\epsilon<|z_1-z_2|\) 实际保留 208 行，而不是把近碰撞样本强行延伸到会碰到另一个奇点的尺度。
- 15 个配置覆盖 \(\tau\in\{.15,.25,.5,.75,.9\}\)、常规分离和 4 组近碰撞位置、不同 \(c_0,c_j\) 与三种有孔有限区间输入。每行用分别实现的 Fraction 端点分段算法和浮点端点-cell 算法计算全域 halo；75/75 的端点误差小于 \(2\times10^{-10}\)。这两套程序共享一维分段线性原理，故只是实现回归核对。
- 对每个公共 halo，按分量计算 \(w_i(H)/w_i(E)\)。混合比值是以 \(w_i(E)\) 为权的凸组合；75/75 数值满足它不超过分量公共-halo比值的最大值。这是代数恒等式的数值回归检查，不是新的极值定理。
- 对 \(m=2,3,4,8\) 记录了临界映射 \(\lambda=\rho^m\)、\(\Phi_m=1-m(1-\rho)\)，其中 \(\rho=(m-1)/m\pm10^{-3}\)。上侧 \(\Phi_m>0\) 可定义 halo；下侧 \(\Phi_m<0\)，不构造虚假的负阈值 halo。
- `matlab_outer_independent.m` 实际经连接的 MATLAB R2026a Update 4 运行：276 个局部有孔输入（24 组混合权），全部有限且凸组合检查通过；最大归一化局部差为 3.18116。脚本随后仅删除了 `clearvars`，避免影响用户 MATLAB 会话；没有因此重复运行。Code Analyzer 只有性能/格式信息提示，没有 error 或 warning 级数学/执行问题。
- `nd_blowup_check.py`：在 \(n=2,3\)、五个 \(\tau<n\)、七个尺度及半径 \(\epsilon/2,\epsilon\) 的两类固定有理方体上，共保存 140 个局部质量检查。中心 \(\ell_\infty\) 奇性质量使用精确公式 \(c2^n n r^\tau/\tau\)；远离奇点的背景和另一个奇性分量由 12/24/48 阶 tensor Gauss 交叉比较。
- `nd_analytic_verify.py` 不重跑 Gauss 节点，而对已经保存的 140 个 Gauss-24 值用闭式残差复核。因为另一个奇点在 \((1.3,0,\ldots,0)\) 且 \(r\le1/8\)，其积分为
\[
c(2r)^{n-1}\frac{(D+r)^{\sigma-n+1}-(D-r)^{\sigma-n+1}}{\sigma-n+1}.
\]
闭式和保存的 Gauss-24 混合质量最大绝对差为 \(1.421\times10^{-14}\)。这验证求积实现，不把阶次吻合称为严格包围。

## 文件

- `outer_mixture_batch.py`：可移植的 15 进程外层主程序，带 checkpoint、PID 与参数版本。
- `outer_mixture_summary.json`：完整可再生结果；正文不复制其全部记录。
- `matlab_outer_independent.m` 与 `matlab_outer_summary.json`：MATLAB 独立交叉核对。
- `nd_blowup_check.py`、`nd_blowup_summary.json`、`nd_analytic_verify.py`、`nd_analytic_audit.json`：二维和三维局部质量检查及闭式审计。
- `audit_runs/correct_outer_11worker/`：正确口径但未实际达到 15 活跃 PID 的预检查，保留用于可追溯性；最终证据采用根目录的 audit-full-pool 结果。

## 限制

这些有限输入、固定阈值和局部缩放只支持候选“有限互异奇点混合权的局部 blow-up 受单分量控制”的数值图景。它们不证明全体集合的精确 Tauberian 常数、全域上界、极限存在，或任何高维 halo 结论。

## 归档复现

依次运行 `python3 outer_mixture_batch.py --workers 15 --audit-full-pool`、`python3 nd_blowup_check.py`、`python3 nd_analytic_verify.py`；MATLAB运行`matlab_outer_independent.m`。主外层几何仅需标准库，ND程序需要NumPy。`run_summary.json`保留最终原运行摘要；上述全表和检查点均可再生，不入Git。父代理为归档版本补了缓存配置/参数版本检查、计算与复用计数、当前worker PID及失败退出码；新目录实际运行与缓存续跑均核验，续跑不计作新数值证据。

ND闭式审计使用已有结果，不重复求积；140个box质量既不是140个独立权，也不是高维halo下界。MATLAB和Python几何共用同一端点分段原则，属于跨实现数值核对。完整误用模型的试跑仅在本机scratch/rejected，Git不保存其错误脚本或结果。
