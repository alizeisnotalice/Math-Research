# 第84轮计算与认证契约

本轮的全维定理来自 [obstruction84.md](obstruction84.md) 的解析证明；有限计算检查恒等式、保守常数和下界公式。所有 Orlicz 分母均为真实 $\Phi_N$ 分母。文件中的 `log_actual_modular_upper` 是 $\log(\int F+\int F^N)$，是 $\log J_N(F)$ 的上界，不能当作精确模块。

## 精确部分

- Python `Fraction` 对4个正指数使用81项 Taylor 正级数、几何尾和逐次平方，向外舍入至 $10^{-80}$；11个有理不等式门槛均通过。Mathematica用独立 Wolfram 实现逐项重放这些有理区间。
- 12个维数 $N=2^{10},\ldots,2^{21}$、两种 slab。精确生成 $(N/2)!$，最大整数为19458756 bit；只保存高位、位长、SHA256及定向对数区间。对数复用第79轮 `sliver79.py` 的有理正级数，依赖只维护一份。最终定向对数区间最大宽度为 $9.5420484\cdot10^{-23}$。
- 中央 Irwin–Hall 概率用整数交错和精确相消，再约分为有理数；$N=2,4,\ldots,1024$、$p=2,3$ 共20例。输出50位向外区间与精确分子/分母哈希。Wolfram独立整数计算20例全部通过；这是小维数复核，全维 slab 体积下界另有偶单峰卷积证明。
- Wolfram 15独立检查9个符号门槛：真实平均的二次损失、两种导数、最优均值参数、外点范围、源质量、N阶矩原函数、对数下界导数、$22/7-\pi$ 的正积分恒等式。

## 数值交叉部分

Python/mpmath与Wolfram均使用80、160位两级精度，保存65位结果。所有记录值落在相应定向区间内；保存值的精度加倍差为0，跨软件最大差为 $4.043\cdot10^{-58}$ 以下。这是保存位数下的一致性，不是对浮点误差的区间证明。浮点锐化值和真实模块上界另列，不把它们登记为认证新定理。

Wolfram配置12核，结果实际出现全部12个内核编号。Python配置12进程，证书生成实际8进程接到工作，slab生成实际4进程，最终重放实际6进程。短作业不保证各进程都分得任务。没有人为填充内存或空转负载；最大阶乘采用约2.4 MB整数，避免构造高维网格。

独立审查还短批调用本机 MATLAB R2026b，稳定对数版本复核原始粗分母公式的5个值及球体常数。已实际执行的稳定命令见 [audit84.m](../../code/audit84.m)。先前直接幂/指数诊断有上溢或下溢的字段未用于证据，MATLAB不参与任何精确证书。审查所启动 MATLAB PID10076 已执行 `exit` 并确认退出；Wolfram并行核通过 `CloseKernels[]` 退出，脚本进程正常结束。

首遍 Wolfram 核验中，对 $u/(1-u)^2=u/(-1+u)^2$ 使用结构相等比较，导致一个符号门槛假失败。修为假设 $0<u<1$ 下化简等式后，全部核验重跑通过。没有更改数学公式；失败产物不归档为完成证据。

## 可复现命令

在 `strong-maximal-endpoint/` 目录执行；`PYTHON` 可替换为带 mpmath 的 Python 3.12，`OUT` 为仓库外路径。

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r084-replay
"$PYTHON" code/verify84.py "$OUT/python"
"$PYTHON" code/slabs84.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/verify84.wl "$OUT/python" "$OUT/wolfram"
"$PYTHON" code/check84.py "$OUT" "$OUT/validation84.json" --project-root .
```

可选 MATLAB 稳定复核：`matlab -batch "run('code/audit84.m')"`，批处理结束即退出，不保留常驻引擎。

最终验收 [validation84.json](validation84.json) 保存代码、输入及历史依赖 SHA256，精确重放全部阶乘证书和20个概率。随机种子不适用：本轮全部确定性，无优化求解器、无随机抽样；没有概率尾截断，小量全用精确和或对数表达式。认证范围不含固定维数预算有限性，也不含新的 $C_N$ 上下界。
