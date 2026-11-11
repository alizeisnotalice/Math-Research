# 第89轮计算、误差与证据契约

严格定理与数值范围分开登记：[core89.md](core89.md) 的所有物理输入渐近结论由证明给出。有限计算核对高斯比较、强 $L^p$ 根和显式费用公式，**没有认证物理定理的有限起始维或新的 $C_N$ 数值**。

## 算法、精度与参数

$N=2^k$，$k\in\{16,32,64,128,256,512,1024\}$。35个斜坡案例取 $K/\sigma=1,1.25,1.5,1.75,1.9$；28个分位数案例取 $\epsilon=1/N,1/\sqrt{NL},\sigma\rho/2,3\sigma\rho/4$；49个超额源案例取 $P=1/4,1,4$、临界方案及二次率 $c=1/10,1/5,3/10$。确定性计算，随机种子与优化求解状态不适用。

[experiment89.py](../../code/experiment89.py) 用本机Python3.12/mpmath1.3的12个进程，在80/160位目标下增加随 $k$ 变化的保护位。斜坡根直接以标准化位移求解，避免从两个约 $\mu N$ 的量相减；有保护Newton与导数残差检查。分位数按对数尾求根；$L^p$ 方程以 $v=p(b_p-\mu)$ 和稳定 `expm1`/正幂级数归一化计算。所有输出保留65位。

[verify89.wl](../../code/verify89.wl) 用Wolfram15实际启动12内核；斜坡根用高精度二分，$L^p$ 根直接计算未展开原方程及解析导数，与Python的重写不同。额外保护位为 $\lceil\log_{10}N\rceil+80$，足以抵消大项消去；另有12项符号门槛。每个任务限额256MiB、120秒，结束调用 `CloseKernels[]`。

## 极小数、边界与验收

- 大整数维数只是标量公式参数。分位数和斜坡残差按实际小尾尺度检查；软件差采用相对误差，没有以 `max(1,abs(x))` 隐去极小概率差。
- 若 $\log s<-10000$，`s` 字段保存 `null`，同时保留 $\log s$ 及明确标志。计算没有把该字段当作零阈值。
- Wolfram仅在 $\log s<-[2L+(wp+50)\log10]$ 时省略 $\log(1-s)$；保留被省略的 $N\log(1-s)$ 绝对误差上界的对数 $\log2+L+\log s$，验收要求其小于 $-(wp+40)\log10$。这是有界近似。
- 临界 $K$、$P$ 和 $c$ 边界按原公式计算；有限维超额费用差可为正，斜坡标准化位移也可未进入渐近符号，不将这些误报为失败或有限物理结论。
- [check89.py](../../code/check89.py) 重放全部112个Python案例，核对双精度、独立软件、残差和边界，并记录源码、输入、历史证明SHA256。接受相对容差 $10^{-45}$，实际最大差小于 $4.85\times10^{-65}$；输出65位的精度稳定不声称无限精度一致。
- 利用既有精确有理区间工具，另包络 $\mu$ 与 $2\sqrt\mu$。这是常数的算术认证，不是有限维端点认证。

## 资源和复现

[math_guard.py](../../code/math_guard.py) 监控自有进程树与macOS报告可用内存，低于1GiB或超出本批次RSS预算即停止本批次。该读数为估计，不能保证每一瞬间的OS硬预留。正式Python/Wolfram批次分别约2.05/13.43秒；峰值约372.27/1832.09MiB，可用内存最低分别约8683.52/7864.32MiB；均正常退出，无停止标志。本轮未启动MATLAB。

在项目根目录复现，输出放仓库外：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r089-replay
"$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 3 --log "$OUT/python/run.log" --report "$OUT/python/resources89.json" -- "$PYTHON" code/experiment89.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel "$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log "$OUT/wolfram/run.log" --report "$OUT/wolfram/resources89.json" -- /usr/local/bin/wolframscript -file code/verify89.wl "$OUT/wolfram" 12 256
"$PYTHON" code/check89.py "$OUT" "$OUT/validation89.json" --project-root .
```

最终原始摘要：[Python](python/experiments89.json)、[Wolfram](wolfram/verification89.json)、[Python资源](python/resources89.json)、[Wolfram资源](wolfram/resources89.json)。验收器再次实际启动12进程重算，不把输入自报的 `pass` 当作全部证据。
