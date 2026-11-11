# 第92轮计算、误差与资源契约

严格物理归约和模型障碍分别由 [core92.md](core92.md) 与 [独立审查](audit92.md) 给出。数值只核对正态门槛秩位移、正部矩和尺度，不认证新的端点数值、未知历史常数或有限物理起始维。

## 参数与独立算法

$N=2^k$，$k=16,32,64,128,256,512,1024,2048,4096$，$B=1/4,1,4,16$，共36案例。每例用80/160位目标，输出65位。$\mu$ 使用既有有理包络核验，并在各软件内独立求根。

先解 $\bar\Phi(z)=\sigma z/\sqrt N$，再解
\[
 \bar\Phi(z_B)=\tau+BL/N,
 \quad r_B=\sigma\sqrt N(z-z_B).
\]
直接形成近分位的差，再放大为秩位移；未从两个约 $\mu N$ 的秩相减。Python/mpmath用两点局部求根和稳定对数尾，额外保护 $\lceil\log_{10}N/2\rceil+80$ 位。Wolfram独立用单初值解析数值根与更多保护 $\lceil\log_{10}N\rceil+80$ 位。目标是核对 $r_B/B\to1$ 的量级，有限参数不被要求已经等于极限。

对 $D=0,1,4$、$a=z+DL/(\sigma\sqrt N)$，计算
\[
 N\{\mathbb E(G-s)_+/N-\tau\}=DL+\sigma\sqrt N[\phi(a)-a\bar\Phi(a)],
\]
以及原正部平方矩的无量纲尺度
\[
 \mu\{(a^2+1)\Phi(a)+a\phi(a)\}/L.
\]
这些是同方差、无额外位移的高斯诊断模型。平方矩并非真实逆 $\Phi_N$ 期望的完整二阶修正，$D$、$B$ 也不是物理定理或第21轮常数的拟合或认证。

## 显示精度、残差与验收

大维数时，$z-z_B$ 和门槛尾概率的差远小于65位输出精度，因此展示的 $z,z_B$ 或不同 $B$ 的尾概率可相同。严格信号在保护位内计算的 $r_B>0$，按 $B$ 严格增加；不得从两份展示值重建这个微小差。

[check92.py](../../code/check92.py) 实际12进程重放全部Python案例，核对两精度和独立Wolfram；相对容差 $10^{-45}$，不能以 `max(1,abs(x))` 掩盖极小量误差。分别核验求根残差、参数域、概率与位移随 $B$ 的方向、同精度输出稳定性和资源读数，并登记源码、输入、历史依赖与核心/审核文档SHA256。65位输出相同只是该表示下的稳定性，非无限精度声明。

Wolfram另核对12项符号恒等式：正态密度、尾函数、stop-loss和正部矩导数，基底/超额恒等式、近一模块导数、逆二次系数和秩导数。所有有限计算均确定性，随机种子及优化求解状态不适用。$B_P$、$C_H$ 与起始维依赖解析证明；本轮无新增定向有限端点认证。

## 资源与复现

[experiment92.py](../../code/experiment92.py) 实际12进程，[verify92.wl](../../code/verify92.wl) 实际12个Wolfram15内核，单任务256MiB、120秒，结束 `CloseKernels[]`。总监控预留约1GiB报告可用内存，Python/Wolfram批次RSS预算3/6.5GiB。采样只监控自有进程树，是估计保护，不是OS硬预留；本轮未启动MATLAB。

在项目根目录复现，输出放仓库外：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r092-replay
"$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 3 --log "$OUT/python/run.log" --report "$OUT/python/resources92.json" -- "$PYTHON" code/experiment92.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel "$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log "$OUT/wolfram/run.log" --report "$OUT/wolfram/resources92.json" -- /usr/local/bin/wolframscript -file code/verify92.wl "$OUT/wolfram" 12 256
"$PYTHON" code/check92.py "$OUT" "$OUT/validation92.json" --project-root .
```

最小原始摘要：[Python](python/experiments92.json)、[Wolfram](wolfram/verification92.json)；资源：[Python](python/resources92.json)、[Wolfram](wolfram/resources92.json)；最终 [完整重放](validation92.json)。数学归约、模型内最优和一般真实常数仍分开登记。
