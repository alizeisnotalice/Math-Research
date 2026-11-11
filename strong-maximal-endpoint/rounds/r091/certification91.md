# 第91轮计算、误差与资源契约

[core91.md](core91.md) 的全空间结论由严格证明与 [audit91.md](audit91.md) 支持。下列计算认证有限二维几何、价格区间，并独立核对高维标量尺度；不认证物理定理的有限起始维、未知常数或一般锐对数系数。

## 精确二维输入与分母

基集为 $[0,1]^2\cup([1,2]\times[0,1/2])$，内部不交。周期 $\varepsilon=1/m$，$m=8,16,32,64$，低相位份额 $p=1/4,1/3,1/2,3/4$。低高度为1、高高度为5/4。每个盒子的全局低份额精确等于 $p$，原分母为
\[
 J_2=\tfrac32\{p+(1-p)\tfrac54[1+\log(5/4)]\}.
\]
定向有理对数给出该分母和每个矩形实际局部模数的区间。不能将含对数的分母直接当作有理数。

同一接收点 $(3/2,101/200)$ 到基集的 $\ell^\infty$ 距离为 $1/200$。16个主案例各检查4个实际矩形：

- $[1,2]\times[0,51/100]$；
- $[1,2]\times[0,13/25]$；
- $[11/10,17/10]\times[1/7,13/25]$；
- $[7/6,19/11]\times[2/9,53/100]$。

前两个覆盖完整周期，均是严格成功见证；后两个检验非零边界误差，逐例记录成功状态，不要求它们全部成功。实际平均为 $5b/4-\ell/4$，其中 $b$、$\ell$ 分别是基层和低层的全矩形密度，严格余量与原模块各自计算。

另独立复现审查中的单位方块例：$p=1/3$、$m=6,12,24,48$、$R=[1/8,12/11]\times[1/7,5/7]$。四个密度偏差精确为 $11/48960,11/257040,0,0$，平均精确为 $41393/39168,217327/205632,539/510,539/510$，均严格超过一。

Python使用标准库Fraction，对相位条带的正部平方二重原函数作四角相减。Wolfram独立使用有理多边形与两条斜线逐次裁剪，再以鞋带公式求面积；没有读入Python结果作几何输入。所有68个矩形的几何量要求精确相等，误差均不超过解析 $2K_0\varepsilon/b_x$ 界。这个粗界只证明包含性，不把有限例的更小偏差误说成普适锐率。

## 高维标量、尾概率与价格下界

$N=2^k$，$k=64,96,128,192,256,384,512,768,1024,1536,2048,4096$。这些仅为标量公式参数，不生成 $N$ 维网格。使用80/160位目标，输出65位。

计算两种尾层：$t_0=\sqrt{L/2-4\log L}$ 与 $t_1=\sqrt{L-8\log L}$，比较
\[
 \log\frac{\bar\Phi(t_0)}{LN^{-1/4}},\qquad
 \log\frac{\bar\Phi(t_1+L^{-1/2})}{N^{-1/2}L^{5/2}}.
\]
第二个比例增长用于单尾扰动，不能直接套用第90轮较大的全CDF误差。全部用对数尾，避免极小概率的补数消去；检查概率随给定参数序列递减、根残差、精度加倍稳定性。

对 $\tau=\sigma z/\sqrt N$，顶端原子给出的**价格比下界**是 $\tau\Phi_N(1+2\tau)/\Phi_N(1+\tau)$，不是整个期望价格。为避免从两个巨大值直接相减，Python先算
\[
 g_1=\log\Phi_N(1+\tau)-d,\quad g_2=\log\Phi_N(1+2\tau)-2d,
\]
再算 $\log\tau+g_2-g_1$，即该下界对数减去 $d$。Python用稳定 `log1p`，额外保护位为 $\lceil\log_{10}N/2\rceil+80$；Wolfram直接用原始对数，额外保护位为 $\lceil\log_{10}N\rceil+80$，支付极小量加一及其后巨大费用消去的精度损失。未生成巨大精确指数幂。

## 验收、资源与复现

[experiment91.py](../../code/experiment91.py) 实际启动12个本机Python进程；[verify91.wl](../../code/verify91.wl) 实际使用12个Wolfram15内核，另有12项符号门槛。每个Wolfram案例限制256MiB和120秒，最后关闭并行核。

[check91.py](../../code/check91.py) 用新的12进程池重放全部参数，忽略进程编号后精确比较Python结果；几何与Wolfram逐分式相等，标量以 $10^{-45}$ 相对容差比较两软件和两种精度。原分母的65位输出比较时允许其序列化舍入误差，严格价格包络来自定向对数。脚本、4份输入JSON、历史证明与本轮核心/审查均记录SHA256。确定性计算，随机种子与优化求解器不适用。

预留约1GiB系统报告可用内存，使用 [math_guard.py](../../code/math_guard.py) 监控自有进程树RSS，Python预算3GiB、Wolfram6.5GiB。采样读数是估计，不是OS硬预留；正式资源见 [Python](python/resources91.json) 与 [Wolfram](wolfram/resources91.json)。早期Python批次末尾漏记部分进程编号而未生成完成JSON，已修正并全量重跑，不计入通过证据。Wolfram首版的矩形参数导出少加列表括号，验收器发现坐标列表被截为单值；已修复导出并完整重跑，几何算法本身未变。旧文件不作正式证据。本轮未启动MATLAB。

在项目根目录执行，输出放仓库外：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r091-replay
"$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 3 --log "$OUT/python/run.log" --report "$OUT/python/resources91.json" -- "$PYTHON" code/experiment91.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel "$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log "$OUT/wolfram/run.log" --report "$OUT/wolfram/resources91.json" -- /usr/local/bin/wolframscript -file code/verify91.wl "$OUT/wolfram" 12 256
"$PYTHON" code/check91.py "$OUT" "$OUT/validation91.json" --project-root .
```

[Python最小摘要](python/experiments91.json)、[Wolfram最小摘要](wolfram/verification91.json)、[最终重放](validation91.json)。有限核验不替代同一例外集、所有矩形和任意原输入的解析量词。
