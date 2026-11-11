# 第93轮计算、误差与资源契约

严格范围见 [core93.md](core93.md)。有限计算验证独立的Gamma充分前提及尾多项式；不能据此给一般 $H_N$ 上界、加权定理起始维或精确命中分位数值。

## 21个精确Gamma前提

$N=2,3,4,8,16,32,64$，$p=9/10,19/20,49/50$。以 $b=\log(2/p-1)$，只用 $r=0$ 的Gamma增量，
\[
 P(\Gamma(N,1)\le Nb)\le(Nb)^N/N!<1-p
\]
足以严格证明 $s_p>Nb$。不需要求取最优侧比向量。Python的 `sliver79.logi` 给 $b$ 的有理向外区间，把上端代入正幂作精确Fraction比较；保存严格正余量的40位十进制有理包络、位长、SHA256和对数区间，避免把几千位分式直接输出。

Wolfram独立使用更粗但完全有理的 $b<2(1-p)/p$，仍对21例证明正余量；该上界与Python的不同，验收分别验证，未要求它们数值相等。实际Gamma CDF、$b$ 和幂上界的对数另用MP80/160位独立核对，仅作数值诊断。

## 140个尾多项式比值

$P_{N-1}(T)=\sum_{j=0}^{N-1}T^j/j!$，$a=0,1/4,1,4$，$T=0,a/2,a,a+1,a+4$。在98个 $T\ge a$ 案例，归一化真实尾比就是有理 $P(T-a)/P(T)\le1$。其余42个 $T<a$ 案例是 $e^{T-a}/P(T)\le1/P(T)$；后者是有理上界，不能当作等式。Wolfram保存的数值是实际指数比值，验收按实际公式独立计算后比较。

35个 $a=0$ 边界的实际比值等于一，其他105个实际比值严格小于一。其中21个 $T=0<a$ 案例的粗有理上界仍等于一，不能将上界与实际指数比值混淆。初版Wolfram把 $P(0)$ 的常数项写成 `0^0`，被完整验收发现；已显式写为 `1+Sum[...j>=1]`，全部重跑通过，旧输出不计入证据。

## 复放与数值精度

[experiment93.py](../../code/experiment93.py) 实际12个spawn进程；[verify93.wl](../../code/verify93.wl) 使用12个Wolfram15内核，每项256MiB、120秒，结束关闭并行核。[check93.py](../../code/check93.py) 再启动12进程重放161案，去掉PID/耗时后精确比对；检查紧凑包络、原分式hash、每种独立上界和实际数值，记录代码、4份输入和历史证明SHA256。

数值目标80/160位，各自另加保护位，保存65位。双软件相对容差 $10^{-45}$，没有以单位尺度掩盖极小量；正式最大差约 $4.22\times10^{-65}$。精确保证来自有理界及数学不等式，不来自保存65位相同。确定性正式批次无随机种子和优化器。

## 额外12核探路的范围

[extra93](extra93/README.md) 保存3个实际完成批次的最小汇总与原始执行源码。Gaussian批次为768个标量窗口。两个共同随机样本批次合计23,324个唯一参数/种子；每例2,000,000对抽样。每条增量的模型是 $Z_r=E+cB$，只覆盖共同动作、两个权重层的独立子族，不覆盖所有 $r_i,w_i$。

两经验CDF的DKW误差半宽为 $\sqrt{\log(4\cdot10^8)/(2S)}+2/S$，额外项支付线性分位插值的秩余量；两个边缘分布的并界每job至多 $10^{-8}$，同时覆盖其全部4个尾概率。整个有限网格并界为 $2.3324\times10^{-4}$。这不是严格的通用量化导数证书。种子由参数repr的SHA256前8字节生成，两个批次未重复已算680例；归一化差商不是无穷小导数。

原始JSONL不进仓库。原始本机源码保留当时路径；复制到 `/tmp` 同名位置可按相同软件重跑，再运行归档的流式汇总。额外脚本仅承诺本机复现，其他路径需调整并重新登记hash。最小摘要保留所有输入hash和最有区分力的10个参数/种子。后续第94轮探路不混入本轮结果。

## 内存和正式复现

所有批次由 [math_guard.py](../../code/math_guard.py) 保护，约1GiB系统报告可用内存余量；Python正式批次RSS预算3GiB、Wolfram6.5GiB。监控是自有进程树和macOS可用百分比的采样估计，不是OS硬预留。正式资源见 [Python](python/resources93.json)、[Wolfram](wolfram/resources93.json)。自有内核已退出，本轮未启动MATLAB。

在项目根目录，输出放仓库外：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r093-replay
"$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 3 --log "$OUT/python/run.log" --report "$OUT/python/resources93.json" -- "$PYTHON" code/experiment93.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel "$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log "$OUT/wolfram/run.log" --report "$OUT/wolfram/resources93.json" -- /usr/local/bin/wolframscript -file code/verify93.wl "$OUT/wolfram" 12 256
"$PYTHON" code/check93.py "$OUT" "$OUT/validation93.json" --project-root .
```

[Python证书](python/certificates93.json)、[Wolfram摘要](wolfram/verification93.json)、[最终验收](validation93.json)。
