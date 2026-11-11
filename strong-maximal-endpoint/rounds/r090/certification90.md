# 第90轮计算、误差与资源契约

本轮物理定理由 [core90.md](core90.md) 的源变换、重排、纯化和接收体积证明给出；计算只认证标量二次价格余项，核对量级。没有有限维物理起始点、最佳误差常数或新端点值的数值认证。

## 精确算术

维数 $N=2^k$，$k=16,24,32,48,64,96,128,192,256,384,512,768,1024$。使用第85轮已认证的 $\mu$ 有理区间和现有平方根、定向对数函数。先包络 $\sqrt{\mu\log N}$ 再除以精确整数 $2^{k/2}$，避免将极小帽高交给固定绝对网格的平方根工具。

在 $x=2\sqrt{\mu L/N}$，三阶导数为
\[
 f_N^{(3)}(x)=\frac{2+(N-1)H(\log(1+x))}{(1+x)^3},\quad
 H(u)=\frac{2u^2+7u+7}{(1+u)^3},\quad
 H'(u)=-\frac{2u^2+10u+14}{(1+u)^4}<0.
\]
它为正且递减，积分Taylor公式给出 $x^3f_N^{(3)}(x)/6\le R_N(x)\le(7N-5)x^3/6$。将帽高上下端、对数上端代入相应单调因子可得有理上下界；对巨大分式复用第88轮192个高位的严格对数包络，保存紧凑对数区间。该区间是解析余项的证书，不是假定数值落在窄误差带。

有限CDF网格误差尺度 $h=LN^{-1/4}$ 的恒等式 $[L^2/(h\sqrt N)]/h=1$ 精确成立。网格节点数以有理上界 $\lceil10\sigma\sqrt L/h\rceil+5$ 控制，验收严格小于 $N$；未分配这些节点或高维数组。

## 独立浮点核验

[experiment90.py](../../code/experiment90.py) 使用本机Python/mpmath的12进程。正余项直接以
\[
 R_N(x)=\tfrac{x^3}{2}\int_0^1(1-t)^2 f_N^{(3)}(xt)\,dt
\]
求值，避免大项相减；80/160位目标另加 $\lceil\log_{10}N/2\rceil+40$ 位。独立 [verify90.wl](../../code/verify90.wl) 用原始对数价格直接减去二次项，另加 $\lceil\log_{10}N\rceil+80$ 位抵消消去。两软件还求正态分位根、$h$、$d/N$ 和相对尺度，记录65位输出。

Wolfram15实际启动12内核，每项限制256MiB、120秒，结束关闭并行核；12项符号门槛检查三阶公式、正性/递减性及尺度恒等式。初版有一项把 `FullSimplify` 输出与未化简多项式作结构比较而误报，改为将等式整体化简后全量重跑通过；不计失败批次。

[check90.py](../../code/check90.py) 重新启动12进程重放13案例，核对有理界、两种精度、独立软件相对差与求根残差，保存代码/输入/历史依赖SHA256。浮点容差 $10^{-45}$；数值只作交叉核验，完整验证结果见 [validation90.json](validation90.json)。确定性计算，随机种子和优化求解状态不适用。

## 资源和复现

预留约1GiB可用内存，Python批次RSS预算3GiB，Wolfram6.5GiB，且有每项限制。监控使用macOS报告的可用百分比和自有进程树采样，属于估计保护，不是OS硬预留。最终批次资源见 [Python资源](python/resources90.json) 与 [Wolfram资源](wolfram/resources90.json)；均正常退出、无停止条件，本轮未启动MATLAB。

在项目根目录执行：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r090-replay
"$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 3 --log "$OUT/python/run.log" --report "$OUT/python/resources90.json" -- "$PYTHON" code/experiment90.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel "$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log "$OUT/wolfram/run.log" --report "$OUT/wolfram/resources90.json" -- /usr/local/bin/wolframscript -file code/verify90.wl "$OUT/wolfram" 12 256
"$PYTHON" code/check90.py "$OUT" "$OUT/validation90.json" --project-root .
```

原始最小结果：[Python](python/experiments90.json)、[Wolfram](wolfram/verification90.json)。保存的65位相同不构成无限精度误差声明；解析余项上下界和物理定理的保证分别来自各自证明。
