# 第85轮计算与认证

本轮严格结论在 [core85.md](core85.md)，全部渐近量词由证明承担。本机证书检验熵常数、统一浓缩余量及匹配构造的实际组合费用，不把有限计算当作渐近证明。没有生成或认证新的有限 $F_{N,L(N)}$ 全输入，$L(N)$ 的存在用先固定维数、再空间极限的合法对角选择；不能宣称一个未计算的具体有限 $L$。

## 精确部分

1. 用160位 mpmath 根给候选端点，随后只用有理指数区间证明 $\mu-1-e^{-\mu}$ 的异号括号，宽度 $10^{-60}$。根函数严格递增，故括号为真正认证，不依赖浮点求根准确性。指数区间复用第84轮81项正级数/几何尾/逐次平方算法；Wolfram以独立实现重放符号。
2. 精确认证 $h_2((\mu-1)/\mu)=\mu-1+\log\mu>0.52$、指数隙 $\log\mu>0.24$。另以 $\epsilon=1/100$、$\delta=1/20$、$b=1/100$ 验证统一尾率
   \[
   b\delta-\frac{2b^2}{(1-b)^3}>1/4000,\quad
   2\epsilon^2=1/5000,\quad
   \log\mu-\delta-\mu\epsilon>9/50.
   \]
   这些宽裕常数独立支持目标预算率下的指数抑制，不给出未知渐近项的有限维起始位置。
3. 在12个维数 $N=2^{10},\ldots,2^{21}$，对校准 $c_N=\mu-\sqrt{\mu\log N/N}/2$、$r_N=e^{-c_N}$、$p_N=r_N/(1+r_N)$ 构造向外区间，令 $k_N$ 为 $Np_N+2\sqrt{N\log N}$ 的区间上端向上取整。全部 $k_N<N/2$。平方根用整数 `isqrt` 向外，指数用有理区间。
4. 精确计算 $\binom N{k_N}$ 并保存位长/高128位/SHA256；最大为1606849 bit。由相邻项比率可直接认证
   \[
   \binom N{k_N}\le V_N:=\sum_{j\le k_N}\binom Nj
   \le \binom N{k_N}\frac{N-k_N+1}{N-2k_N+1}.
   \]
   $V_N$ 正是Hamming源分配的实际模块倍数，未改用更强的振幅线性费用。对数区间复用 `sliver79.py`；无需存储百万位整数的十进制展开。

## Wolfram 与精度复核

Wolfram15运行12内核，独立求根、80/160位计算并复算全部精确二项整数的位长、高位及对数区间。9个符号门槛覆盖：增量MGF、均值及导数、熵式、临界熵恒等式、统一二阶指数矩、模式概率式、二项相邻比和有理尾率。

Python配置12进程，生成实际10进程接到工作；最终重放实际进程数见验收JSON。两软件保存65位计算结果并比较精度加倍稳定性；所有相关结果落在有理定向区间内。对数二项式的定向区间很窄；$\log V_N$ 的上下界还包括上述有限几何因子的间距，不能把后者也称为 $10^{-22}$ 宽的精确值。

所有计算确定性，无随机种子、优化求解器或尾部数值截断。Hoeffding直接给丢失平均 $\le2N^{-8}<N^{-2}$。本轮没有启动MATLAB；Wolfram并行核在结束后全部关闭。首个工作目录验收只因依赖哈希路径不在临时目录而停止，数值门槛已过；正式验收从最终项目代码目录完成，并登记实际依赖路径。

## 复现

从项目根目录执行，`OUT` 在仓库外：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r085-replay
"$PYTHON" code/verify85.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/verify85.wl "$OUT/python" "$OUT/wolfram"
"$PYTHON" code/check85.py "$OUT" "$OUT/validation85.json" --project-root .
```

[validation85.json](validation85.json) 登记最终源码、输入和既有证明SHA256。单向算子的精确常数来自既有r010/r028/r083接口，本轮没有重新对所有单角点常数做无必要计算。
