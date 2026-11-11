# 第88轮计算、误差与资源契约

严格范围见 [core88.md](core88.md)：同一个仿射积源、整个接收patch上的反射中心预算上下匹配；另给任意输入的充分残差公式。数值未给所有输入选择源片/proposal，未登记新的 $C_N$ 上下界。

## 精确证书

- $N=2^{10},\ldots,2^{21}$，两种深度slab的预算下界保留精确 $(N/2)!$、$\pi<22/7$ 和真实模块上界 $J(F)\le3(21/20)^N$。既有定向对数/正指数函数只维护一份。
- 整个patch的显式dyadic grid在第 $j$ 层有精确整数节点数 $(\lceil128/\sqrt{\delta_j}\rceil+1)^N$，$\delta_j=2^{-j}/4$。有限项至 $j=3\log_2N+20$，用有理区间计算真实 $\log\chi$；其余项用严格几何尾控制。全无限构造费用的对数因此有向外区间，并非只计算峰值。
- 对 $x\le-300$，用 $e^x<10^{-100}$ 包络logistic及log-sum-exp，依据 $e^3>10$；大正参数同理。没有在极小尾上使用浮点减法或生成巨大指数整数。
- 外扩相容中心的**单位矩形**体积，在 $N=2,4,\ldots,256$、$\delta=N^{-2},1/10,1/2$、$r=0,1/4,1$ 共72例，用81项正级数
  \[
  I_N(s)=\sum_{k\ge0}\frac{s^{N+k}}{(N-1)!\,k!\,(N+k)},\quad s=\log(1+\delta),
  \]
  和明确几何余项认证；相对级数尾小于 $10^{-100}$。大有理数的对数先以分子/分母的192个高位构造严格包络，避免对数级数直接放大万位分式。

## 独立软件与数值边界

Python12进程与Wolfram15的12内核独立80/160位核对。Wolfram另加80位输入保护，检查12项符号门槛、12维dyadic费用和72个中心体积。Wolfram数值logistic在参数绝对值超过400时使用有明确误差界的近似，log-sum-exp中相对差小于−400的数值项同理；合计相对遗漏不超过 $170e^{-400}$。这些只是独立数值核对，精确包含区间由Python负责。

保存65位结果的一致性只说明该输出精度下稳定，不是无限精度误差声明。中心积分先缩放到 $[0,1]$，再加回尺度对数，避免体积下溢。随机种子不适用，全部确定性；没有优化求解器。

## 内存保护与已处理的失败

按用户最新指令使用多线程并预留约1GiB系统可用内存。每个Wolfram子任务设置额外分配上限256MiB、90秒限时；[math_guard.py](../../code/math_guard.py) 监控自有进程树RSS和macOS报告的可用内存，触线只停止本批次。它是监控保护，不是OS级硬预留。

最终12核批次的监测峰值为2341.34MiB，系统报告的可用内存最小11468.8MiB，未触发停止条件；见 [resources88.json](wolfram/resources88.json)。全部并行核及本批次脚本均已退出。本轮没有启动MATLAB。

首个Python原型尝试把巨大有理级数余项输出成字符串，触发4300位限制；改为保存已证明的紧凑上界并优化大分式对数后，完整重跑通过。Wolfram原型混合精确 `N Log[M]` 与近似量，存在将指数化简成巨型整数幂的风险；Luna只读审查确认了该风险。最终脚本在构造每个对数项时就数值化，并设置内存保护。早期无完成文件的批次不计入通过证据。

## 复现

在项目根目录，输出放仓库外；默认12核，保留总内存监控与单项限额：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r088-replay
"$PYTHON" code/experiment88.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel "$PYTHON" code/math_guard.py --reserve-gib 1 --max-rss-gib 6.5 --log "$OUT/wolfram/run.log" --report "$OUT/wolfram/resources88.json" -- /usr/local/bin/wolframscript -file code/verify88.wl "$OUT/wolfram" 12 256
"$PYTHON" code/check88.py "$OUT" "$OUT/validation88.json" --project-root .
```

内存紧张时将Wolfram命令最后的12改为较小核数；保持单项上限及预留监控。最终验收 [validation88.json](validation88.json) 完整重放精确输出，保存代码、输入、既有证明与资源报告SHA256。
