# 第86轮核验与复现契约

全维结论见 [free_rectangles86.md](free_rectangles86.md) 的解析证明。计算只核验任意端点统一矩界、异常边费用与选参，不把渐近起效维数假设为已认证。第85轮的真实源和对角空间截断继续作为已证明依赖，不重复生成全量输入。

## 定向常数

复用第85轮经异号认证的 $\mu$ 区间。Python有理运算对 $a=10^{-2},10^{-3},10^{-4},10^{-6}$ 认证
\[
g_a=\min\{\mu-1,\mu-m(r_*-a),\mu-m(r_*+a)\}>0,
\quad D_a=\log\frac{\mu+a}{\mu-a}>0.
\]
对 $a=10^{-6}$，另认证 $2D_a<1/40000$。结合第85轮的严格浓缩余量，可以按证明的选参顺序取足够小固定 $a$，再令 $N\to\infty$。源成本 $j_N$ 的有限值仍须保留；未求出 $O_a(\sqrt{N\log N})$ 异常坐标界的数值起效点。JSON中的示例指数是充分大维数的解析选择，不能作为任意有限 $N$ 的数值证书。

所有有理区间由最终验收器完整重放。区间对数与正指数复用已有代码，不重复归档旧依赖。没有人为消耗内存的操作。

## 任意边端点测试

统一写 $U\sim\mathrm{Unif}[1,d]$（端点无序），$W=(-\log|U|)_+$。18种端点覆盖：短向内、跨原点、临界最优、跨过 $-1$ 的长边、接近零或一、向外长边。每种测试 $b=1/100,1/25,1/10,1/4$，共72例。

- Python12进程配置；实际4进程分得短工作。mpmath80/160位给闭式均值及MGF，SciPy用变换后的独立积分复核，最大双精度绝对误差小于 $3\cdot10^{-15}$。
- Wolfram15实际12内核，对同一72例以80/160位独立数值积分，输入另加40位保护，保存45位结果。分别与闭式、Python及精度加倍结果比较。
- 8个符号门槛核验两类统一尾界、短/长内向均值、由尾界推出二阶指数矩、$3/8$概率资格代数、Hardy $L^2$ 分部积分恒等式、cap2平方费用。全端点统一证明在正文，有限采样不能替代它。

第一遍积分虽通过数值门槛，但 $d\approx1$ 的相减损失触发输入精度警告；为输入增加40位保护后重跑全组，最终仅归档修正后的结果。随机种子不适用；没有优化求解器，没有数值尾截断或重要性抽样。

独立审查使用本机MATLAB作公式sanity check，实际成功命令保存在 [audit86.m](../../code/audit86.m)。该会话未查询版本，故不补造版本信息；PID11609退出后经进程检查消失。本轮MATLAB不参与认证。所有Wolfram并行核已关闭，脚本进程退出。

## 复现

从项目根目录运行，输出放仓库外：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r086-replay
"$PYTHON" code/verify86.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/verify86.wl "$OUT/python" "$OUT/wolfram"
"$PYTHON" code/check86.py "$OUT" "$OUT/validation86.json" --project-root .
```

可选MATLAB命令：`matlab -batch "run('code/audit86.m')"`，批处理结束即退出。

[validation86.json](validation86.json) 保存最终代码、输入和历史证明的SHA256及验收门槛。认证不包含新的 $C_N$ 上下界、一般残差分解或一般固定物理角点列预算。
