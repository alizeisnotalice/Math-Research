# 第87轮核验契约

本轮的主要效率定理由 [corner87.md](corner87.md) 的解析证明及第22轮已验收的统一参数定位共同给出。它不是拟合有限维斜率。关键区别：第22轮给 $D\le C\log N$，但未提供本轮可直接代入的具体统一 $C$ 和数值起效维数；本轮也未补造它们。

## 严格几何与定向常数

取固定 $e=1/50$、$r_*=\mu-1$，认证 $[r_*-e,r_*+e]\subset[1/4,3/10]$、统一曲率下界 $133000/65910>2$、异常边均值缺口至少 $e^2$。直接有理区间给外框膨胀 $R_e$、cell比率 $P_e$ 和球体积计费常数 $C_{\rm geom}=R_e^2/[r_-(1+r_-)]$，向外压缩至 $10^{-60}$ 格点，不保留冗长的等价分式。

Wolfram15独立核验10个符号门槛：外框长度比、cell比、Jacobian/源盒体积、均值二阶导数、对数比较曲率、临界曲率、精确有理曲率、长跨零边的Exp随机支配、Gaussian球体积驻点及有限条件式的单调导数。180位根独立落在Python几何区间对应的正确范围。

## 有限条件计算的范围

若另已知全部成功动作满足 $D\le\gamma\log N$，正文(10)给出明确有限维 $K_E$ 上界。测试 $N=2^{12},\ldots,2^{23}$，$\gamma=1/100,1/10,1,10$，共48个输入。35个满足该简化公式的组合/单调性门槛，13个明确排除。排除输入不是认证值，也不能当作反例或通过结果。

对35例，Python完整有理对数区间与Wolfram80/160位分别结算同一个保守上界表达式；保存65位，最终验收检查每个结果都在定向区间内并完整重放所有输入。这个区间包围的是**条件上界公式的数值**，不是对 $K_E$ 的双侧区间，更不是实际 $C_N$ 的新认证。

Python配置12进程，实际分得任务的进程数见JSON；Wolfram实际12内核。全部确定性，无随机种子、求解器或概率尾部数值截断。本轮没有启动MATLAB，Wolfram并行核均在结束后关闭。

第一遍Wolfram脚本将双变量列表用于 `Table` 迭代器，导致语法求值中止，未生成完成证据；修成显式 `Map/Function` 后重新运行通过。常数压缩后再完整重跑，以最终代码/输入为准。未把未完成产物计入验收。

## 复现

从项目根目录执行：

```sh
PYTHON=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r087-replay
"$PYTHON" code/verify87.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/verify87.wl "$OUT/python" "$OUT/wolfram"
"$PYTHON" code/check87.py "$OUT" "$OUT/validation87.json" --project-root .
```

[validation87.json](validation87.json) 登记代码、输入和关键历史证明SHA256。复用已有有理对数/指数函数，只维护一份依赖。本轮未重新计算第85轮的近Q源，也未声称认证具体有限 $L(N)$ 或一般锐上界。
