# 第27轮：完整小齐次次数弱型渐近

主稿paper/small_homogeneity_weak.tex的SW.1–SW.9是对所有输入的解析结论，不由有限原子实验替代。固定范数K、正有界零次齐次轮廓a，hτ=ρ^(τ−n)a，a_*为角向本质下确界，Ea为锥概率均值，则τW_{K,hτ}→nEa/a_*。匹配上界由局部/中尺度/大尺度见证给出；匹配下界在低角集合选非零Lebesgue点并平滑原子。对contrast≤R的角向族和所有n维范数还一致有相对误差O_{n,R}(τ^(1/(n+1)))，未声称该余项阶尖锐。

## 计算的独立角色

`three_scale_guards.wl`：精确核验中尺度径向积分、包含球标量优化和极限，失败非零退出。原空间核为ρ^(τ−2n)，乘径向Jacobian后才是r^(τ−n−1)dr。

`cube_shell_batch.py`：120位Decimal、15进程、36个唯一(n,b,τ)，K=[−1,1]^n，n=1,2,3，a=1于x1>0、a=b于x1<0，b=1,2,5，τ=.1,.01,10^−4,10^−6。原子在(1,0,…,0)，质量2^n(1+τ)^n。光晕的加权壳体积

A(s)=[min(s,r+1)+b min(s,r−1)] [2min(s,r)]^(n−1)，r=1+τ。

程序按r−1,r,r+1分段解析积分∫s^(τ−n) A'(s)ds，独立计算输入成本，并比较给定半径的包含球下界及三尺度上界。签名包括源字节hash、参数、解释器、平台；不匹配缓存移入stale。36项fresh与36项replay已实际完成，最小统计和样例在run_summary.json；全表/日志/检查点仅在本机保留，可复算，不重复提交。

`cube_shell_independent.m`：MATLAB R2026a独立对数变量积分，36行；与Python的τR最大绝对差3.553e−15。最初直接端点积分的warning输出已废弃。`cube_shell_guards.wl`在Wolfram15计算9个(n,b)的极限和36行正性。三个程序确认显式输入公式，均不是完整弱型上界的计算证书。

## 单参数唯一性路线的严格反例

同一双原子族Rτ(δ)在全区间0<δ<1/3上未必只有一个驻点。`critical_fraction_certificate.py`以整数50次根和Fraction区间，对τ=11/50认证G=(1+δ)N'−N在δ=17/100,18/100,8/25,33/100的符号为+,-,-,+。必要有理证书在同名JSON；连续性给两个互异开区间内至少各有一根，未声称恰有两根。该反例防止把局部边界层唯一性误升级为全局约化。

子族的局部高阶展开未进入正文或本归档；本轮以完整算子匹配首项为主，不增加缺乏全域价值的高阶表达。准备的多原子探索器没有执行，未虚计探索评估数。

## 复现

在本目录执行：

```
python3 cube_shell_batch.py
python3 critical_fraction_certificate.py
wolframscript -file three_scale_guards.wl
wolframscript -file cube_shell_guards.wl
```

MATLAB运行cube_shell_independent.m。Python只需标准库；程序默认15个工作进程，首波屏障验证实际并发。所有有限输出均与严格证明分别标注。
