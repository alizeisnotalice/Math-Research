# R048 main.tex 独立严格审查

日期：2026-11-05。审查者：Sol，独占 main 审查目录。

**裁决：PROOF-VALID，所审 main 的数学陈述与传递无需修补。** 主定理的几何输入为已经发布且经两独审确认的 R046 全输入八阶上包、R041 每个充分小参数的实际有限八阶下界；新 `eighth.tex` 的整模块证明另由独审负责，本报告不重复审查，也不将尚未冻结的模块当成缺陷。最终完整组装引用及编译检查应在该模块冻结后补做。

## 1. 冻结对象和本次范围

完整读取：

`/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/scratch/r048-main-typeset/main.tex`

实际 SHA256 为：

`b976cf01d89b73ff2ec540a19aebc5efbc44848cd2d2fc6bf271667eb4eb63f5`。

与指定冻结值相同。并对照旧 repo main 的完整 diff、未变 global/sharp、R046 正式报告及审报的结论接口、R041 下界与代数审报的 J_n 接口。此前 main 的精确约化和七阶系数已经独审；本次仍核对更新是否改变其数学条件或传递强度。

沿用 math-paper-writing 的量词、结论强度及计算证据规则。未改源、未改仓库、未调用新代理、未重复编译，仅写本报告。固定维数 n≥2，所有常数可依赖 n；m≥1 不预先有界，可随 ε 变化。

## 2. P7 基线下正八阶结论与归一化

main 39–53 行摘要、117–178 行主定理及 446–468 行证明使用一致对象：

\[
 a=v^n=2^nt,\quad v=2t^{1/n},\quad
 P_{7,n}(t)=L_n(t)+K_nav^6-(n+26/9)K_nav^7.
\]

由此 av^r=2^(n+r)t^(1+r/n)。main 129–134 行将 av^6、av^7 换成 t 幂次时，两个 2 次幂和负系数完全正确。正 Θ_n(av^8) 与正 Θ_n(t^(1+8/n)) 在固定维数下等价。没有遗失规范化因子 2^(n+8)。

main 161–162 行的

\[
 J_n=\frac{(n-1)n^3(2609+1971n+81n^2)}{26244}
 +\frac8{27}\mathbf1_{\{n=2\}}
\]

与 R041/R046 接口一致。n≥2 时第一项已经严格为正，第二项非负，故 J_n>0 是解析事实，不由有限维数数值外推。

R041 对每个充分小 t 提供实际有限输入，具有 P7+J_n av^8+O_n(av^9) 的严格 halo 比值。因此 C_n 的相同归一化余量的下极限至少 J_n。R046 的全部输入上包给其上极限有限。main 166–172 行恰写成

\[
 0<J_n\le\liminf R_n(t),\qquad\limsup R_n(t)<\infty,
 \quad R_n(t)=\frac{C_n(1-t)-P_{7,n}(t)}{av^8}.
\]

这给全部充分小 t 的正双侧界：可取下界常数 J_n/2，并取某个固定的有限正上界。main 明确 Θ 为正双侧量级，并明确没有断言 R_n(t) 存在极限，未把具体构造的 J_n 称为最优八阶系数。

二维 K_2=4/27，(2+26/9)K_2=176/243；main 174–181 行的二维式为同一主结论的特例。J_2=15694/6561 的具体有限族系数仅用作下极限界，叙述准确。

## 3. ordinary j=0,…,7 展开的误差必须且确为 O8

main 138–144 行 b_(n,j) 与旧七阶系数完全相同，范围仍为 0≤j≤7。其一般 floor 公式继续由

\[
 L_n(t)-1=\sum_{k\ge1}v^{kn}(1+v)^{-kn}
\]

严格导出：v^(n+j) 的二项式次数 r=j-(k-1)n≥0，所以 k≤1+floor(j/n)，系数为

\[
 (-1)^r\binom{kn+r-1}{r}
 =(-1)^{j-(k-1)n}\binom{n+j-1}{j-(k-1)n}.
\]

因此全部低维碰撞仍被包括，j=6、7 的真实修正分别相加，其余 b 不变。

将解析的 L_n 截至 v^(n+7)，其剩余项为 O_n(v^(n+8))。这与主结论的正 Θ_n(av^8)=正 Θ_n(v^(n+8)) 同阶。两者相加一般只能给 O_n(v^(n+8))，因为完整 L_n 自身尚有八阶 Taylor 项；原余量的正性不能自动转移给普通截断多项式后的余量。

main 147–149 行准确采用

\[
 C_n(1-t)=1+\sum_{j=0}^7b_{n,j}t^{1+j/n}
 +O_n(t^{1+8/n}),
\]

没有把 ordinary 余项写成正 Θ8、o8 或未知八阶精确系数。main 453–456 行的证明也明确将两项八阶误差合并。主定理的正 Θ 只相对于完整非线性基线 P7，位置正确。

b_(3,5)=0、最初归一化极限以及第七系数全维非零的旧结论仍有效：新误差阶严格高于所有保留阶，没有改变这些系数或符号。三维第七系数为正，其余维数为负的限定保留。

## 4. 精确共享方体约化的条件未变

main 62–87 行保持同一半径、同一中心方体、公共交集 F 归一化、0<|F|<∞，并允许单个 E_i 无限测度。main 191–198 行的有限分支 λ_m^*<λ<1 和 m≥2 临界发散分支仍为原陈述。

实际字符串比较确认，307–336 行参数引理与 342–403 行精确约化证明均与旧 repo main 相同。此前独审的条件不被这次更新无声扩大：上包含依靠同见证总缺陷，反向等分使用可允许的无限测度 E_i，移动中心保留严格余量；临界发散仍在严格水平定义下成立。

对 0<ε≤1/2，参数引理给

\[
 \varepsilon\le s=\sigma_m(\varepsilon)
 \le-\log(1-\varepsilon)\le2\varepsilon,
 \qquad s\le\log2<1.
\]

故所有 m 同时在有限约化分支，有效 ε 区间不依赖 m，约化是

\[
 C_n^{(m)}(1-\varepsilon)=C_n(1-s).
\]

## 5. 全 m 的 ordinary O8 与 quartic σ 截断

令 p8=1+8/n≤5。ordinary 标量余项由第 3 节满足 |r_n(t)|≤C_n t^p8，对全部充分小 t 成立，常数和阈值只依赖固定 n。main 488–492 行对同一个标量余项代入 s∈[ε,2ε]，得到

\[
 \sup_{m\ge1}|r_n(\sigma_m(\varepsilon))|
 \le C_n2^{p8}\varepsilon^{p8}.
\]

这里没有从逐 m 的 O 常数推断一致性。

参数引理统一给 s=εA_m(ε)+O(ε^5)。两底数均在 [ε,Cε] 内，对八个固定 p=1+j/n∈[1,9/2]，均值定理给 O_p(ε^(p+4))。有限求和是 O_n(ε^5)。因为 n≥2，所以 p8≤5，这确为 O_n(ε^p8)。

最关键边界为 n=2，此时 p8=5，参数误差与目标余项同阶，仅可并入 big-O；不能称其为 o(ε^5)。main 498–499 行明确写“包括二维的相等阶情形”，且不据此识别八阶精确系数，强度正确。

本次还用本机 Python 标准库准确有理算术做了独立边界核验：m=2 时 σ 的 ε^5 系数是 7/128；n=2 的 b_(2,0)=4，故仅线性项替换的 ε^5 误差系数就是 7/32≠0。这直接展示 quartic 截断误差不能普遍升级为 little-o5。该计算只核参数截断边界，不核八阶全输入几何。

因此 main 480–482 行的全 m ordinary O_n(ε^(1+8/n)) 推论成立。三维 ε^(8/3) 项仍来自 192s^(5/3) 的修正 320a_m，等于 160(1-1/m)；八阶误差远高于这个保留阶，没有影响该示例。

## 6. 自然 σ 下的统一正 Θ 与 inf_m 下极限

定义前述同一个标量 R_n(t)。第 2 节的 liminf/limsup 给固定 c_n,C_n'>0、t_n>0，使所有 0<t<t_n 满足 c_n≤R_n(t)≤C_n'。取 ε<t_n/2，则全部 m 的 s≤2ε<t_n。精确约化因此同时给

\[
 c_n\le
 \frac{C_n^{(m)}(1-\varepsilon)-P_{7,n}(s)}{a_sv_s^8}
 \le C_n'.
\]

其中 a_sv_s^8=2^(n+8)s^(1+8/n)，main 519 行分母正确。因 m 一致的区间已明确，main 503–514 行的正 Θ 的两个常数及起效阈值都不依赖 m。

更精细的下极限传递也成立：给任意 δ>0，由标量 liminf≥J_n，可选 t_δ，使所有 0<t<t_δ 有 R_n(t)≥J_n-δ。ε<t_δ/2 时，这个不等式同时覆盖全部 s=σ_m(ε)，所以先取 inf_(m≥1) 仍保留 J_n-δ，再对 ε 取下极限并令 δ→0，恰为 main 511–512 行。

这是同一个标量函数在一个共同小区间上的控制，不是交换逐 m 极限与上、下确界。这里保持精确自然参数 s，没有用 quartic 近似替换 P7(s)，所以第 5 节二维的同阶参数误差不会破坏正 Θ 或归一化下界 J_n。

## 7. 删除旧弱命题及基础依赖

旧 `prop:nonexact` 和 `eq:nonexact` 是只给六阶下极限的弱下界，其结论已被 R041 的实际八阶有限构造覆盖。它不提供 global/sharp 上界或精确约化所需的基础引理；保留原周期孔下界 `prop:lower`、有限见证、同时等分、参数引理和上界模块，数学依赖没有被删除。

实际扫描新 main 与不变 global/sharp，`prop:nonexact`、`eq:nonexact`、`cor:nonexact-multi`、`cor:uniform-planar`、`planar:`、`skew:`、`seventh:` 出现次数均为 0。因此没有旧弱命题被正文引用而失去定义，也没有新证明通过该命题循环回到待证下界。

main 180–181 行用实际有界输入的正下界直接承接新主定理，main 447 行引用的新接口拼写为 `eighth:theorem`，与作者承诺一致，并非 `eith:theorem`。新模块按约定不再保留旧弱命题补证即可。

## 8. 实际静态检查和计算证据

因新模块尚未冻结，本次只统计完整读取的 main/global/sharp 三文件，而不虚报全链计数：

```text
main/global/sharp labels: 138 refs: 159
duplicate labels: {}
unresolved (new module intentionally absent): ['eighth:theorem']
new main inputs: ['sections/global', 'sections/sharp', 'sections/eighth']
unchanged block parameter lemma True
unchanged block reduction proof True
```

唯一待供接口是明确预期的 `eighth:theorem`，不因暂未合入模块判为错误。完整 label/ref 总数与编译验收须在最终模块冻结后核对。

首次尝试当前 `python` 环境的 SymPy 未成功（`ModuleNotFoundError`），随后立即采用已可用的 Python 标准库 Fraction 完成准确核验，没有安装依赖或把失败调用当作计算证据。实际 stdout 为：

```text
arithmetic: Python standard-library Fraction (exact)
J2: 15694/6561
J2 equals expected: True
m2 sigma epsilon^5 coefficient: 7/128
n2,m2 leading b20*s substitution error coefficient: 7/32
n2 p8: 5
```

这些计算仅核抄录分数及 quartic 参数截断的边界阶；R041/R046 的集合几何、全输入上确界和本次新模块独审仍由相应严格证明负责。

**最终签署：main 冻结 SHA `b976cf01…b63f5` 的正八阶 P7 结论、ordinary O8 系数展开、全 m 传递、自然参数正 Θ 及 inf_m 下极限、旧弱命题删除均 PROOF-VALID，无必要数学修补。新 `eighth:theorem` 的完整冻结模块及最终静态组装/编译由后续接口验收绑定。**
