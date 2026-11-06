# 第17轮可接续推导（未审定，勿作已证成果）
2026-11-06。r015/r016已完成；完整历史见仓库报告。
外部实际运行：sm_soft_shift_probe.m、soft_shift_probe.json（MATLAB约27秒，n256…1e7）。

## 目标
将r016的真正分位数收益升级为实际指数软轮廓的下界评分收益。
定义A_c(q)=E exp((S_n−q)/n)1S_n≥q，唯一qsoft(c)>0满足A=1。
候选c_n=mu−sigma ρ/2，ρ=√(logn/n)。期待qsoft(c_n)−qsoft(mu)=(mu−1)logn/8+o(logn)。
J=e^-q I_n，I_n与c无关，故同量级改善可转到实际下界评分。
严格可行通过qsoft−1/n、固定n有界截断再深层极限实现；不宣称Q全局最优或C的新logn系数。

## 证明路线
1. A严格递减连续，从A(0)=E exp(S/n)>1至0；n>1确保MGF有限。
2. 写d=nm(c)−q。未截断全均值exp(nK_c(1/n)−q/n)=1+d/n+O(L/n)。
在d=(1±eps)sigma√nL，通过FGL粗尾比较可证qsoft位置d~sigma√nL均匀于c=mu+tρ,|t|≤T。下侧d大：A≥full−F>1；d小时F≫ρ且低尾带指数权重至少F−E(q−S)+/n≥F−O1/√n，故A<1。
3. root处低尾加权与普通CDF差≤E(q−S)+/n=oρ（固定Vmax负中心mgf，r014(4)）。所以F(qsoft)~d/n~sigmaρ。
4. r016统一Fourier引理仍可用：N=n−1，a=q−c，负倾斜λ使N倾斜均值=a，则λ~ρ/sigma。
F_N(a)=M^N e^λa ∫_{−∞}0 e^λy g_N(y)dy。
用||g_N'||∞=O1/N及g_N(0)~const/√N比较积分与g0/λ，relative error O1/(λ√N)=O1/√L。故f_N/F_N~λ。
CDF Taylor保证F_N(a)/F_n(q)=1+O(ρ²)，f_n(q)/f_N(a)=1+Oρ²。于是f_N(a)~f_n(q)~L/n。
5. signed kernel derivative的精确等式对任意积分目标成立，给∂c A=np_c[A−B_N(q−c)]，B_N(x)=E exp((S_N−x)/n)1S_N≥x。
∂q A=−A/n−f_n(q)。因此softq'=np_c(A−B)/(1/n+f_n)在root。
6. B'=-B/n-f_N, B''=B/n²+f_N/n-f_N', B'''=-B/n³-f_N/n²+f_N'/n-f_N''。
用r016 f_N'/f_N=λ+O1/n及全局f_N''(s)/f_N(a)≤C e^{λ(s−a)}(λ²+1/n)，且B_N(s)≤untruncated raw positiveMGF≤C e^{|s−a|/n}。
Taylor A−B=(c−m)B'+E(Z−c)²/2 B''+R。
R=O(f_N(a)(λ²+1/n)+1/n³)。root B_N(a)=O1，f_N~L/n，所以B'=-f_N(1+O1/L)，B''=-λf_N(1+o1)。
于是A−B=−f_N[(c−mu)+mu λ/2]+o(f_N ρ)。
分母f_n+1/n=f_N(1+o1)。得到qsoft'/(nρ)→−p_*[t+sigma/2]，积分即所需收益。

## 本机数值
指定偏移soft q gain/logn：n256 −.00999261；1024 −.00389254；4096 .000564023；1e4 .00287604；1e5 .00753211；1e6 .01094897；1e7 .01357600。
理论候选极限(mu−1)/8=.034808…，小n负收益必须保留。
数值使用正倾斜精确Gamma-binomial公式，稳定logmgf=-log1p(-1/n)+log1p(p_c*expm1(c/n))。
q以55步二分，binomial仅众数±12sd（未定向认证）。
logI用稳定Stirling表达.5log(2πn)−1−nlog1p(-1/n)+log Q(n,n−1)+1/(12n)−1/(360n³)。
后续应Wolfram60/100位独立复核正倾斜公式与I_n，明确Stirling余项，Analyzer检查，再写r017并独立审核才提交。

## 更远线索（纯预测）
固定mu的soft q可能nmu−d_n+(3−mu)L/6+oL（第三累积量mu(3−mu)）。logI~L/2+const，故固定mu的实际ramp score可能nmu−d_n−mu L/6+oL。
动作偏移改善(mu−1)L/8后，score可能nmu−d_n−(mu+3)L/24+oL。
这些尚未证明；要用足够精确的一阶Edgeworth尾或参数导数，不可直接登记。C/Q的真正多轮廓或反馈logn项仍未知。
