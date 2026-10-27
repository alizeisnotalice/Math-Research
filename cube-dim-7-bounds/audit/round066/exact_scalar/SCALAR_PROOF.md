# R066 · 同一完整Phi的精确log价格推论（conditional）

## 身份与结论

本附录只审核准确scalar代数，不重审operator链。冻结作者源为
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md`，SHA256
`4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1`。

固定L=log n、T=nL−nlogL、gamma=317/1000、sigma=83/1000、x=(917/1000)L及t=n/x。**条件于该源的完整算子定理被接受**，本附录证明

C_n≤Phi66:=F(A_n(t),B_n(t)−gamma n)<T−(23/100)n

对所有充分大的整数n成立。没有给出显式N。这里Phi66与作者原式完全相同；较强的易读显示推论来自保留精确log价格，不是另一kernel、另一缺额或同一Phi函数的新优化。

## 一条准确对数包络

令u=83/917∈(0,1)，并令

S_j=Σ_{k=1}^j (−1)^(k+1)u^k/k。

积分有限几何恒等式给

log(1+u)−S_7=−∫_0^u t^7/(1+t)dt<0，
log(1+u)−S_6= ∫_0^u t^6/(1+t)dt>0。

所以S6<price:=log(1000/917)<S7。本轮新Fraction核验实际给

S7<4333/50000。

这是严格有理上包，不从显示小数重建价格，也没有调用作者的log/16项验证程序。于是准确有

gamma−4333/50000−23/100=17/50000>0。

## 完整根式的余项支付

作者原scalar界已经含有精确价格，尚未作sigma/(1−sigma)的放松：

B_n(t)≤T−nlog(1−sigma)+1=T+price*n+1，
A_n(t)≤n exp(−sigma L)，B_n(t)≤nL。

在原Phi66正域中，B_n−gamma n≤B_n，故

Phi66=(B_n−gamma n)+A_n+2sqrt[A_n(B_n−gamma n)]
 ≤T−c_star*n+n R(L)，

c_star:=gamma−price=gamma+log(917/1000)，
R(L):=exp(−L)+exp(−sigma L)+2sqrtL exp(−sigma L/2)→0。

这保留A、cross及1/n全部正cost。任取固定实数c<c_star，余量c_star−c>0，因而在完整operator定理自身的eventual有效域中再取充分大n使R(L)<c_star−c，即得Cn≤Phi66<T−cn。这里不声称端点c=c_star本身。

具体地，准确包络给c_star>23/100+17/50000；取R(L)<17/50000便得上述23/100推论。它严格强于作者原9/40显示推论，但比较的是同一完整cost的两个scalar上包，不宣称Phi66自身被进一步降低。

原一般L1输入、全部levels、所有正共同radius的桥接完全依赖R066 full审，本附录不以scalar核验替代那项责任。作者原source和receipt逐字不改；这次新receipt只记录本附录的新七项有理计算与条件性结论。
