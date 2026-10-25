# audit_I028：I027三个真实R负矩cap资格

2026-10-26。**通过。** 同一真实形状概率R下有

|gamma|严格可用上cap|
|---|---:|
|1/5|E_R d^(-1/5)<926346056693/500000000000|
|3/14|E_R d^(-3/14)<2060603220817/1000000000000|
|1/4|E_R d^(-1/4)<1515346060487/500000000000|

只有gamma=3/14与既有同指数Q作改善比较：新有理cap严格低于旧完整球下端。另两项是新的真实矩输入，不与不同指数的旧cap混比。**本批没有新I上界或C5界**；可把这些矩输入送入另行资格审查的SIP/majorant候选，不能省略其点态和全域验证。

固定句柄 `audit_I028_fixed_001`，PID **35290** 在计算开始前打印并保存在结果中。一次批次仅三个新点，实际约0.000542秒；不重跑022/024几何或旧Gamma公式，不调用MATLAB，不提交、不改其他锁目录。

## 1. 输入锁与两个不同Z

027五件白名单及全部来源逐项核hash，manifest为 `27f69351a186acfe18cf3071943933350a88de912d482e60a3dd41e60b96c867`。024 manifest和保存的Zflag准确有理区间另行核对，不能用旧显示小数重建。

模型保持16条未单位化body方向、独立alpha=3/2 stable、Lebesgue。原始S_i为独立正beta=3/4 stable，A=sum S_i v_i v_i^T，T=sum S_i，M=A/T，d=det M。R以T^-beta倾斜原概率再推送形状，

\[
 Z_{tilt}=ET^{-\beta}=1/[16\Gamma(7/4)].
\]

这是真正的概率归一化。它与024的几何代数常数
`Zflag=sum_f theta_f R_f^(1/beta)`不同；本批没有混用。

## 2. 原始detA负矩：最优权与gamma无关

已审022分解给det A/W=sum theta_f F_f，W=16^5，每个F是五个独立stable因子的规范化乘积，三个合并线性和的尺度积为R_f。对任意gamma>0，所有原始负矩存在且

\[
 EF_f^{-\gamma}=m(\gamma)^5R_f^{-\gamma/\beta},\quad
 m(\gamma)=\Gamma(\gamma/\beta)/(\beta\Gamma(\gamma)).
\]

对任意确定flag权q，以同一熵AM–GM和普通Holder得到log预算中的可调项

\[
 \gamma\bigl[\mathrm{KL}(q\Vert\theta)
 -(1/\beta)\sum q_f\log R_f\bigr].
\]

因此唯一最优 `q*_f=theta_f R_f^(1/beta)/Zflag`不依赖gamma。固定几何和Zflag可直接承接到三点，不需要重新优化或生产。得到

\[
 E(\det A)^{-\gamma}
 \le16^{-5\gamma}m(\gamma)^5Zflag^{-\gamma}.
\]

该原始S概率下的cap对全部gamma>0有效。后面的R转换还有额外条件，不能把两者适用范围混同。

## 3. R转换与Harris方向

由det M=det A/T^5及R的倾斜定义，精确有

\[
 E_Rd^{-\gamma}=Z_{tilt}^{-1}
 E[T^k(\det A)^{-\gamma}],\qquad k=5\gamma-\beta.
\]

若k>=0，T^k对每个S_i递增；det A由非负Cauchy–Binet单项式组成，对每个S_i递增，所以其负gamma次幂递减。原产品概率上的反向Harris给上界方向

\[
 E[T^k(\det A)^{-\gamma}]
 \le ET^k\ E(\det A)^{-\gamma}.
\]

方向可独立证明：单变量两副本的积差 `(f(S)-f(S'))(g(S)-g(S'))<=0`；逐坐标条件化归纳保持一增一减。两非负函数分别取min(N,·)，先用有界不等式，再单调极限。右端由独立矩公式有限，因而没有循环假设联合矩有限。并未假定T/P独立，也未在倾斜R下使用产品概率。

正矩存在需 `0<k<beta`，其Mellin积分给

\[
 ET^k=16^{k/\beta}\Gamma(1-k/\beta)/\Gamma(1-k).
\]

在k=0，T^k=1，恒等处理即可。故本R-cap路径的准确范围为 `3/20<=gamma<3/10`。gamma<3/20时两个函数同向递减，这条Harris比较不再给所需上界；gamma>=3/10时所用正矩发散。后者不证明真实R负矩必然发散，只说明本条乘积分离不能产生有限cap。

## 4. 尺度及固定公式

保留三个来源的16幂：倾斜逆归一化16、正矩16^(k/beta)、detA负矩16^(-5gamma)。准确总指数

\[
 1+k/\beta-5\gamma=5\gamma(1/\beta-1)=5\gamma/3.
\]

因此

\[
 Q_\gamma=\Gamma(7/4)16^{5\gamma/3}
 \frac{\Gamma(1-k/\beta)}{\Gamma(1-k)}
 m(\gamma)^5Zflag^{-\gamma},\qquad E_Rd^{-\gamma}\le Q_\gamma.
\]

三点的k、k/beta、净指数依次为

|gamma|k|k/beta|净16指数|
|---|---|---|---|
|1/5|1/4|1/3|1/3|
|3/14|9/28|3/7|5/14|
|1/4|1/2|2/3|5/12|

它们严格位于Harris正矩条件内部。gamma=3/20的符号端点还原024公式；本批未数值复跑这一旧Gamma基线。

## 5. 不同数值路径与完整球比较

作者把净16幂合并、使用正参数Gamma比。本批保留上述三块尺度，倾斜逆常数用 `16 beta Gamma(beta)`，负矩使用未移位Mellin式，并把正矩比例改写为

\[
 \frac{\Gamma(-k/\beta)}{\beta\Gamma(-k)}.
\]

这是Gamma递推给出的同一正值；本次两个负Gamma参数都在(-1,0)，不触及极点，数值包围也核其比值严格正。k=0端点不使用这个带极点的表示。

以Arb512位算正prefactor，再用Zflag上端生成Q下端、Zflag下端生成Q上端，方向随负gamma次幂单调性确定。三个新完整有理区间均**包含于作者保存的完整有理区间**，不是只比较显示中心；并严格低于表列cap。50位外向显示分别为

\[
 Q_{1/5}\in[1.85269211338500696138601415910284796776225617324733,
             1.85269211338500696138601415910284796776225617324734],
\]
\[
 Q_{3/14}\in[2.06060322081644261752091448451310557804816335235458,
              2.06060322081644261752091448451310557804816335235459],
\]
\[
 Q_{1/4}\in[3.03069212097398446859442799463492845657412770031681,
             3.03069212097398446859442799463492845657412770031682].
\]

同指数旧Q来自锁定round011结果，完整球约2.45763054393722847；未重新计算其旧Gamma公式。新Q3/14的**严格有理cap本身**已低于旧球下端，比较方向充分。

Gamma/Arb底层库仍共享；独立性是尺度分解、不同正矩Gamma参数、准确条件/端点方向及全部区间包含，不是跨库Gamma认证。Q公式的双侧包围不是实际E_Rd^-gamma的下界。

## 6. 可消费资格、运行及白名单

三项cap可作为同一真实R下的矩约束供后续新候选使用。旧majorant族排除中的抽象分布必须满足旧负矩约束；cap变小后必须重新核可行性，不能自动迁移旧排除，也不能反过来宣称新cap已解决majorant或I。

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/audit_I028/check.py --source-root /path/to/audit-cont-repository --model-root /path/to/model-repository --run
```

默认只核hash。固定PID35290批次512位、三个新点、约0.000542秒，未复算几何、旧Q公式、旧g、I或SIP。精简白名单仅registration.md、check.py、result.json、本报告，未提交。
