# Fine log-score与coarse count的精确接口

本稿写于实际运行 `correlation_probe.py` 与 `mixed_tilt_probe.py` 后。仅在真实product query通道上推导恒等式，不调用Gaussian/CLT独立结论。

## 1. 基本量与相关系数

固定一个表、一个共同L。令a(x)=z_+(x)+z_-(x)为单位质量query密度，Q(dx)=a(x)dx，p(x)=z_+(x)/a(x)，η=∫z_+。以binary源类别C与query X的真实Bayes coupling解释p=P(C=1|X)。假设以下score矩有限，记

`ℓ=log a`，`μ=E_Qℓ`，`u=η(1−η)`，`f=Var_Q(p)/u`，`c=Cov_Q(ℓ,p)`。

全期望给E_Qp=η，且

\[
c=\operatorname{Cov}(\ell,C)
=u\{E(\ell\mid C=1)-E(\ell\mid C=0)\}. \tag{1}
\]

n坐标独立复制时，V=Σ(ℓ_i−μ)与posterior mean count Σp_i的相关系数恰等于单坐标corr(ℓ,p)；同一个V与**潜在真实源count** S=ΣC_i的相关系数为corr(ℓ,p)√f。这些是有限n恒等式，不是高维分布的独立性。

进一步，条件于完整labels后各query独立，且score条件均值差由(1)给c/u。因此

\[
E(V\mid S=s)={c\over u}(s-n\eta),\qquad
E(V\mid S\in E)={c\over u}E(S-n\eta\mid S\in E). \tag{2}
\]

它描述count gate导致的fine score**一阶均值位移**，没有控制其条件方差、尾分布或共同尺度最大值。

## 2. 对count后验高阶chaos的线性正交性

用39轮的记号g=(p−η)/√u，

`v_k=binom(n,k)^(-1/2) Σ_{|D|=k}∏_{i∈D}g(X_i)`。

对k≥2，展开E(Vv_k)时每项至少含一个只出现一次的g坐标，其均值0，故

\[
E(Vv_k)=0\ (k\ge2),\qquad E(Vv_1)=\sqrt n\,c/\sqrt u. \tag{3}
\]

所以线性fine log-score仅与count posterior的第一阶成分有耦合；它不意味着exp(V)、尾indicator或多个L的最大值与高阶后验无关。若E为count段，令L_E=P(S∈E|X)/P(S∈E)，则由(2)或正交展开

`E_Q[V L_E]=(c/u) E_prior[S−nη|S∈E]`。

这里Q是完整单位质量product query倾斜。另行条件于safe query事件后，坐标独立性及E g=0可能丢失，不能免费沿用。

## 3. 精确混合mgf与真实类别漂移

对使m(t)=E_Q e^{tℓ}有限的t，定义

\[
\eta_t={E_Q[e^{t\ell}p]\over m(t)}.
\]

用e^{tℓ}/m(t)倾斜单坐标query，Bayes posterior仍是p(X)，但倾斜后的latent binary边际成为Bernoulli(η_t)。坐标仍独立。因此对任意count集合E，

\[
\boxed{E_{Q^n}\left[e^{t\sum_i\ell_i}P(S\in E\mid X)\right]
=m(t)^n\,P\{\operatorname{Bin}(n,\eta_t)\in E\}.} \tag{4}
\]

对中心化score只需两边乘e^-tnμ。此式是有限n精确乘积恒等式，没有PB/normal替换。也可逐坐标写成生成多项式

`E_Q[e^{tℓ}(1−p+wp)]=m(t)(1−η_t+wη_t)`

后取n次方及对应正系数求和。`mixed_tilt_probe.py`已对n3、t=-1/0/1/2及有理正权重/后验参数直接枚举验证(4)，四例完全相等。

η_t是这个额外query倾斜中的类别边际；原来源的实际prior仍是Bin(n,η)。不能把(4)中的η_t直接替换进原先未倾斜的source质量或importance公式。

在可微分换序条件下，η'_t=Cov_{Q_t}(p,ℓ)，特别η'_0=c。当前真实macro来源a被正的常数上下包，故所有实t均可用，所需导数换序由有界性保证。要认证某段fine saddlepoint范围下的类别漂移，可直接对真实一维核计算η_t的区间；仅报告c的小数值不足以替代整个t区间。

## 4. 零协方差不能升级为独立

实际Fraction小例取Q在三点等概率，p=(9/20,3/5,9/20)，正score权重w=(4/5,1,5/4)。p均值η=1/2，而log w=(-h,0,h)，所以Cov_Q(log w,p)=0严格成立。按w倾斜后

`η_1=Σwp/Σw=609/1220 ≠1/2`。

因此即使一阶相关严格为0，混合mgf里的类别参数仍可改变，fine score与latent count也未必独立。这是抽象有限binary channel反例，用于说明逻辑限制，未把它冒充某个已构造cube来源。

## 5. 已运行真实核诊断的适用范围

38来源L7/6四表的连续点核Gauss2积分，corr(ℓ,p)约−.05502、.02521、.03764、.001003；真实潜在count的相关再乘√f≈.126。对t∈[-2,4]所测η_t，最大标准化count均值位移绝对值约.0757。它们是有用的近正交诊断，但积分和完整FP尚未向外认证。

新general-R来源在R4/8/16/32、L10/9/7/6/11/9的180个bin-channel读数中，最大|corr|由约.063增至约.091，未达到严格正交，也没有从有限R增长拟合极限。R之间midpoint相位网格不同，来源族不嵌套；不能把变化读成单个来源上逐次加表的单调性。

后续若以(4)为分析入口，应保留实际η_t、真实prior段质量、同一个L和来源权重。该接口没有证明full tail、Gaussian collapse、weak排序或Cn新对数收益。
