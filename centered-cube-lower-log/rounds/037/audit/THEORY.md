# 分位数分类：真实PB似然比与R级点态构型

先运行本机Fraction/真实box权重小例，再提出下述解析结论。quantile_probe.py保存12组n/R控制；没有Gaussian、尾MC或原匹配来源假设。结论是点态似然倍率可以随R增长，同时明确其尚不产生weak增益。

## 1. 质量与后验的准确表达

K偶，plus为前K/2格。第b表的每格箱质量1/K，故产品来源f_b下plus计数S先验为Bin(n,1/2)。令I_b为固定整数分位数段，pi_b=P(S∈I_b)，保留空段和实际格点质量。joint来源为F=sum_b 1_(S∈I_b)f_b，质量sum pi_b=1；prior为bar f=sum pi_b f_b，也质量1。

固定同一query x、同一L，记z_i,b^±为窗口内plus/minus箱质量除L，s_i,b=z^++z^-、A_b=prod_i s_i,b。当A_b>0时p_i,b=z^+/s，T_b=P_(PB(p_1,b,...,p_n,b))(S∈I_b)。则

    A_LF=sum_b A_b T_b,       A_L bar f=sum_b pi_b A_b.         (1)

只对pi_b>0定义ell_b=T_b/pi_b，omega_b=pi_b A_b/sum pi A。准确有

    A_LF/(A_L bar f)=sum_b omega_b ell_b.                      (2)

没有额外pi，也不能把表平均移入各坐标。空段T_b=0，直接跳过该项即可。

Q_b,L在query空间的密度为A_b，质量1。Fubini及cube核保积分给

    E_(Q_b,L)T_b=pi_b,       E_(Q_b,L)ell_b=1,
    Q_b,L(ell_b>u)≤1/u.                                      (3)

ell是条件gated来源cube平均相对unconditional表平均的真实Bayes因子。它最大可达1/pi_b；均值恒等式说明大倍率仍需支付其query概率，但(3)自身不是maximal weak上界。

对L=q+ell≥1，整周期plus/minus各q/2、余质量0≤E≤1，故1/4≤p_i,b≤3/4。这个范围并不把PB分位段概率锁在常数倍pi内；plus计数分段不是模类别的非零字符，不满足第36轮模分类的指数退化机制。

## 2. 已运行的合法R级共同尺度构型

取K4，theta_b=(b+1/2)/R，plus格原子offset theta_b、minus格offset 1−theta_b，delta=1/(64KR)。所有箱严格留格。令b*=R−1、L*=1+2delta。固定query坐标为top表任一plus原子减1/2（mod1）。额外窗口全长2delta，完整包含该箱，并与其他表箱不交。

每坐标窗口已包含一个整周期。top表多得到一个plus箱质量1/K，其他表无额外质量。因此

    p_i,b*=3/5,   p_i,b=1/2 (b≠b*),
    A_b*/A_b=(5/4)^n=:a.

设u=P_(Bin(n,3/5))(S∈I_b*)、pi=pi_b*。由(1)，真实同一L比值为

    r=(a u+1−pi)/(a pi+1−pi).                                (4)

这些不是假设的PB参数：脚本对每表、三个有理偏移逐个计算periodic box积分，全部准确核对。n1023结果为：

| R | 实际top先验pi | top PB概率u | 同一L joint/prior |
|---:|---:|---:|---:|
| 4 | .24578737 | .9999999961 | 4.0685573 |
| 8 | .11739307 | .9999999321 | 8.5183898 |
| 16 | .05897370 | .9999995054 | 16.9567021 |
| 32 | .03030695 | .9999976674 | 32.9956574 |

倍率不是恰R，因为pi不恰1/R。若额外窗口改成对所有表都给相同plus质量的完整格段，所有A_b、PB分布反而相同，sum_b T_b=1，joint/prior准确为1。这给出另一处真实抵消。

## 3. 正体积plateau与完整M量词

每个top-plus原子的额外窗口中心偏移允许区间为[−delta/2,delta/2]，**全宽delta**。任意相交的中心区间才为3delta，二者不能混用。K4共有K/2=2个plus原子，定义Q为每坐标都落在这两个完整包含区间之一的product集合。各区间互不相交，准确有

    |Q|=[(K/2)delta]^n=(2delta)^n=(1/(128R))^n.                 (5)

因此(4)在Q的全部内部同时成立，且用同一个L*。固定L*的平均值F*也在Q上常数，其限制到Q的尾系数sup_lambda lambda |Q∩{A_L*F>lambda}|准确为F*|Q|，不超过

    2[((5/4)/(1+2delta))/(128R)]^n.                           (6)

这是**仅固定L***的局部贡献帽，不能用于全尺度M。

事实上半周期把top表plus原子准确映成table0的minus原子，因为1−theta_0=theta_b*。Q就是table0、所有source cell均为minus的箱支持。该cell计数S=0归入I0。joint密度在Q内等于峰值H=(64R)^n，所有cube平均均≤H，而足够小L达到H。因此对全正尺度M，

    MF=H  a.e. on Q,
    sup_lambda lambda |Q∩{MF>lambda}|=H|Q|=2^-n.              (7)

n1023/R32时(6)的log帽约−8280，而(7)的log准确为−1023 log2≈−709.09。后者仍只计算Q的贡献，绝不以上界整个来源的weak。对有限L网格未必能在Q全部达到H，不能把(7)原封不动声称有限grid等式。

prior各表箱支持两两不交，故其密度峰值为pi_max H；Q上prior局部密度为pi_0 H。全尺度上进一步有

    1/pi_max ≤ MF/Mbar f ≤1/pi_0  on Q.                       (8)

实际n1023/R32给下端约20.0971；因此固定L的32.9957也不能直接充当两个max的比值。

## 4. 不存在一般常数点态帽，但没有weak增长定理

这个例子确可随R增长，而非只到32。取n=2R^4+1。分位数定义给pi_top≤1/R。对Bin(n,1/2)用Chebyshev，可取top cut≤ceil(n/2+sqrt(nR))，除n趋于1/2；故Bin(n,3/5)落入top段的u→1。a=(5/4)^n≫R，由(4)

    r≥u/[1/R+1/a]~R.

所以这类合法来源在L>1上的相对点态帽不能统一为常数。

全尺度点态倍率也可达同阶。令t_m=C(2m,m)/4^m，递推t_(m+1)/t_m=(2m+1)/(2m+2)以及准确不等式

    [(2m+1)/(2m+2)]²≤(3m+1)/(3m+4)

归纳给t_m≤1/sqrt(3m+1)。n=2R^4+1时最大Bin(n,1/2)原子≤1/sqrt(3R^4+1)=O(R^-2)。各分位数段质量与1/R之差至多该原子质量，因此pi_max≤1/R+O(R^-2)，由(8)得Q上MF/Mprior≥(1−o(1))R。脚本实际核对了65个该有理平方递推，不以有限检验代替上述归纳。

但源密度本身已有完整归一化抵消：joint是高度H、支持体积1/H的平顶密度，W(F)=1。prior在R份互不相交、各体积1/H的支持上高度pi_b H，故

    W(bar f)=max_(k=1..R) k pi_(k) ≤1,

pi_(k)为降序先验；在上述n/R族它也趋于1。L1恒给Mprior的weak≥1。点态R倍可以由峰高度与支持体积互换支付，既不证明joint最大weak更大，也不排除其他query尺度和区域产生真正收益。仍需在原律下控制整个共同尺度尾，而非把(6)/(7)的局部区域当全局weak上界。
