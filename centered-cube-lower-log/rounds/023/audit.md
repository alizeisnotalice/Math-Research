# 第23轮：整样本 winner 混合重要抽样审查

结论：所述整样本混合恒等式成立；round023/probe.py 的 winner 选择位置、混合分母和共用子集逻辑符合该恒等式。理论权重界为 r exp(n max F_j)，不是无条件的 r。当前结果仍是固定正权重求积原律的有限样本估计。

## 混合密度与估计量

设 sum_x w_x=1、w_x>0、Y_j(x)>0，Z_j=log Y_j，F_j=log sum_x w_xY_j。令 q_j(x)=w_xY_j(x)exp(-F_j)。先为整组样本均匀选择一个 J，再从 q_J 独立抽取 n 个坐标，则相对 P^n 的密度为

D(x_1,...,x_n)=(1/r)sum_{j=1}^r exp(S_j-nF_j)，S_j=sum_i Z_j(x_i)。

因为所有 Y_j>0，混合律与原律具有相同支持。取 A_I={max_{j in I} S_j>h}、lambda=exp(h)，严格有

lambda P^n(A_I)=E_mix[exp(h)1_{A_I}/D]。

它不等于每个坐标重新选择 winner 的提议律；后者密度是 product_i[(1/r)sum_j exp(Z_j(x_i)-F_j)]。二者不可交换。条件于单一 winner 的 multinomial(n,q_J) 计数抽样与 n 个独立坐标完全等价；之后仍须计算全部 S_j 和全部 r 项的混合 D，不能只使用已选 winner 的项。

任意嵌套子集 I 都可使用同一批 r 分量混合样本。每个事件仅修改其指示函数；分母 D 一直保留全部 r 项，不能改成该子集自己的混合密度。所有坐标必须来自同一来源向量、相同 n/a 和共同 h。这样每个子集的估计均无偏，且由于事件嵌套、权重非负，各复制及其样本均值都随子集扩张单调。该单调性本身是结构事实。

## 权重界与非零 F

若事件 A_I 发生，选取其中一个 S_j>h 的 j，则

D>(1/r)exp(h-nF_j)，

0<=exp(h)1_{A_I}/D<r exp(nF_j)<=r exp(n max_{j in I}F_j)。

因此可统一采用 B=r exp(n max_{1<=j<=r}F_j)。即使子集只有 k 个尺度，该共同提议律下的界仍有 r，而非 k。只有相关 F_j 均 <=0 时可安全放松为 r；仅“接近零”不够，必须控制 nF_j。F_j 非零时把它强行改零也会破坏混合密度恒等式。

当前结果记录 nF 范围为 [-9.22414408418276e-10,2.855092216532759e-9]。由这些浮点展示数得到 B 约为 65.00000018558099，数值上接近 65，但不能据此声称严格上界恰为 65。理论 F_j 应取所用正权重律的精确归一化常数；当前代码 F_j=log1p(a sum wX_j) 的代数表达正确。实际浮点 q 再归一化及 n 倍舍入误差尚未认证，机器计算不能据此称为严格算术无偏证明。

## 固定 h=30 与 pilot 口径

只读代码显示本轮 h=30、r=65、复制数 256 在抽样前冻结；每复制重新独立选一个 winner，之后生成一整组计数。3/9/17/33/65 的索引确为嵌套子集，且分母保持相同 65 分量。这不是用本轮样本再挑最优阈值；各固定子集估计目标明确。

若先用独立 pilot 选择 h，再冻结并用独立 final 样本评估，则条件于 pilot，仍对被选定的 lambda P(A) 无偏。如果使用同一批样本搜索 h 后报告最大估计值，就不能把最大值当作无选择偏差的固定阈值估计。固定阈值的重要抽样恒等式也不直接估计对所有 h 的总体上确界。

对固定 I，令 U_m=exp(h)1_{A_I}/D，M 个独立复制的标准误通常估为 sd(U,ddof=1)/sqrt(M)，代码使用此式正确。不同子集的估计存在共同样本相关性；比较增量应使用配对差值。所列 SE 与普通正态区间不是严格概率下界，尤其小子集命中数少时。当前 order=2 求积实验没有支付连续来源误差，也不能从有限尺度点估计推出尺度数趋于无穷的结论。

## 可复现有限小例

本机 Python / NumPy 2.5.4 实际运行以下程序，枚举 3 状态、n=3 的全部 27 个有序样本：

```python
import itertools
import numpy as np
w = np.array([.2, .3, .5])
Y = np.array([[1.6, .8, 1.1], [.9, 1.3, .7], [.7, 1.1, 1.2]])
Z = np.log(Y)
F = np.log(w @ Y)
q = w[:, None] * Y / np.exp(F)
n, r, h = 3, 3, .3
subsets = [(0,), (0, 1), (0, 1, 2)]
direct, correct, wrong = np.zeros(3), np.zeros(3), np.zeros(3)
for xs in itertools.product(range(3), repeat=n):
    xs = list(xs)
    S = Z[xs].sum(axis=0)
    p = np.prod(w[xs])
    proposal = np.mean(np.prod(q[xs], axis=0))
    D = np.mean(np.exp(S - n * F))
    assert abs(proposal - p * D) < 1e-15
    # Deliberately wrong denominator for this global-winner proposal:
    wrong_D = np.prod(np.mean(q[xs] / w[xs, None], axis=1))
    for k, ids in enumerate(subsets):
        hit = S[list(ids)].max() > h
        direct[k] += np.exp(h) * p * hit
        correct[k] += proposal * np.exp(h) * hit / D
        wrong[k] += proposal * np.exp(h) * hit / wrong_D
print(F, direct, correct, wrong, sep="\n")
```

F=(-0.06187540371808741,0.09531017980432493,0.02955880224154443)，刻意包含非零归一化量。直接原律目标为

(0.14038531598790435,0.7356730501289217,1.1068842222123227)，

正确混合 IS 枚举为

(0.14038531598790435,0.7356730501289218,1.106884222212323)，

最大 float64 差为 2.220446049250313e-16。若误用逐坐标混合的分母，则变成

(0.15002840127748215,0.7515987860155989,1.1282384172769178)，

明确改变目标。这是独立数值交叉核验；严格恒等式依据是前面的密度推导。
