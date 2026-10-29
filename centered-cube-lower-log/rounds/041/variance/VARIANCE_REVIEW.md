# 当前uniform-bin方差改进独审

结论：**数学推导通过**，适用于当前真实source的连续uniform隐藏query-bin及已列均值误差前提。它改进的是variance上帽，不改变来源、真实方差或原40轮结果。未读取/变换tail payload，未重跑MC，未生成新CI。

已只读审查 `../round040_source/CANDIDATE_VARIANCE.md`、`variance_probe.py`、`CANDIDATE_VARIANCE_RESULTS.json`，source/script哈希匹配。另实际用Fraction复核32个保存有理case的 `old=3new`、`new=Lip²h²/12`、`Var≤new`，其中10个非零affine例等号；24个common-D数值例最大Var/new≈.9938831564，仅作普通FP控制。340组非负有理产品参数及16个Cantelli恒等式小例也通过；没有新随机样本。

## 1. 当前range确实是Lip×h

当前h=δ/6≤δ，原range的min(1,h/δ)未饱和。N_partial使用窗口端点的(δ+h)/2扩大候选带，覆盖整bin中所有可能非零的peak斜率；不能只看bin中心。class内部gap证书保证每端≤1，总N≤2。

class核的a.e.导数分别由常数背景、D-mask两端和N个峰端点控制，故

`Lip=[(1−p)+(p/K)N/δ+(p/K)2Dsup]/L`，

原range_cap恰为Lip×h。窗口积分在class边界仍绝对连续，a.e.导数帽给整个bin的Lipschitz帽。N0只删峰项，背景和D项必须保留。

独立uniform X,X'的差平方恒等式给 `Var z(X)≤Lip²h²/12`，affine可达。**一般振幅范围不能直接除3；两点端点hidden law也不能使用1/12。** 当前counts提议只改变bin边际，原query在已分配各坐标bin后的隐藏位置仍独立uniform；当前来源的坐标对称性使同一counts可用任一canonical bin assignment代表。未把跨bin提议概率当作bin内密度。

## 2. 单核到Minkowski帽

同一个真实均值下包不变，每class relative cap恰为旧值/3，两个class取max仍/3。对非负t_i，展开正乘积有

`∏(1+t_i/3)−1=Σ_{非空D}3^{−|D|}∏_{i∈D}t_i≤(1/3)[∏(1+t_i)−1]`。

重复counts只是重复坐标，不影响此式。正系数多项式与分层Minkowski沿用原证明，表均值上包u_b保持不变，因而新variance帽 `Vnew≤Vold/3`。fixed0、joint、prior均满足；没有增加表数费用。不能拿 `(1+Rold)^(1/3)−1` 当上帽，其通常在反方向。

## 3. Lower变换及后处理边界

对g=(mean_lower−λ)_+及旧lower `b=g²/(Vold+g²)`，采用仍合法的Vold/3给

`b_new=3b/(1+2b)`。

g0约定b0；Vold0、g>0时b1，两端也正确。φ(b)=3b/(1+2b)在[0,1]单调，所以 `max_L φ(b_L)=φ(max_L b_L)`。这一步保留同一个共同尺度；不需要每L细节即可转换保存的conditional lower最大值。

payload存的是x=λWB，必须先利用logD恢复B，再用λWφ(B)构造新贡献；**不能直接φ(x)，也不能φ(旧CI端点)**。下溢丢失的信息不能恢复，保守保留0合法。数学上新lower≤真实条件tail≤旧upper，故旧Markov/source contribution帽仍可用。

旧upper含union、Markov截断和不同L的分解，不能盲目用同一个φ改写。未来要收紧upper须另保留/重算局部结构。本次没有应用任何payload变换；41轮应对新贡献重新冻结统计口径并重算对应CI，明确这是同一数据的确定后处理而非新独立确认。原40轮504项CI及结论保持原样。

全部解析改进仍以正确的均值/variance包为前提；原source三角函数、BLAS、erfc及经验FP guard的完整向外认证限制未因除3而消失。
