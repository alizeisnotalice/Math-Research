# audit_multibasis060：056多符号基81点闭盒资格独审

通过。正16-body Cauchy混合的全空间log-Hessian下界，可对同一H分别使用16个合法符号共轭基；逐点最大后得到分象限凹的下包围函数。anchor在闭盒内时，每坐标“下offset、0、上offset”的至多81点检查足以控制整个闭盒，无须全局凹性。

独立整数回读056固定的162个原split父盒，全部通过。2592个旧角点、13122个新分片顶点、209952个signed-basis整数曲率均实际核验；每盒的新最小margin及worst offset与056逐条准确相等。最小margin为
\[
\frac{33983968010587124693717903513538128433237870479571267959365256633789859952279}
{115792089237316195423570985008687907853269984665640564039457584007913129639936}
>\frac{293491}{1000000}>0.
\]
这给指定162个闭父盒的区域资格及批边界折叠的数学依据，不是全域g、全局覆盖或C5改善。本审不修改052树、控制器、队列、STOP或缓存。

## 1. 锁输入、有限运行与独立性

056原件在`Math-Research-cube5-astra-I/cube-dim-5-bounds/upper_multibasis056/`，manifest SHA
`0f3261b09dabc9557903547c0670125392e2bf106a8f4e7c3d952005b3978ba4`。其6件白名单逐hash/字节核对，其中report为`7aee0b3a6f685f5e722f2abd95cbf303357ef7ef4d4dfe8659c0b34e487c5303`，results为`564d454c0d0f2631843c1fcdff12f1b34596ba4697c177d1c6a2a6d69effa395`。

准确只读batch0002状态锁`b79c134fdcec283c120d464cf78f131c669c905ac5cdacda309a34f07a9fb29b`及登记`585a5de78b8f159ca03f138013894041dba350d144c3217b7bbda09fd5c7c9e9`；模型epsilon=1/64、target=1/10000、x1=1，四个物理tail梯度。读审043/047及046端点实现；全部输入hash见本result。

句柄`audit_multibasis060_integer_001`，实际先armed报告 **PID54933 / exec34340**，再解闸一次，约**0.524秒**完成，退出0。预算162盒/120秒；每16盒外置checkpoint `/tmp/audit_multibasis060_checkpoint.json`，未超时或重启。状态与登记计算后hash仍相同。未调用056或主checker函数。

本审只用标准库整数与Fraction，无SymPy、Arb或其它矩阵/区间库。新H、gradient、Green、kernel、log、exp、Gamma调用均为0。与作者共享的数值输入是冻结H/导数区间及Delta下端；没有独立重算其超越值。Delta来自046对准确H下端点/target作Arb log再取下端；读审该方向正确，但本轮不把源代码读审等同于独立log重算。基础模型归一化及缓存核资格沿043/046/047/052既有输入。本轮新的独立部分是理论推导、16基整数身份、全部有向比值端点交叉乘法，以及不同整数路径的162盒复核。

## 2. 从正混合到全空间抛物面

选同一16条body方向中的任一可逆五列基N，余11条时间1 Cauchy和为Y；无需Y有密度或有限矩。写
\[
q_t(z)=|\det N|^{-1}\prod_{j=1}^5\frac{t}{\pi(t^2+(N^{-1}z)_j^2)},\quad
H_\epsilon(x)=\int_\epsilon^\infty E q_t(x-tY)dt.
\]
对空间方向u，设a=N^-1u、w=N^-1(x-tY)，则score及log二阶方向导数为
\[
\sigma=\sum_j\frac{-2a_jw_j}{t^2+w_j^2},\qquad
\tau=\sum_j\frac{2a_j^2(w_j^2-t^2)}{(t^2+w_j^2)^2}.
\]
有\(|\sigma|\le\|a\|_1/t\)，且\(\tau\ge-2\|a\|_2^2/t^2\)。后者可由r≥0时
\(2(r-1)/(1+r)^2+2=2r(r+3)/(1+r)^2\ge0\)直接核对。q及其前两次空间导数由与x,Y无关的常数倍\(t^{-5},t^{-6},t^{-7}\)控制；epsilon>0故时间可积。支配差商/导数给H属于全空间C²，且0<H<∞。

将正可积权q/H归一为\([\epsilon,\infty)\times\Omega\)上的概率测度ν，score有界，全部期望有限，得到
\[
u^T\nabla^2\log H(x)u=E_\nu\tau+\operatorname{Var}_\nu\sigma
\ge-2\|N^{-1}u\|_2^2/\epsilon^2.
\]
这是有界score的方差，不是重尾Y的方差。沿固定直线x+ru积分，二阶余项权重\(\int_0^1(1-r)dr=1/2\)，故
\[
\log H(x+u)\ge\log H(x)+\nabla\log H(x)\cdot u
-\|N^{-1}u\|_2^2/\epsilon^2. \tag{1}
\]
公式任意有限位移、所有空间位置均成立，不依赖投影间隔、奇异集距离或线段内节点互异。每个基先固定并完成整条线积分，之后才能比较多个下界；没有对分段最大函数求Hessian。

## 3. 合法符号基、共同物理gradient及最佳曲率

固定N的列为043所列五条body方向。令M=6N^-1，其整数行是
\[
\begin{pmatrix}2&1&1&1&1\\1&-1&-1&-1&2\\1&-1&-1&2&-1\\
1&-1&2&-1&-1\\1&2&-1&-1&-1\end{pmatrix}.
\]
独立整数核N M=M N=6I、detN=-48。对s∈{−1,1}⁴，D_s=diag(1,s)把16条body方向置换，故N_s=D_sN仍是同一过程中的五个不同方向组成的合法基，|detN_s|=48、N_s^-1=N^-1D_s。重新选基仅改变正混合分解，不改变H、anchor或其**物理梯度**；gradient无需乘s。

本审用全部16个Ms=M D_s的tail列平方和，准确得到
\[
Q_s=(N_s^{-T}N_s^{-1})_{tail}=I_4/4-ss^T/36.
\]
整数Gram分子是9I−ssᵀ；s方向特征值5/36，正交补特征值1/4，严格正定。式(1)对u=(0,v)各给一个合法下界，取最大得到
\[
\log H(x+(0,v))\ge\log H(x)+g\cdot v-\psi(v)/\epsilon^2,
\quad\psi(v)=\min_s v^TQ_sv=\|v\|_2^2/4-\|v\|_1^2/36. \tag{2}
\]
三角不等式给\(|s\cdot v|\le\|v\|_1\)，s按v符号取值即可达到；零坐标任选符号。max下界的基可依赖最终v，不需要一组基沿位移路径变化。对所有v，\(5\|v\|_2^2/36\le\psi\le2\|v\|_2^2/9\)，且\(\psi\le v^TQ_{(1,1,1,1)}v\)。同号位移可能无改善，不能宣称每点严格获益。

## 4. 区间端点与81点覆盖全部闭边界

设真实anchor H∈[H_L,H_U]、物理导数d_j∈[d_L,d_U]、H_L>0。保存g_j∈[g_L,g_U]的包围资格，本审对四个端点组合准确核
\[
g_L h\le d\le g_U h\quad(h\in\{H_L,H_U\},d\in\{d_L,d_U\}).
\]
两边都是矩形上的仿射不等式，所以覆盖所有h,d，除以正h得到整个比值区间。这处理正负导数，不以midpoint或同侧端点比替代区间除法。对任意v的安全下界是
\(L(v)=\sum_j\min(g_{j,L}v_j,g_{j,U}v_j)\)：v_j>0取下端，v_j<0取上端，v_j=0两者均为0。Delta_L≤log(H_L/target)≤log(H/target)。因此
\[
\log(H(x+(0,v))/target)\ge\Phi(v):=\Delta_L+L(v)-\psi(v)/\epsilon^2. \tag{3}
\]
本轮直接比较准确有理Φ>0，无新的exp或H球乘法。若另转成显式H下界，H_L与exp下界均为正，应用下端乘下端；不能将不确定符号的gradient乘法按同侧端点处理。本轮线性项已按min正确处理。

−ψ并非全局凹：在(±1,5,0,0)取平均，中点的−ψ低于两端平均，准确差−1/18，本审准确核验。在每个**闭**象限s内，L固定为线性、ψ=vᵀQ_sv，所以Φ为凹二次式。设anchor在闭box内，offset为l_j≤0≤u_j。该box被16个闭正交矩形覆盖，每个坐标分别为[l_j,0]或[0,u_j]。凹函数在闭矩形上的值不小于顶点最小值，全部顶点的并集恰为
\[
\prod_{j=1}^4\{l_j,0,u_j\},\qquad\text{至多81点}.
\]
因此这些点Φ>0给整闭盒H>target。零坐标、坐标轴、边、面、相邻象限公共边界均包含；anchor在box边界或宽度为零时去重即可，仍覆盖全部闭胞元顶点。不要求拼接函数可微。本批162个anchor都严格在盒内，所以每盒恰81个不同点。anchor盒外不在本定理假设内；不能套本版本。

Hε的全空间正混合定理覆盖box内所有位置，包括零/重ray投影、排序边界及S1/S2。数值只有四个tail导数，故资格限于x1=1仿射平面，不能据此宣称五维数值盒。G≥Hε使相同区域下界适用于完整G，即使G在某些点无穷，也不需要遗漏时间段的上界。

## 5. 不同整数回读路径与结果

只读取指定state的162个splits，逐条核path为空、root对应物理盒\(\prod[r_j/32,(r_j+1)/32]\)、anchor匹配cache key且在盒内、epsilon/target/输入fingerprint一致。缓存H下端大于target，五个区间有序；梯度端点交叉乘法均通过。旧16点最小margin和最差bits与失败证书准确一致，均≤0；旧失败仅说明旧包围未过，不是H小于目标。

每盒把全部offset分母统一为V，再把Delta和gradient清到共同整数分母C。内循环完全整数：线性项使用清分母后的有向端点；对每个候选v，分别计算16个signed inverse的五行整数平方和，取最小曲率，再准确核等于\(9\sum v_j^2-(\sum|v_j|)^2\)。候选集由16个闭象限矩形的角点并集生成，再与三值坐标乘积集合相等交叉。该路径没有调用作者ψ函数或复用其Fraction角点计算。

162盒均新margin>0、new≥old；每盒minimum、worst offset、point count、improvement与056准确相等。最小新margin及逐盒结果保存在result。每个闭父盒都有\(H_{1/64}>1/10000\)，没有选择新anchor、重算球或读取后来batch扩范围。

## 6. 批边界转换资格与范围

**数学上可转换**：这些证书控制整个指定coarse-root闭盒，独立于之后如何分裂该root。若当前树仍使用相同物理root盒、epsilon=1/64和target=1/10000，在同步批边界可将该root当前全部后代叶/待办记录替换为已通过父盒。原二分子盒并恰为父盒，新父盒证书没有遗漏边界；因此覆盖账本的区域并集不会缩小，反而把原未证后代一起认证。缓存数值无需删除或重算。

操作时须按物理盒及node前缀核对身份，不能只按162这个计数；只折叠这些root的全部后代，保留其他根。须先确认无该后代的活动/inflight任务，或由主在已有同步协议中完成收据归并，再原子更新节点/队列/计数，避免并行worker回写破坏账本。本审没有检查当前活动树状态，也没有实施转换；这些是将数学资格接入当前控制器的边界条件，不是额外数学假设或新的授权请求。

这162个root不是全部目标域。未被其覆盖的区域仍需原有证书或新有限工作；本审不将部分批资格升级为全域g或C5<32。

## 7. 精简文件、hash与复现

仅4件白名单：registration.md、check.py、result.json（含162条准确结果）、本report.md。没有复制52 MB旧state、生产核、缓存或日志；无提交/推送，无MATLAB，052及旧锁只读。

登记SHA `b71d3fed5e8bb0169e26aa10bd5173c715d60d21616e3b39708340af182f5c55`；代码 `0a701e6993afc742ca1e6663d566932898570da9e829bdd90ba1841a5d7d2638`；结果 `f02927fe338e213e5359f925ceb5733b9d39e2219b492f546c799209d8052f1d`。

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/audit_multibasis060/check.py \
  --author-root /path/to/Math-Research-cube5-astra-I \
  --upper-root /path/to/Math-Research-cube5-upper \
  --state-root /path/to/frozen/upper_cover052/batch0002 --run
```

两个project root也可直接指其cube-dim-5-bounds子目录，默认无--run只核源hash。--arm先打印PID再等换行；上列普通复现无需该闸。源manifest/layout及锁定state必须保留；本审不在别处寻找更新替代输入。
