# 第025轮：第一轴period3位置位移的首两方向

固定精确C格 `a1=1+u^3/2`、四tail间距 `s=(1+u)/2`、`u^5=2,v^5=3`。单位原子质量不变，仅令首轴位置

```
p_j(epsilon)=j*a1+epsilon*z_(j mod3).
```

首有限批准确计算 `z=(1,0,−1)` 与 `−z` 两个一侧方向。两份导数的精确Q(u,v)系数相等，且严格负：

**D_plus = D_minus ≈ −0.11835180380160859195153514754135355395291228452062**。

这两个方向的range为2；若归一range1，对应值应再除2。本轮没有正方向，没有非零epsilon候选，没有执行完整位置fan；结论只排除这两个方向在该1024长度有限functional上的局部一阶改进，不是Cn或一般位置输入的上界。

## 查重、不同机制与真实对称

两个root先只读检查period3 position/displacement历史，未找到本族的已运行排除。相关旧009是**period2**交替首轴位置，delta=1/32、1/8、1/4三个点；Astra001/002是插入一个额外小质量卫星轨道。它们不涵盖本period3单位质量原子位移。旧024是质量扰动，未重跑其union或已通过C基线。

位置族的首轴平移诱导相位循环；反射x→−x并重标j→−j则诱导 `z'_j=−z_(−j)`。因此真实群为 `G={P_even,−P_odd}`，不能直接复用质量族的纯S3排序基本室。`(1,0,−1)` 与其负向量在G下各有3点orbit，两个orbit不相交，所以本轮分别计算；相同导数来自实际准确结果，而非错误的对称假设。

下一fan的资格可改用两条半边

```
±((1-t)/2,t,−(1+t)/2), 0<=t<=1/3.
```

理由是range1、零均值方向构成六边形的边界；G把6边分成两个3边orbit。群中带负号的max/min交换把每条代表边的t映到−t，所以每代表边只需一半。这个资格已用准确有限群运算核对，但本轮没有提取、登记或计算该fan的事件。

## safe starts、互异原子、period mass与共同半径

对两方向取 `0≤epsilon≤1/8`。相邻首轴原子间距至少 `a1−1/4>0`，所以索引顺序严格保持，首末原子确为block位置极值；不同首轴索引给不同位置，配合正tail间距，Cartesian原子均互异。每原子质量恒为1。

`p_(j+3)=p_j+3a1`，每3a1周期恰有3个不同单位原子；五维周期质量恒为3，没有质量导数项。任意长度向量ell有 `N=Πell` 个不同原子，其唯一共同立方体半径为

```
r=N^(1/5)/2, 32*r^5=N,
```

与epsilon无关，五坐标同时使用这个r。没有因位置变化而独立调节其他轴半径。

完整 `ell∈{1,2,3,4}^5` 共1024项。N≤1024使r≤2；首轴位移幅值≤1/8，又2a1>3，因此j≤−2时接收区间上端 `p_j+r<0`，j≥5时下端 `p_(j+ell1−1)−r>3a1`。所以安全start范围为 **j=−1,...,4**，每相位两个平移像、每方向完整6144 raw候选。没有只套旧42盒。

尾四轴不随epsilon变，质量也不依赖tail相位，因此canonical奇偶tail block仍支配其他tail starts，四tail反射仍合法。接收域为 `[0,3a1]×[0,s/2]^4`；物理周期比准确为 **16V/3**，不除格胞体积、不用旧32反射。

## 新静态四维section与一次准确到限接续

先从全部1024长度构造固定tail四维盒。213项有正tail体积，其余811项固定为空或零宽，不会因首轴位移而出生。初始tail cut每轴26段，共456976个cell。登记每进程最多200000新cell、无wallcap，首进程PID74134正常到node上限：10.4052秒、峰值532185088字节，完整保存11/26 first-tail slabs和193336个cell。没有fault或硬kill。

这次真实尺寸暴露了可证明的局部冗余。未裁切首轴接收宽度是

```
2r−(ell1−1)*a1 + epsilon*(z_first−z_last).
```

213项中，162项的常数宽度严格负。因为 `|z_first−z_last|≤2`，可用全部严格负宽度的正间隙下界给出一个统一正的epsilon上限，使这些项仍不可能产生首轴receiver。其准确margin保存在 `projection_certificate.json:valid_epsilon_upper`。所有常数零宽或正宽项均保留，共51项；未用严格正宽测试排除潜在出生盒。

随后先登记最小修复：将已保存fine histogram的mask投影到这51项，体积系数准确聚合，保留其整个 `[0,x_cut]×[0,s/2]^3` 成功区域；只对未完成的 `[x_cut,s/2]×[0,s/2]^3` 使用51项的粗cut积分。**193336个成功旧cell零重算**。PID74488用3.0557秒、峰值258359296字节补算36864个新coarse cell，得到2224个完整静态section。

mask投影只是删去已证明在统一小epsilon窗口不可能接收的长度，不改变该窗口的finitefunctional。全1024 raw身份、213正tail盒及162个准确跨度排除均保存。另一个独立scalar核验从1024长度重新生成8192个tail端点保护式，确认213/51/162分类，见 `inventory_verification.json`。

## 只做一维union的精确一侧导数

对每个固定tail section，其mask保留共同长度向量身份。相应首轴区间是

```
[max(0,(j+ell1−1)*a1−r+epsilon*z_last),
 min(3a1,j*a1+r+epsilon*z_first)].
```

全部端点是准确仿射函数。使用jet `(Q(u,v)常数位置, 有理速度)`，按一侧epsilon>0先比位置、同位置再比速度，完成max/min裁切、空集判定和并集。真实首轴位置顺序已有正间距保护；tail sections固定，所以每section并长度及全finitefunctional在epsilon中都是有限分段仿射函数。

并脸同值异速会产生一阶窄条，必须保留；同值同速高阶差在本位置模型中其实没有（端点准确仿射），也不会制造虚假事件。完整raw检查中，两方向各179个可行raw、6组同位置异速度tie，没有新增常数零宽出生盒。空tail或固定零宽tail始终零体积，不是假设其影响消失。

每个section分别用一维区间merge和独立cut-coverage mask积分求首轴并长度导数，两引擎在全部2224 section准确相等。结果每方向均为：2064个section速度0，160个section速度−4。因此

```
D = (16/3)*Σ_section weight*speed
  = −(64/3)*Σ_(160 affected sections) weight < 0.
```

这些weight均是正的真实四维section体积。严格负既有该几何和式解释，也有完整Q(u,v)系数及有理根隔离上端<0的证书。

没有计算旧全C基线。epsilon=0的原C值由完整003 source和3周期/四tail反射身份继承；这里只评估新位置方向的导数。

## strict level与collar/L1接口

本轮r是等质量半径，闭receiver边界取值不影响Lebesgue体积。对任一固定允许的epsilon，从有限receiver并集中先删边界，再把所有block共同半径统一减小delta>0，取 `r_strict=r−delta>0`。此时 `32*r_strict^5<N`；各坐标仍使用同一r_strict，向内receiver并体积在delta↓0时趋向名义finitefunctional。不能将十个面独立移动以获得虚假严格余量。

截取M^5个完整源period，在各索引方向加固定collar（例如8个索引足够覆盖所有本批start、最长4块和有界位移）。输入质量为 `3M^5+O(M^4)`，接收体积为 `16V_strict*M^5`；M→∞时比趋向 `16V_strict/3`。不同period中的位移图案一致，collar代价只在边界，period mass没有改变。

再以同质量1、sup-norm支撑半径eta的非负光滑团替代每个互异原子，选eta小到 `r_strict+eta<r`。receiver用共同半径 `r_strict+eta` 可完全捕获相应N个bumps，而 `32*(r_strict+eta)^5<N`，故严格level1保持，L1质量不变。先M→∞、再delta↓0给出名义有限ratio的下界传递接口。主代理保留正式proof/collar终审；本轮两方向失败，没有以该接口宣告新Cn。

## 实际运行和最小复核

所有数学进程顺序运行，nice10、单线程、六个数值线程环境变量均为1，无程序wall时限。除第一次预登记node到限外没有fault、盲重启或无限搜索。

| 对象 | PID | 秒 | 峰值字节 |
|---|---:|---:|---:|
| 首fine static批，到有限node上限 | 74134 | 10.4052 | 532185088 |
| 准确投影及未完成区域接续 | 74488 | 3.0557 | 258359296 |
| plus一维双union | 74628 | 6.7993 | 39518208 |
| minus一维双union | 74907 | 6.6675 | 38633472 |
| 完整static独立pair-factorization审计 | 75438 | 1.1865 | 159858688 |
| 完整1024 scalar资格审计 | 75576 | 0.1313 | 23773184 |

最小审计将tail轴倒序，先构造两份2D mask histograms，再用27556个mask对重建**整个**2224-key四维histogram，与最终每个精确系数相同。这是一次明确登记的证书复核，不是新候选或5D/参数union。它还核对51个tail marginal、41784个逐block互异原子索引、两个6144 raw计数、完整section prefixes、同一准确D、signed-D3两个3点orbit不相交。旧023/024全白名单哈希未变。

该完整final certificate复核使31MB fine prefix不必归档。原fine checkpoint及其压缩副本仍保留本地供恢复；白名单只保留其source/hash/成功prefix回执，完整最终tail histogram及其小型独立验证器足以审核本轮数学结论。原始成功prefix没有删除或改写。

完整结果审核命令（均不运行新的参数或5D union）：

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round025/verify_saved.py
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round025/verify_inventory.py
```

source、登记、成功prefix、修复及审计对象的完整SHA在各registration及白名单中。本轮只写025，未提交、推送或使用MATLAB；不是倒填012，交主代理作为下一013检查点材料。

未解决：完整位置fan尚无事件登记或计算；两个中方向负导不能排除两代表半边的内部。finitefunctional之外的长度、其他输入和全尺度问题也未排除。
