# R039 任意pair优化与三cube驻停例独立审计

**ACCEPT** `upper/round039/joint_cover/PROOF.md`，实际SHA-256为`b0719f3c04f0fb88fe794763e36e916636be58a0572db139db67ffb6e979c7d1`。未发现需要修正的数学或scope错误。

结论严格限于有限原witness支付证书：固定原L1输入及λ>0、有限正体积开盒、固定正预算0<h_i≤M_i/λ、非负q及Σh_iq_i=|U|。pair可激活原零系数，因此完整极大团列表必须来自所有allowed-positive原witness。单独的mass comparison仍使用h_i=|V_i|。本轮证明不提供新C_n界、算子下界、全局pair可达性否定或迭代有限终止。

## 1. 四affine分支与≤8候选

完整极大团按包含i,j的四种情况划分。外部q固定后，每类max的外部支付是常数；代入q'_i=x、q'_j=(P−h_ix)/h_j恰得到正文四条affine支，没有遗漏11类或把空类当成真实团。因h_i,h_j>0，合法域是闭区间[0,P/h_i]，P≥0。

最多四条直线有至多六个不同斜率pair交点，加入两端共至多八个候选。max-affine是连续凸分段线性函数，非恒定段内部不能产生未被端点/交点捕获的独立极小值；常值极小段也包含域端点或breakpoint。因此此候选表严格足够。平行/重合直线不需新交点；P=0时仅一个状态。

取得m后，各支≤m的交集就是完整最优区间。负斜率除法时反向得到x≥(m−b)/a，正文符号正确。该区间非空，因为m已有合法候选取得。

## 2. lex(L,Σq,E)精确tie规则

pair系数和为P/h_j+(1−h_i/h_j)x：h_i>h_j时严格递减，故取最优load区间右端；h_i<h_j时取左端。只有h_i=h_j时S在全区间恒定，需要再最小化E。

E_pair=h_ix²+(P−h_ix)²/h_j，导数为(2h_i/h_j)[(h_i+h_j)x−P]，二阶导数严格为正，故裁剪P/(h_i+h_j)至最优区间即为所需解。即使一般h不等时E更大，S的严格下降也优先，正文没有混淆旧lex(L,E)与新merit。

全局比较(newL,newS,newE,i,j)或同一旧状态下的S/E差等价。最优merit不增，接受严格下降只保证单步改进，不自动给有限终止或全局最优。正文已明确这些限制。

## 3. active峰团的必要donor条件

非平凡的合法pair守恒在正h下必有一个系数降、另一个升。若某个旧active peak团不含donor，它的支付至少保持旧峰：没有receiver时不变，有receiver时上升。这对任意有限步幅成立。因此donor属于所有active peak团是立即严格降L的必要条件，交为空时所有pair都无法立即严格降L。

这只是必要条件，不是充分条件；也不禁止plateau pair。若另有合法向量L更低，沿两向量线段由凸性给任意足够小正步长的严格下降方向；交空时该方向不可能只动两个坐标，所以至少动三个。该结论仍不排除把多坐标改变拆成经过plateau的pair路径。

权重Cauchy给E≥(Σh_iq_i)²/Σh_i，等号当且仅当所有q相等。故constant-q在固定支付下是唯一E最小点；这准确解释旧E tie在plateau上可能阻断路线，而不说明没有路线。

## 4. P3真实选择、最优mass与25%gap

对任意n≥1，中心0、±5e_1、原半径1的三个Q两两严格分离。witness半径4时中心与每叶差5<8而双叶差10>8，交图为P3，所有边均为正体积交。额外n−1坐标为共同区间，体积依次v=2^n、h=8^n，p=v/h=4^(−n)。

τ=3p使每个候选至多受到两个旧p项，eligibility≤2p<τ；任意固定正b均可full-select。因Q分离，每个increment为v，无rejected-flow项，q=(p,p,p)是真实原流程的基线。

mass primal在两个极大边cell各放h，总量2h，双叶各得h、中心得2h；dual双叶z=1、中心0，每极大边dual和1，目标2h。原始和对偶显式等值证明M(V)=2h。支付最优3p/2也可直接看：任一合法q满足Σq=3p，两个边load平均为(3p+q_0)/2≥3p/2；q_0=0、双叶各3p/2达界。

uniform q时，中心-叶pair保留对应边支付2p，叶-叶pair保留双叶和2p所以较大叶至少p；故任何pair都不能严格降L。全部h相等使S固定；uniform q又是唯一E最小点，故新三层merit也无严格pair下降。初始L=2p与最优3p/2之差p/2，确为初始25%。这是finite certificate gap，不是算子下界。

两步路径的精确merit为

* 初始：(2p,3p,3hp²)。
* 中心向一叶付p/2后：(2p,3p,(7/2)hp²)。
* 中心向另一叶付余下p/2后：((3/2)p,3p,(9/2)hp²)。

第一步E增而L不变，会被该严格merit规则拒绝，但作为合法支付路径存在。故驻停与“所有pair序列不可能改进”截然不同。f=(1+ε)1_(∪V_i)、λ=1使三个witness同时真实。

## 5. 实际最小核验与缓存观察范围

本机Python Fraction实际核验新增P3的严格交图、三个pair的最优merit、两步路径、显式mass primal/dual；另检查h_i大/小两种S优先而E上升的plateau、四类全在场、仅11类、P=0。只涉及新例子的有限代数，没有生产pair扫描或旧LP求解。

另只读三个已冻结缓存的q_after和已缓存clique列表，实际SHA及本轮peak观察均与原文一致：

| 缓存 | 已缓存极大团 | 当前峰团 | 峰团交 | q恒定 |
|---|---:|---:|---|---|
| C/2/5 full_disjoint | 87 | 3 | 空 | 是 |
| C/1/6 full_disjoint | 71 | 24 | 空 | 否 |
| A/3/15 partial | 37 | 24 | 空 | 否 |

这仅核验原文§6新增只读观察及源hash；缓存团的完整性继承既有核验，本审计未重新枚举几何图/团、原union、selection或mass LP。三缓存完整SHA与实际核验结果保存在`verification.json`。峰交空只决定这些当前状态无立即pair严格降L，不决定是否有S/E plateau改进，更不提供注册batch统计。

复现本审计最小检查（第二路径参数可省略，则不读缓存）：

```sh
PYTHONDONTWRITEBYTECODE=1 python3 cube-dim-7-bounds/audit/round039/joint_cover/verify.py cube-dim-7-bounds/upper/round039/joint_cover/PROOF.md cube-dim-7-bounds/upper/round038/compute
```

本审计只写本目录报告、复核脚本和精简receipt，未做Git修改、n=7工作、MATLAB、生产LP、79-case复跑或新uniform C_n结论。
