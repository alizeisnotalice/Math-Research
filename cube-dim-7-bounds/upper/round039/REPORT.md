# 第039轮：共同质量支付与任意配对的准确接口

本轮没有新的统一Cn数值界；当前已归档上界仍为R036，R042更强算子链已审并将另轮封存。以下是对原有限输入/见证的严格接口与反例，不是Cn下界。固定Lebesgue、centered轴平行cube、共同半径、一般L1及全部正尺度。

## 对任意共同f支付一个star簇

设同一真实history的若干witness S_i=aQi有共同交点，R为最大半径、j为最大项。原真实increments互不交且I_i⊂Qi，而所有Qi的并集落在以共同点为中心、半径(a+1)R的cube中。因此λΣv_i≤[(a+1)/a]^n∫_{S_j}|f|。它对任何满足这些原见证的共同f成立，远强于只给一个spike例；原pointwise W仍可能约nlogn，说明W热点不是必须付到整个输入质量的成本。

对多个簇J，取θ_J≥0覆盖每个selected index，且其最大半径代表witness的点态重叠Σθ_J1_{S_rep(J)}≤K，则λΣselected v_i≤c_a K||f||，c_a=(1+1/a)^n。每个簇使用原history的increments，不重新计算子history。这个代表重叠是必须支付的成本，不能免费加多个star的近1常数。若对全部原histories统一构造K再接原coverfactor，才会产生Cn改进；本轮没有该uniform K定理。

有限完整LP不能只保maximal cliques并随意换代表：用全部positive-volume arrangement cells C及其中每个代表j，列J_{C,j}={i:C⊂S_i且ri≤rj}。任意真实common-star簇可在共同open交集避开faceplanes后扩到这样的列，保留代表j不变，且只改善coverage。因此列族无漏、方向正确。完整证明joint_mass_dual/PROOF.md及独审audit/round039/joint_mass_dual均ACCEPT。

原centres作固定pivot的最大star greedy可能失败：R_i=3^-i,y_i=.75R_i e1,r_i=R_i/a，适当a满足a>max(3,2t)、(N−1)a^-n≤τ。smallQi不交且全eligible，stars共同含0；逐次以当前最大star自己的centre作pivot却每次只移走自己，代表重叠达N，改用自由pivot0则一簇K=1。任意N需要相应a，不能说fixed a/τ无限N。这是该greedy的反例，不是Cn障碍。

## 任意pair的完整max-affine及三cube反例

保持原支付Σh_iq_i=U（预算正、共同f固定），选任意i,j。完整cliques按含两点的00/10/01/11四类分组；各outside最高价给至多四条关于x=q_i'的affine lines，q_j'=(P−hi x)/hj。区间端点和相异斜率交点最多8个有理候选，准确求最小load m；再交各line≤m得最优区间。在此区间若hi>hj取x上端、hi<hj取下端最小Σq；相等再clip P/(hi+hj)最小E=Σh q²。每次lex(L,Σq,E)严格下降才接收，不是原E-only势。

三个原cube半径1，a=4、中心0及±5e1（任意n≥1），tau=3p,p=4^-n，fullselected合法且increments真不交。witness图P3，h=8^n，初q全p，真实L=2p。两edge各mass h primal与两leaf dual=1给M=2h，最优有限load=3p/2，差距占原load25%。所有pair均不能立即strict降L；等h付款使Σq恒定，uniform q又全局最小E，所以连新三层merit都pair驻停。但三点方向或允许E先上升的两步plateau pair路径能到最优。它反驳局部势的万能性，不是Cn下界，也不排除其他pair路径。joint_cover/PROOF.md及独立audit/round039/joint_cover为ACCEPT。

## 六个新单步扫描：实际运行与准确证书

Terra仅对R038三个已审LP失败终态各做nonedge/all-pair一轮单步扫描，总6个新任务；没有重跑79例、原atoms、旧LP、原输入搜索。完整算法/6输入hash/停止条件在registration.json启动前冻结。六个chosen pair均L不变、Σq严格下降、E增加：C/S2/j5两域(0,7)，C/S1/j6两域(24,55)，A/S3/j15非边(66,79)、全pair(3,11)相邻。这说明旧E-only拒绝的plateau可由新势打开；没有证明以后迭代必降load或全局最优。

所有分支、交点、opt interval与chosen payment均为Fraction，原完整有理证据在compute/scan*.json。独立审计复用R038已审input/cliques hash，重新推导chosen pair全部四类和端点/clip，验证付款与真实cost、严格Σq比较全部PASS；未重扫所有pair，所以独审只认证chosen pair，不认证global-best选对。登记cc2424...及proof b0719...身份吻合。实际6workers、wall约12.42s，队列快项自然结束，不空转；资源及PID见compute/summary。没有n7或下界工作。

最小复核audit/round039/joint_cover/verify.py为新P3/Fraction代数，audit/round039/six_scan/audit_six.py核6个chosen证书（保留原主机输入路径）。全参数定理由解析证明而非有限样本承担。公开比较沿已核原始来源，本轮为内部接口，不称发表、全球最好或新的Cn记录。来源hash见MANIFEST。goal继续active，新selector/共同mass和R042正向算子改进独立推进，旧论文修订另轮交付。
