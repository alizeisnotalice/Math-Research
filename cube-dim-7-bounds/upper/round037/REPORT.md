# 第037轮：真实点态负载障碍与独立拆分支付

本轮没有新的Cn上界。当前严格上界仍为R036：loglogn≥128时，Cn<T−n sqrt(loglogn/logn)/131072，T=nlogn−nloglogn。本轮的严格新结果说明当前点态负载框架的限制，并构造能够越过集中支付停滞的局部机制。原算子始终为Lebesgue、中心轴平行立方体、坐标共同半径、一般L1、全部正尺度；下界暂停，七维专项停止。

## 真实选择历史已达到leading1

固定n≥3、t>S≥2、τ>0，β=1+1/t。令ν为原corner未捕获包络限制在S<||z||<t的测度，M=ν(total)，M0=nlog(t/S)，Γ=(n−2)S/(16t)。完整入射势T(w)=∫1_{||z−w||<βm(z)}dν(z)满足∫T dν≤M0²e^-Γ，故删去T≥τ/4的点损失≤4M0²e^-Γ/τ。

这次结论恢复到实际有限原实例：对任意ε>0，存在a0(n,t,S,τ,ε)，使每个实a≥a0都有一组真实小立方体，全部原eligibility成立、半径非增、实际increments恰等于不交Qi，weights均a^-n；0在每个aQi中，最后r=1项是真正正权唯一最小项且未被任何前序bQi捕获。负载满足

W(0)≥B_n(t)−1−nlogS−4[nlog(t/S)]²e^-Γ/τ−ε。

完整证明general_cover/PROOF.md及独审../../audit/round037/real_history/proof_review.md均ACCEPT。有限恢复用good open set内紧盒、R_l略大于盒内m(z)、统一超平面边界层控制与网格积分上包；网格的all-source负载（包括未来及自身）小于τ，因而每个真实prefix更小。a0允许依赖上述参数，不冒称整个维数窗口共用一个a0。

对固定C≥0、t∈[n/(logn)^C,n(logn)^C]、τ≥1，取S=(logn)^(C+3)，得真实历史W(0)≥nlogn−O_C(nloglogn)。因此仅提高该原算法在这个窗口内的统一点态W估计，不能降低首项系数。这是**框架负载障碍，不是Cn下界**；不排除联合使用覆盖/f见证质量、其他选择机制或窗口外参数。

general_cover/CLARIFICATIONS.md明确I_cap是原帽壳积分，非R035的二点容量。同一高负载history有共同star见证；小support尖峰可同时strict witness且固定n,Mf,λ下的水平集比值趋1。这是固定level观察，不声称整体weak norm在homothety缩小下趋1；整体norm尺度不变。

## 集中支付的精确失败与新的拆分机制

R035集中规则把clique C支付搬到最大预算anchor，必须有s_C=Σ_Cq−Σ_Chq/max_Ch>0才接受。若所有正q都在全局最大h上，则对任意clique及任意子clique都有s_C=0，与枚举pool大小无关。

最小解析反例可取任意n≥2的四个r=1小cube，a=4、中心±4e1与±4e2。small cubes不交，witness图为K_2,2，h=8^n、p=4^-n。b=7、τ=3p/2、tie序0,2,3,1给真实选择0,2,3并拒1。拒绝cube两半分别被E-memberships {0,2,3}和{2,3}覆盖，比例flow准确给q=(7p/6,0,17p/12,17p/12)，L0=31p/12。集中法全停，而mass LP的四edge各h/2 primal与两个非邻顶点dual=1给M=2h，最优支付load=2p；严格gap为7p/12，即原load的7/31。

将q0均分给非相交witness 0和1，新的每条edge负载均2p，恰达到该有限LP最优。这个四cube例不是注册N=48/64/80/96 batch中的案例，不能单凭它声称注册batch H失败；它反驳该半gap主张的全有限族推广。

一般独立module F可保持总支付Σ_F hq不变，把各q改成t=Σ_Fhq/Σ_Fh；真实load不增，严格下降条件见joint_cover/SPLITTING_PROOF.md。一般独立集需支付外部边界Δ_F：L'≤max(Lrest,L−(qmax−t)+Δ_F)。对非邻pair，四clique类缩为三条仿射线，精确交点clip给固定pair支付下的真实最优；这是新可证局部接口，尚非全几何统一Cn界。独审../../audit/round037/joint_cover/README.md为ACCEPT，正预算和支持图scope已收紧。

## 实际本机批次及注册差异

Terra实际15进程完成240几何：223完成，17按该实现的32次生成上限不可用；383支付模式原程序标202 pass、79 fail、102 gap不足或未分离。wall约35.13s，累计worker CPU约119.86s，实际采样7个worker近满，不能声称持续15核满载。没有MATLAB或n7任务。

生产实现与书面SPEC不完全一致：tie hash及并列选择、A direction slot顺序、C重复sign及retry规则有差异。原输出没有逐件存全部原对偶向量，且102类合并。**这些标签只是实际实现的未逐项审核诊断，不能升格为同一注册batch的严格频率结论。**原SPEC、IMPLEMENTATION_REVIEW及原sourcehash均保留；没有事后修改规则伪称预注册一致。

只有一个冻结实际48-cube A/S0/j6 partial实例另获完整独审：真实整数网格重建union/increments/两个拒绝类及比例flow、合法selected [0,2,3]、独立168节点Bron–Kerbosch完整40cliques、可行mass dual与显式改进支付均通过。所有正q在全局maxh，故任意集中clique全停；LP上端严格低于原load。其48cube改进支付未证全局最优，不与四cube的exact optimum混淆。证据见../../audit/round037/joint_counterexample/。

## 准确比较、复核与后续

内部Cn门槛不变；新旧比较只涉及已证明的实际支付成本。四cube全有理31p/12→2p与48cube的完整分数均来自原件，未从展示小数重造。点态W框架下界不是Cn下界。公开对照沿已核原始来源，未称新结果发表、全球最强或检索穷尽。

最小复核audit/round037/joint_cover/verify.py验证新解析四cube和clip；audit/round037/joint_counterexample/audit.py独立重建实际48cube，均标准库Fraction。原生产代码、最小实例、源hash及240条紧凑状态见compute/。15完整JSONL持有后续所需79个失败q缓存，仅留本机原独立clone，不删除、不重复归档；它们的hash已封存。历史绝对路径保留源身份，不称跨机器自动可携带；复核会写结果时须在可丢弃副本执行。解析proof承担全参数结论，有限样本不替代它。

R038新非边拆分直接消费79个已缓存q，不重新扫描旧atoms/选择/基线LP。有限family成本改善不自动代表Cn下降。R039/040继续全几何共同mass-dual、允许改变选择的机制；所有本机生产由Terra，Sol推导/独审，主保留证明责任。旧论文仍有效较弱，未重新编译。goal保持active，其他项目和下界状态不改。
