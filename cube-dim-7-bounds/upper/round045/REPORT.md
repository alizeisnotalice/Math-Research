# 第045轮：逆列加权重算行预算，得到负线性修正

本轮新的严格内部原算子结论：存在绝对整数N，使每个n≥N满足

    Cn ≤ F(A_n(t*), B_n(t*)−2^-40 n)
       ≤ nlogn−nloglogn−2^-41 n，
    d=2^-24，L=logn，x=(1−32d²)L，t*=n/x。

一般Lebesgue、中心轴平行立方体、全部坐标共同半径、一般L1及全部正尺度均保留。N尚未量化，L≥2^240仅是下述完整cost比较的标量门槛，绝不是已认证operator门槛。R042在其原显式域L≥2^108继续有效，其他旧界也保留各自域取min；当前有限维七维界没有用本结论替换。首项系数仍1，不声称全球最优、公开发表或Cn的真实渐近等式。

## 新机制与准确参数

已完整审过的主机制对每个固定0<d≤2^-24成立：h=2^-20/(Ld)，theta=exp(−12Ld²)，k=nh theta，beta=1+1/t，zeta=n^-3，kappa=beta−zeta。目标半径r=1+u、来源半径rho=tq=t(1−v)，u,v在[d/4,3d/4]，匹配|u+v−d|<h，捕获||z−w||∞<kappa rho。环境lambda是原角点未捕获cap D上的2^-n m(z)^-n dz；actual mu仍为原算法同一history的真实增量，mu≤lambda。

对每个source q设R=q+d−h，a0=(kappa−1)rho，g=R−a0=d−h+zeta rho。所有匹配targets有r∈[R,R+2h]，归一交长ell_R(z_i)=min(1,(R+kappa rho−|z_i|)/(2R))，G(z)=∏ell_R。新来源权重

    b(z)=min(1,theta/G(z))，0<b≤1。

逐点旧几何column≤(4nh/R)G；因此新column≤4k，不需要统一高列cap H来除旧energy。关键是同时重新推导row，不声称保留旧unweighted row的固定比例。b≤1使同一实际target j的真实eligibility仍只支付tau，不多一个theta^-1。

## 全匹配半径的inverse row

定义J(w;q)=E_source[1_capture/ell_R]。精确边带原函数给闭式(见主源§4及独立LEMMA)；在参考target radius r=R时E_W J=1严格成立，ElogJ在[−6g³/rho²,0]、VarlogJ≤12g³/rho²，最大面pin坐标另保留其O(g/rho)贡献。不能将source边界熵Ld³误当作这一新row的代价。

以每个actual target r的q0=r−d+h作baseline，n−1自由坐标的log和均值及方差是O_d(L²/n)，趋0；外边带数Kw=1+Bin(n−1,2d/r)≤4nd的坏概率也趋0，原target cap补集指数小。故可取同一个Borel目标好集A，lambda(A)≥nd/16。精确J关于q局部Lipschitz，其a.e.导数满足

    |d_q logJ|≤64/rho·[d²+1{|w|>r−2d}]。

在kink处使用绝对连续性积分，不声称处处可微。对同一w求和且遍历全部匹配q，变化≤1280Ldh=5/4096<1/8，从而所有匹配rho统一满足1/2≤∏J≤2，完全不需有限q网格外推。

完整source face系数为U=(∏J)/(2ell_pin) Σ_i N_i/J_i，N_i=1+1{|w_i|<a0}；保留全部2n source faces和target固定坐标后，n/8≤U≤8n。按capture/G归一的inverse-row概率下，每个free坐标的E(−logell)≤d²/rho，pin≤d，总和≤3Ld²。Markov因此给inverse-row的G<theta坏概率≤1/4；这不是原始source坏概率。故soft row系数在[3n theta/32,8n theta]。

## 同history缺额与原算子传递

匹配径向log长度在[h,4h]，完整soft rows在[3k/32,32k]。来源原cap补集的完整质量O(nd exp(−Ld/32))，由b≤1按原质量扣除；相对k为(d/h)exp[−Ld/32+12Ld²]→0，域内d/32−12d²>0。故最终ambient row下界3k/64，上界32k，column≤4k。

选tau≤3k/128，同一缺额D=lambda−mu满足Eλλ−Eμμ≤32k[D(A)+D(P)]≤32k·Delta，且Eμμ≤tau lambda(A)。两个径向集合不交，不能把缺额双算，得到Delta≥3nd/65536；保守取delta=nd/65536。对所有a_alg≥8n^4，角点环境在A/P上精确同乘ca=alpha^-n；energy乘ca²、degrees和目标质量乘ca。实际target item半径<3、source>4且kappa+2/(a_alg+1)<alpha beta，使每个加权捕获都进入同一真实target的一次资格测试。新b≤1保留actual row≤tau，缩放不削弱缺额。

取tau=sqrt[A_n(t)(B_n(t)−delta)]，x=(1−32d²)L使A_n≤n exp(−32d²L)、B_n≤nL，tau/k≤sqrtL/h·exp(−4d²L)→0。把uniform finite-a缺额插回原finite cover和真实witness积分，仅在标量上令a_alg→∞，得到Cn≤F(A,B−delta)。任意L1及全正尺度由真实witness有限覆盖任意紧严格水平集、Lebesgue内正则性传递；没有support假设、尺度截断、选择器极限或不同history拼接。该覆盖路线不使用Green核、小时间积分或n≥5奇异子空间，不把四维线外连续证明迁移至高维。

## 端点与完整新旧cost

完整cost的正损失已付清：B_n(t)≤T+64nd²+1，A+2sqrt(AB)=o(n)，T=nL−nlogL。参数族全链ACCEPT后，独立参数审核确认d=2^-24属于原闭域，全部d-dependent geometry、good-set、inverse tail、cap、eligibility及finite-a guards都能共同成立；并非只替换展示系数。端点delta=2^-40n，64nd²=2^-42n=delta/4，剩余o(n)可付delta/4，故得到开头的T−2^-41n。作者原展示d=2^-30只给T−2^-47n，亦保留完整准确预算作为较弱实例。

另一完整独立版本layered_source选d=2^-64、K=2^20、B=4K、H_R=2^24，用omega=G^-1及Phi≤Kxd²截断；normalized column≤4nh，highPhi按原source质量付能量，真实row正确付Q tau。它给eventual T−2^-97n，两条全链各有独审。它的跨半径完整pin-face二阶矩解析证明不依赖R046/R047数值。

三份完整cost的准确参数为：

|实例|d|sigma（x=(1−sigma)L）|实际扣减Delta*n|展示c（T−c*n）|
|---|---|---|---|---|
|I端点|2^-24|2^-43|2^-40 n|2^-41|
|I作者展示|2^-30|2^-55|2^-46 n|2^-47|
|S独立方案|2^-64|2^-106|2^-96 n|2^-97|

比较附录使用完整Phi=F(A(t_sigma),B(t_sigma)−Delta*n)。原cap给Phi_old≥T−n(4/L+Delta_old)，并对R042给Phi42≥T−n[4/L+2^-31 L^-1/3]；旧fixed-t家族共同lower为T−n[4/L+2sqrt(logL/L)]。Taylor全尾和准确Fraction在L≥2^240证明：Phi_I端点<Phi_I作者<Phi_S<列明且自身适用的每个旧完整cost。old unit-core跨全部t族仅eventual比较。12个旧q读取原准确certificate，最大74275164113/10^12，未从小数重建。

这个比较不能变成Cn下界，也不能把新unknown N设为exp(2^240)，更不能在尚未认证的新域删除旧有效界。主链与比较分别审过后，conditional附录中的Cn激活条件现已满足“主定理接受”这一项；每个参数自己的有效起点仍保留。

## 审计、实际工具及公开范围

geometry全链、one-coordinate/AC子链、另一layered-source全链，以及完整cost和端点guard均有不同Sol独立审阅。主代理也重建inverse law、clip支付、同history mu bridge及完整scalar，保留最终证明责任。审前原源中的candidate文字冻结保留作为实际溯源，当前接受身份由本报告和最终audit负责；未无声修改原hash。

本轮实际工具是作者及独立审阅者各自必要的Fraction/SymPy/Arb小核对，保留原receipt与hash，没有重跑其通过基线。R046的4365高精度积分和R047的61440跨半径组合属于另外登记的新有限诊断，不替代本轮全参数证明；R044的支付LP也不自动成为Cn改进。最小代码在各source/audit目录，部分原脚本按当时独立clone绝对路径锁源；跨机器复核应在临时副本适配来源根目录，不能把推荐复现记作新运行。

公开比较仍限已读原件：Stein–Strömberg 1983与Aldaz 2017 v2的固定epsilon界等，内部新结果不称公开发表。另作两条有界最新检索并核官方摘要，[2026年Spector–Stockdale预印本](https://arxiv.org/abs/2609.05377)的O(sqrt(n)logn)明确针对Euclidean balls，不迁移为cube界；[Aldaz v2元数据](https://arxiv.org/abs/1605.05596v2)与原归档一致。新增HTML实际下载hash及读取范围在literature_refresh.json，未声称穷尽全部最新文献或新读两篇全文。

当前paper/main.pdf仍是R041已终审实编译24页稿，证明R042；R048正在把本轮新主定理完整写入并将另行全文审稿及编译，未提前声称当前PDF已包含本轮。继续一般高维覆盖；下界与七维专项不恢复，原资料包只读，其他项目未修改。R049下一步尝试显式量化本轮真实门槛，尚未接受，绝不提前拿候选基点作operator域。
