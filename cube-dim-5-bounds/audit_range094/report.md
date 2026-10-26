# audit_range094：091 新比值与 092 全阶段独立有理审计

**完整通过。** 对冻结 092 的 56,139 节点、4,098,147 个正时间段及全部 14 个完整前缀，094 独立确认新 range、mean、量化凸弦、旧证书交集和全域层级付款均合法。没有发现阻断该完整 I 上证的解析或实现 gap。

实际严格简上端为

\[
\boxed{I=\int_{[0,1]^4}G(1,y)^{5/4}\,dy\le\frac{6010267}{10^{11}}.}
\]

它严格小于旧完整 cap 70301737/10^12。内部准确上端的显示值是 UH≈5.669574350580058×10^-5、I_upper≈6.010266982568090×10^-5。已接受 g=1/10000 的配对 cost 上端约44.93722703599，准确32门为 false，**正式 C5 上界仍32，没有新 Cn 结论**。上证未过门不能推出实际 C5≥32。

## 1. 冻结范围与独立性

092 最终 manifest `05cd878504800f4a66673d9a408666c6480b92a1db5828de01f6ba5ec0cd801c`、report `951366667b8d5064b06348677565c1bf65d9d10340a71092fdc78b2c424bfa64` 的全部65件字节在运行前锁入 whole_scope；没有以动态 progress 作为最终证据。094只在 astra-cover 的本目录落盘。

091 manifest `6589e35086bd5082c66c4e94219b0f7c80c393e735f5011eb6b2a6d24d22665b`，ratio.py `2931e5d34b71d93eb11dcdb52d5ea0921f9fdbfa0e5e814a60f4ff64eed0de3d`，092 catalog `b33176d60b4c7942262169f14cf3bff7db40ad255527e5661831a7f9cc642e74`。另完整核088/082/080/085/076/075/087冻结 manifest 文件；准确来源清单见 whole_scope 和结果中的 manifest_file_counts。

本轮独立的是091解析推导审阅、89项新代数比值、全部节点新range/mean的准确有理消费、整数四次根凸弦和14前缀完整I合成。原始 H 区间的数值生产仍共享079/080/082既有 Arb证据；085/088已通过的全部M2矩证书明确复用。094没有把“共享端点输入的不同下游算术”称为端点数值独立重生产，也没有重跑旧087成功空间边界/投影实验或088全量moment基线。g沿用主round012完整接受，未新认证g。

## 2. 解析资格

独立逐步审阅见 analysis.md。关键步骤是固定合法 lift 在 K=ker B 上的测度保持平移，单Cauchy比值二次式判别式准确为零，五方向chain lift确实满足Ba=u，四个非负gap之和r给Jensen上界。u=0、相等绝对值、零符号、任意anchor和无穷时间尾均被覆盖。

对于u1=0、r=||u||∞，所得

\[
R(z)=(\sqrt{1+z^2/16}+z/4)^2(\sqrt{1+z^2/256}+z/16)^8,\quad z=r/t,
\]

是全空间密度的双侧比值。R递增，所以真实盒的h∞和band左端a给所有anchor-to-box线段统一常数R(h∞/a)。未改变body16未单位化方向、alpha1、Fourier(2π)^-5、eps=1/512、q=5/4、原anchor或盒。

085的面内二阶常数25/24及Taylor系数25/48同样独立审阅；真实一次项先平均得到D_μ，不误设偏心anchor的μ=0。消费形式为

  C=min(RU,1+||μ||∞/a+(25/48)M2 RU/a²)。

092未保存逐band旧RU，实际仅使用新RU，再与既有aggregate range/mean/self证书取交，完全合法；无需重新求旧exp/asinh，也无需假设192-bit新RU在数值上严格小于未保存的旧RU。完整G短时误差E按076范数经Minkowski合成，没有重复加064tube。

## 3. 89项catalog的不同路径

catalog_registration在计算前冻结来源和有限89项。全部56139节点实际只有3个h∞：2049/524288、2049/262144、2049/131072；与同一73个左端生成的所需集合恰等于catalog的89个x。与091旧六例{1/8,1,2,4,8,16}交集为空，未重跑旧案例。

094不导入生产ratio函数，采用不同缩放

\[
R(z)=\frac{(\sqrt{16+z^2}+z)^2(\sqrt{256+z^2}+z)^8}{16^9},
\]

用224-bit整数平方根和Fraction构成外区间。89项全部完整包含于生产192-bit区间；另准确核生产保存的根区间平方比较、正幂及最终向外量化。此为有限代数证书，适用全box/全时间来自解析定理，不是89点数值外推。

PID80674，0.286523秒；新H、Arb、log、exp、asinh均0。完整89项外端及源码SHA在catalog_results。

## 4. 全56139节点与全部新凸弦

whole_scope在完整092最终manifest冻结后登记，SHA `00358b59811799f208eb5457e65fdb52ecc9bbc83a56e28a676e37ac095fd26a`。代码SHA `00743a5f7b9056d4d498f99a52fde2749f180ddcec127349e153330dfb9eaffd`。固定14块（13×4096+2891）、8worker、无wallcap；每完整块持久checkpoint。

逐节点核原始raw共1,410,133,039 bytes及endpoint-only packets共465,506,822 bytes的SHA，并将point/endpoints/indices/weight/level准确比较。全部73原H端点和正band质量一致；最后H∞=0保留完整时间尾。这是输入复制和身份审计，不新计算H。逐节点从真实box端点和anchor重建h∞、||μ||∞；M2和全部旧self/range/mean逐字段、逐记录SHA绑定已接受088证书。

新range和mean用逐band整数/Fraction向外累加到2^-256，再确定真实准确总和在2^-160网格的floor/ceil bin。全部56139节点与生产的max/min交集量化值精确相等；没有一项模糊bin，没有使用预登记的准确Fraction总和回退。

对每个新li/ui，求li^5/ui^5的纯整数四次根并核四次幂方向，再证明递增凸函数z^(5/4)的真实mean凸弦。全部new_power_lower和new_chord_upper通过，L≤M≤U通过，selfL=max(oldL,weight·newlower)、selfU=min(oldU,weight·newchord)准确一致。旧self更紧时合法保留；没有把点值或mean本身当作q次幂积分。生产记录56138个selfU改善、56139个selfL改善，均由完整记录和交集验证支持。

## 5. 14个完整前缀与全域I

初末roots、全56139节点的全部不可变geometry/edges和55819 physical frontier逐项相同。独立准确核permutation weights、320父节点孩子权重和=16倍父权重、physical四维体积=1。闭覆盖数学沿用已锁087/082同一几何资格，未增加或删除空间点。

用递归DFS重建全部14前缀，区别于生产level逆序迭代。每次更新后明确检查所有未处理节点仍付旧088 selfL/U；有效L=max(selfL,children和)、U=min(selfU,children和)，UH逐前缀不增，所有receipt准确总数一致，末self/eff逐节点等于冻结final_state。数值bestcover与physical frontier不混加。

准确完整合成检查为

  B>E，(B−E)^5≥UH^4，I_upper^4≥B^5，cost_upper^4≥(5B/g)^5。

所有14前缀都通过，最终简cap6010267/10^11也满足cap^4≥B^5，且严格小于70301737/10^12。准确32门5^5B^5<32^4g^5为false。最终截时下证LH≈2.606975979384590×10^-5只属于∫H^(5/4)；没有把上包表达式的某个下端误称成实际I/C5下界。

全量实际PID83145、session12383、exit0，含加载29.3783835秒；全部14块通过。0新H/stencil/jet/锥，0Arb/log/exp/asinh，0新moment基线。逐块实际耗时和恢复指纹汇总在run_receipt。此次数学运行无失败或恢复；第一次尚未登记whole_scope的准备过程因手抄087哈希漏字被预检拒绝，未产生数学输出，按原manifest纠正后才登记启动。

## 6. 交付与复核边界

本目录仅保留本次报告、解析审阅、两份有限协议/登记、两条标准库Python核验源码、89项结果、全量结果和压缩汇总运行收据。14个已完成checkpoint及catalog末checkpoint已并入run_receipt后清理；必要来源使用准确路径/SHA，不复制约1.88GB共享端点。

verify_catalog.py及verify_whole.py是完整可读复现算法；已完成输出存在时拒绝重复成功计算，符合不重跑旧基线要求。全量复核须有whole_scope绑定的同SHA原raw、packet、initial/final state和来源目录；本证书不声称仅靠精简目录即可脱离这些输入重建H。若外部缓存缺失，应恢复同SHA输入，不能静默换anchor或重新产生不同端点。

只有本目录新增，无MATLAB、提交、暂存或推送，旧067/072/074/087及生产/他人文件未改。主保留一般L1传递和最终正式常数责任。本审计已闭合所登记完整092范围，不继续扩张实验。
