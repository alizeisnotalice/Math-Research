**100：混合 H₁/₆₄ 与 H₁/₁₂₈ 的闭覆盖独立审计（2026-10-26）**

全部有限门已通过。保留100804个已接纳H₁/₆₄叶，新194个H₁/₁₂₈叶的970个数值分量经不同整数部分分式路径实际核验；194个新Delta、1024曲率罚项及3104个实际闭角点全部严格通过。最终100998个通过叶、48638个内部split、52360个root形成精确闭覆盖，queue、inflight、fault及新H128障碍均为0。在已锁075/083/099的同模型理论下，这些证书给出完整G在全cube上的下界 (G\ge21/200000)。本目录将完整有限证据交主代理终审；没有更新I、Cn或C5，也没有将固定epsilon的H当作齐次函数。

模型仍是五维alpha1，16个未单位化body方向，共同时间、权重1、Lebesgue测度及Fourier逆因子(2π)^−5。不同epsilon只改变正时间积分下端；每个闭叶分别从 (H_{\epsilon_i}\ge g) 推得同一个 (G\ge g)，不混合梯度、不对分片epsilon求导。新目录audit_mixed100/独占；只使用本机Python，未调用MATLAB、生产kernel或作者的box/certificate/partition函数，未提交或修改生产。旧090、081及既有保护锁共133件保持原字节。

冻结098初态为 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098/initial_mixed_state.json`，SHA `19294a9bd493524fdc23ec2abe2bb3e70577f04b762b08d4fc57e78af5377480`。准备manifest SHA `fe03bfa691b6947a24f13a2f94885962332a77edfce087e2f9e01bc49a3ace15`。新batch0001 manifest SHA `0333c44e9d8bc29f2f1c2ac16ab38bbc935af8ae738f334a767812aed4a61928`，point_manifest SHA `293622876e07bb998227bad5b1a2f14af56a9727d17e3242023a43ca77d7d06c`。完整final_state.json SHA `73bc873f23f7b68f6f70132e590d28ebbe4cd619012ac9c81980f9fa500314af`。所有raw只在收到这些冻结SHA后读取；没有读取活future。

源门逐对象比对100804个原H64叶及原certificate，没有把旧叶重新贴成H128。全部184796个H64记录精确等于原180040、089正complete4752及4个历史H64障碍记录之并。typed cache由明确的epsilon外层字典分隔，初始H128桶为空。190个旧queue和4个counter inflight恰迁成194个新pending，原H64障碍的球与历史不等式完整保留；它们仅成为新epsilon下的待证叶，没有被当作通过或丢弃。

全部194个迁移mapper、闭盒边界、execution task、generic candidate和280个投影守卫独立重建；四个坐标相同的H64观测仍是H128 cache miss。物理盒使用root索引与每四位交织的path位直接构造dyadic整数边界，未调用作者映射。fingerprint098独立再生一致。所有根的叶前沿转成公共dyadic整数区间，排序后精确铺满单位区间，并验证prefix-free、每根Kraft=1、全部活split恰为严格祖先以及每次闭二分。该检查依赖真实路径与区间，未以叶数等式替代覆盖证明。

source_gate.py第一版在源/缓存/mapper/树已检查后，对非命中任务中显式old_H64_interval=null字段检查过严而停止，尚无H/log计算或数值成功前缀。保留其代码和登记，source_gate_v2.py明确核对null并完成有限源门。最终PID89827、10.100845083秒、exit0；source_gate_result.json SHA `7cebdfeab699d05e9f3cfb5ea09149bba4fb429cffc54a07c5846cff7374b13e`。源门与后续完整状态门均逐项验证原记录保留，没有重算旧H64值。

新值的独立路线保持已接纳081 arith.py只读，SHA `c0efda8a678b45fc608c8344dd6384e7d5a625b734fa54e5d1861c7fe8aa510e`。其原始epsilon为1/64，但输入改为2x，再分别将H区间乘16、四个真实物理梯度区间乘32。准确因子来自半群密度恒等式

\[
p_t(2x)=2^{-5}p_{t/2}(x).
\]

正时间积分作u=t/2替换，dt=2du，故

\[
H_{1/64}(2x)=\frac1{16}H_{1/128}(x),\qquad
\partial_{x_j}H_{1/128}(x)=32\,\partial_{y_j}H_{1/64}(y)\big|_{y=2x},\quad j=2,\ldots,5.
\]

这里同时改变输入与截时。arith.py允许首坐标2，给出的尾导数是实际物理导数；梯度因子32包含链式法则的额外2。scale_gate.py实际用Fraction核对194个输入的公共分母、全部投影和时间参数：若D、E分别是x与2x的公共分母，两者整数投影及时间参数同乘r=2E/D，且144E/64=r·144D/128。因子16、32与新罚项1/(16ε²)=1024均精确核对；零H/gradient/log调用。PID89747、0.066420583秒、exit0，结果SHA `db2904086846d8eeeec20c5bfa7747d4b1371030fe585c5673e5f2f3886e0363`。

numeric128.py在每个新输入2x上直接累加真实017 fan的3504个锥。独立部分分式算术只用Python整数与Fraction，log采用768位向外整数atanh级数、范围缩减及显式余项，π采用Machin公式及显式余项，保存512位向外dyadic端点；不导入DD、flint、Arb或mpmath。H与四个梯度的全部970个最终缩放区间严格内含对应生产球。017 fan及075理论是共享的已锁输入；本轮数值实现与生产EPS128 DD/Arb路径不同，081整数源码是已接纳的共享独审依赖，不声称重新独立验证它。

计算前登记194个冻结raw及代码SHA；12个独审worker先只派首点pilot，严格内含后才派其余固定至多16点块。coordinator逐点fsync、落盘后ACK；块与全局成功前缀持久保存，不重复成功点。所有输入首坐标2，旧面点首坐标1；本轮没有重算旧180040点或090的4756个成功H64面点。数值master PID90879、workers90889–90900，6.744841875秒，全部exit0并正常join。numeric index SHA `1d76ae172a8ce1193ea38f50d4f04d1bf9e3e1c220e63711bd52411b2b77059d`。

194个独立H128下端中的最小值位于原最强H64障碍坐标

\[
(1,2095105/2097152,16760833/16777216,134086657/134217728,2097153/1073741824).
\]

该点的独立向外区间满足 (10^6H_{1/128}\in[118.714525719529592399,118.714525719529592400])。旧H64在同点低于105µ的事实仍真实保留；两种epsilon对应不同积分尾，二者没有矛盾。这一中心值本身不替代闭盒角点门。

closed_gate.py逐记录检查raw→point manifest→task→returned→cache→certificate关联、epsilon128、kernel SHA、256位完整3504锥记录、物理导数次序、PID与calls。每个新证书重新用独立768位整数log下界检查 (0<\Delta\le\log(H_L/g))，用四个端点的交叉乘法验证同epsilon的gradient/H包络，在独立物理盒的16个闭角点用整数公共分母精确检查

\[
\Delta+\sum_{j=1}^4\min(\gamma_{j,L}v_j,\gamma_{j,U}v_j)
-1024\sum_{j=1}^4v_j^2>0.
\]

3104个角点全严格正，最小有理余量为

\[
\frac{52396691606863602612757118519868770672424416708790097061833191864106851163589}
{463168356949264781694283940034751631413079938662562256157830336031652518559744}
\in[0.113126665111544,0.113126665111545].
\]

证书机制、leaf_epsilon、target、罚项1024、point key与最坏角点均逐项核对。全局凹下界在闭盒上的角点门来自075/099，覆盖非generic边界、零投影、重复投影及盒交界；只对锚点要求generic，没有删去这些边界。

583个真实journal事件的canonical hash链独立核对并精确重建全部overlay：194dispatch、194received、194merge、1finish。每次dispatch恰为已登记194集合的当前叶，所有record均是新H128，零缓存命中；没有新增split。完整final_state的H64桶每个键/记录、100804个旧passed记录、所有split/root及保留元数据与初态逐项相同，H128桶恰等于194个新记录，passed尾部恰等于194个新证书。新增字段及fingerprint098_final独立再生一致，未只比较状态标签。

最终typed cache为H64的184796条加H128的194条，共184990条。最终100998个闭叶全部有严格证书，全部52360根精确铺满且prefix-free，活split48638个恰为严格祖先；queue/inflight/resume_jobs均空，fault与H128 obstruction均空。生产master90096、workers90116–90127和全部独审数值进程均OS查无。closed_gate PID91717、9.417805833秒、exit0，结果SHA `27ff2f48f9bd6be29ffe4d31fd01f457852b93986954a2237cac64fbf82caea1`。

全局结论的消费顺序是：每一闭叶以自身正epsilon给 (G\ge H_{\epsilon_i}>g)；root集合覆盖closed ordered face；同模型的tail置换和全signed-coordinate/cube-face对称性给全cube边界上的G下界。对非零cube内点，令r=||x||∞≤1，z=x/r，由完整G的−4齐次性得 (G(x)=r^{-4}G(z)\ge g)。原点G(0)=∞。这使用083/099中对extended正积分和全部闭边界的资格，不对H或混合H作同epsilon齐次延伸。

075、083、099与同alpha1的029原件SHA在qualification_result.json中逐一锁定。029的一般L1输入、同一个接收集及全部正共同cube尺度资格保留；它消费完整G的cube下界和完整I预算。本轮只完成g105输入，没有重算或联合I，故没有新的Cn/C5数值主张。qualification_result.json SHA `5b7d617b27bc3a737f7f85e01c22df5279f6b876d60960430d107c7254806b00`，主代理保留最终证明接纳责任。

复核先核manifest.json及其中外部冻结输入和133件旧锁，再读source_gate_result.json、scale_gate_result.json、numeric128_batch0001/index.json、closed_gate_result.json与qualification_result.json。所有数值块保存完整有理端点、输入/代码锁和实际worker回执；每一新证书保存精确margin及record/certificate hash。锁定后只读复核；需要重新运行时另开新审计目录登记，不能覆盖本轮成功前缀。
