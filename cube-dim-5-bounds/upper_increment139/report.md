# upper_increment139：直接夹逼正一维短时增量

**仅准备，尚未运行实际 J。** 固定唯一实际输入：ε=1/512、134 已合格 Cbar 上端、132 已合格 M portable cap、014 完整 H512 积分上端、完整 g=21/200000，与135实际输入相同。本轮只改进对133正一维积分的求值，不重做135/136标量、M优化、旧toy、范数或Gamma。

主已在第015轮完整接纳 C5 上界约31.93801968477627（准确上界以 `round015/primary_upper_acceptance.json` 为准，主提交6eb0625）。该结果是本轮已知参考，不在此重算。139尚无新的 I 或 Cn 结论。

## 1. 唯一积分与共享资格

沿用133的完整面弱次主序、非减凸函数比较与正积分变换。给
\[
D=\frac{Cbar}{2\epsilon^2},\quad W(z)=D+(M-D)z^4,\quad h(x)=x^{1/4},
\]
所需辅助增量是
\[
\boxed{J_{\rm aux}=\frac{5\epsilon^2}{2\pi}\int_0^1z^2h(W(z))dz.}
\tag{139.1}
\]
前因子是5mD=5ε²/(2π)，不是135归一为概率后的5ε²/(6π)。两者分别配合权 z²dz 与概率3z²dz，不能混用。

共享的原证明保证：真实Rε的凸增量不超过 J_aux，且完整
\[
I\le I_{H,U}+J_{{\rm aux},U}.
\]
Cbar要求整个F与0<t≤ε的短时峰值，M要求同ε、整个F的Hsup；IH要求同ε的完整面H^(5/4)积分。133的集合质量/stoploss/非减凸/MCT论证覆盖S1/S2的canonical无限代表，139不删空间点或奇异邻域，也不改变029的一般L1和全部正共同半径接口。

139的夹逼区间只包住固定保守参数 Cbar、M 的 **辅助预算表达式 J_aux**；其下端不是实际R增量、实际J或完整I的下界。输入π区间只用076保存端点，不重算π或Gamma。

## 2. 每面板的正权弦与Jensen夹逼

固定N=128，在面板[a,b]=[i/N,(i+1)/N]置
\[
w_i=\int_a^bz^2dz=\frac{b^3-a^3}{3}>0,
\quad
\mu_i=\frac{\int_a^bz^6dz}{w_i}
=\frac{3(b^7-a^7)}{7(b^3-a^3)}.
\]
令
\[
\theta_i=\frac{\mu_i-a^4}{b^4-a^4}\in[0,1],\quad
A_i=W(a),\quad B_i=W(b),\quad
\overline W_i=D+(M-D)\mu_i.
\]
准确有理恒等式给
\(\overline W_i=(1-\theta_i)A_i+\theta_iB_i\)。

h为非负凹函数。面板上
\(\lambda(z)=(z^4-a^4)/(b^4-a^4)\in[0,1]\)，且
\(W(z)=(1-\lambda(z))A_i+\lambda(z)B_i\)。由弦下界并积分，
\[
L_i:=w_i[(1-\theta_i)h(A_i)+\theta_i h(B_i)]
\le\int_a^bz^2h(W(z))dz.
\]
再对条件概率 z²dz/w_i 用凹Jensen，
\[
\int_a^bz^2h(W(z))dz
\le U_i:=w_i h(\overline W_i).
\tag{139.2}
\]
全部系数、参数与权为有理数；只有四次根需要有向包围。M<D时W随z递减，但凹函数的弦和Jensen仍适用，因此没有单调方向的隐藏假设。

将端点根取可靠下端、均值根取可靠上端，因全部权非负，得到向外的离散和
\[
L_N\le\int_0^1z^2h(W(z))dz\le U_N.
\]
给π∈[πL,πU]，完整辅助预算夹逼必须是
\[
\boxed{J_L=\frac{5\epsilon^2}{2\pi_U}L_N,
\qquad J_U=\frac{5\epsilon^2}{2\pi_L}U_N.}
\tag{139.3}
\]
两端不得使用同一个π端点或反向配对。

## 3. 与135损失的关系及严格误差方向

这是有限面板的一维积分认证，不是重新求135整段Jensen。数学上，Σw_i=1/3、Σw_iμ_i=1/7，故对离散概率3w_i再次用凹性可证明其**准确**上和满足
\[
\sum_iw_i h(\overline W_i)
\le\frac13h((4D+3M)/7).
\]
因此精确面板上包不会大于135整段Jensen；细分的方向也正确。外向根舍入后的实际改善必须由本次新结果与旧保存上端准确比较，不凭此符号关系预先宣布数值改进。

有限误差无需估计某个未知真实核。若u0=min(D,M)>0，h二阶导绝对值≤3/(16u0^(7/4))。每面板Jensen与弦的差可由插值误差界控制为
\[
0\le U_i-L_i\le\frac{w_i}{8}\frac3{16u0^{7/4}}|B_i-A_i|^2.
\]
又|b⁴−a⁴|≤4/N、Σw_i=1/3，所以根为准确值时
\[
0\le\sum_i(U_i-L_i)\le
\frac{(M-D)^2}{8u0^{7/4}N^2}.
\]
192-bit dyadic根区间每端误差至多δ=2^−192，两侧离散和另至多增加2δ/3的宽度；π输入宽度按(139.3)直接向外传播。实际证书保存完整准确区间与其宽度，而不是只显示一个数值。

本轮只登记N=128一次，不自动重复N=1/135或倍增面板。若该固定批上端未改善已接纳I，就报告该结果并停止；不擅自追加新精度、面板或H求值。

## 4. 唯一输入的准备与字段锁

`prepare_inputs.py` 只做字节/JSON字段绑定，不做任何J、D、根或面板的物理数值计算。实际PID **57264**、exit0。原135输入SHA是
`c6eb027efee798f140f6cb7362721b246de257a722b58cf65ba57ea1a6ad426c`。

逐项复制并核对：

- Cbar：134 `condition_results.json:Cbar_upper_256dyadic`；SHA `1cf1aa32b44eb580226692bfe37ab40d2ef96deeda4edc155f16a2f0a48cb024`。
- M：132 `portable_components.json:objects.analytic_stationary.M_upper_portable_cap`；SHA `2ec64713d80f2ecc791a7207a706d035b4c6723c24868cc44383c4d91ff9855b`。只用原准确有理portablecap，不从display构数，不重优化。
- IH：014主接纳 `complete_H512_integral_upper`；SHA `14374097c6d164fcbaca56a4b6d003069e92e5e16ffa990a2ebb6cb470d2c219`。
- π：076 `derived_constants.pi_outward_384_dyadic` 的两个原端点；SHA `caeda104af7973ccb153daf44cfd3711cac332ae262f5960415263a32f102ec3`。
- g与旧I：主015完整接纳记录；旧135保存cap只用于新上端比较，旧135结果SHA `8b0623db708b493cf7d196e5f8dc0f88724d401154a7c888da95987568502686`。不重新运行旧根或标量。

新固定 `inputs.json` SHA：
`fe6dac94e6de5c7fd19d1a1da2b32279945b54b6d9aae428f1c8bd42796315b7`。

完整来源锁与资格元数据已嵌入输入，source_binding只记录复制和来源一致性。根和面板实际运行数均为0。

## 5. 生产与可携checker的有限接口

待主审的 `run.py` 固定一个case、128个面板、192-bit根与257个根原子（129个面板端点、128个条件均值）。单进程、无pool、无wall alarm；实际根函数只用两次整数isqrt，并以准确四次幂检查向外方向。

每个完成根原子与每个面板后都原子写入外置state，保存来源/代码/输入/计划身份与成功前缀。故障保留prefix与独立fault记录并退出，不自动重试。只有显式`--resume`且所有锁、输入副本与未完成state一致时才接续；完成根只校验已保存幂包围，不再次调用isqrt。完成状态拒绝续算；不自动加密或创建新实际case。

最终packet包含输入和计划的**逐字节原副本**、registration、257根区间、128面板证书与结果。`check.py` 不导入producer，不使用isqrt，只用stdlib Fraction核对：

1. 固定预期输入SHA，而不是仅信任包内相互引用的SHA；固定producer SHA，校验plan/code/checker/input/registration/certificate/result绑定。
2. 每个根原子的准确argument、非负dyadic端点与其四次幂包围。
3. 以整数指标差重新推导w_i、μ_i、θ_i，核每个正权面板下上端与准确和。
4. πU/πL前因子、辅助J区间、完整I上端、新成本四次幂与对旧保存上端的比较。

因此不能把packet换成另一组自洽的C/M/IH数据而通过这个固定case checker。它可在没有原H/R生产代码的环境中运行；外部来源的数学资格仍是主、原132/134/029与独立审计的共享依赖，不由有限算术checker重证。

checker验证的是一个新积分证书，不重跑135/136的求值。producer始终标记实际J/R/I下界主张=false、global_Cn_claim=false；任何新的正式界仍须独立checker、主数学接纳和正确归档。

## 6. 仅运行的新符号核验

本机已有Python/SymPy实际运行 `symbolic_panels.py`，PID **57265**、exit0、约 **0.201秒**。只核六个新恒等式：面板质量、条件四阶矩、仿射均值、积分弦权、总质量望远镜与加权均值望远镜。没有给符号赋实际C、M、ε或N=128，没有计算根或实际J，也没有重跑旧toy/范数/Gamma。

新实际J case、物理根原子、物理面板运行数均为 **0**。生产与可携checker只做语法检查，未执行。主完整审阅具体冻结代码与输入后，才可授权这一唯一有限新积分批。

只拥有并新增 `upper_increment139/`；冻结128/133/135/execution135保持不变。低资源单进程，不使用共享MATLAB，不争抢主15核H/PF生产，不启动下界或Astra，不提交、不推送。
