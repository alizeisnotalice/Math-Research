# 174：固定 EPS64 完整 band 曲率的独立审计准备

结论：172 公式正确，可作为同一 Q64 的新增上包，与已资格 old155 q 及真实 H32 的 153 split 上包取 min。下述是解析审计与未绑定接口；没有读取真实 H 数值字段求 q，没有运行 band 乘法、角点或新 log/H/gradient。173 固定对象仍待作者冻结和主审，168/158 的旧基线资格不在本项预先认定。

## 解析链

固定 029/075/129 的五维 alpha1、16 条未单位化 body 方向、共同 Cauchy 时间、Lebesgue、逆 Fourier 因子 `(2π)^−5` 和物理 Euclidean 度量。正 Gaussian 混合给任意实 c

\[
\beta_c(t)=t^5p_t(c)=(4\pi)^{-5/2}E[(\det A)^{-1/2}e^{-c^TA^{-1}c/(4t^2)}]
\]

随 t 非减，并由处处有限的 `p_1(0)` 统一支配。对任意 band `0<a<b≤∞`，令 ν 为归一化的 `t^-5 dt`。β 非减、`f=t^-2` 非增，双重积分的协方差恒等式给 `Eν[βf]≤Eν[β]Eν[f]`。β 有界，f 在 t≥a 有界，ν 有限，所以 b=∞ 直接合法，或由有限尾极限得到同一结论。因 `p_t>0`，band 质量正有限。因此

\[
\int_a^b p_t(c)t^{-2}dt\le(H_a(c)-H_b(c))w(a,b),\qquad
w=\frac23\frac{a^{-6}-b^{-6}}{a^{-4}-b^{-4}}>0.
\]

完整尾 `w(a,∞)=2/(3a²)`，且 `H∞=0`。这是同一 anchor 的完整尾势差，不是对无限 Green 作相减。固定 ε>0 时 `0<Kε≤ε^-2 Hε<∞`，Hε 处处正、平滑。S1/S2、交叠、零/重复谱投影和边界均包含；这些集合只可能使完整零下端 Green 无限，不影响本证明。没有要求 Cauchy 绝对矩、对 β 求导或截掉无穷尾。

原 donor 的 73 个全时间键中，`t≥1/64` 的确有 49 个，含 t=1 及其后完整无穷尾。相邻 band 的端点交集是时间零测集，有限正和恰支付 `[1/64,∞)`。取可靠端点

\[
m_{j,U}=H_U(t_j)-H_L(t_{j+1}),\quad
m_{last,U}=H_U(1),\quad K_U=\sum_{49}m_{j,U}w_j,
\]

则 K≤K_U。合法区间自动给每个 mass upper>0；无需额外要求区间计算出的 mass lower 严格正。若 upper 非正则是源不一致，应停止。

同点原 g H64 与 donor H64 区间首先核准确非空交叠；可以在单点相交，不能静默 clip。不相交表明至少一个源或绑定有问题，应停。两者真实下端可取 `L64=max(gL64,donorL64)>0`，于是

\[
\boxed{Q64=K64/(16H64(c))\le K_U/(16L64).}
\]

正分子上包除以分母下端才保持上包。新符号诊断只验证新的方向残差

\[
\frac{K+d}{16L}-\frac K{16H}
=\frac{dH+K(H-L)}{16LH}\ge0
\]

（d≥0,K≥0,H≥L>0）。它不同于 153 带负号的尾质量比。153 分支应是

`A(7/10−8·donorL32/(15·original_gU64))`，`A=1/(16ε²)`。

这里必须是真实同点 H32 下端及**原 g H64 上端**，不能把 L64 填入该分母。两个分支分别上包同一个 Q64，才可与同 cutoff/point/metric 的 old155 q 取 min。最终只要求 `0<q_new≤q_old`，不沿用旧129的 q≥256、A/4 或 A/6 下限。有限区间误差可能使 band 分支不改善，理论正确不代表实际 strict gain；保留 min 仍安全。

129 的固定锚点 Jensen 资格独立保留：固定 t 的全空间 log-Hessian 下界经空间 Taylor 给

`p_t(c+u)≥p_t(c) exp(∇log p_t(c)·u−|u|²/(16t²))`。

按固定 anchor 概率 `p_t(c)dt/Hε(c)` 对有界指数变量作 Jensen，score 的 `|u|/t` 支配允许完整时间积分交换，得到

`log Hε(c+u)≥log Hε(c)+∇log Hε(c)·u−Qε(c)|u|²`，对所有实 u 成立。

这不是沿 c→c+u 路径冻结一个变化后验 Hessian，也没有把 t 依赖因子提出积分。故 anchor 在接收盒外、线段穿过奇异子空间均合法。只换 q，Delta 和全部四个物理 log-gradient 分量仍来自原 H64 full jet；donor H32/H64 仅供应质量，不能替换任何梯度或 log 分母。

以 directed gamma 的线性 min 加 `−q_new|u|²` 得到的安全函数凹；一个 anchor 的全部 16 个实际闭角点通过即可支付整盒。不得临时为不同角点换赢家。对旧已测的完全相同 source core/receiver，只有 q 严格下降才可重测；新 source core 仍须按完整 typed source、target、Delta、gamma、q 与实际盒去重，不仅按 point/canonical id。

## 两条准备中的权路径

`arithmetic.py` 仅定义函数，没有对真实或人工时间点执行。第一路使用差幂 w。第二路令 `z=1/t=B+du`，`B=1/b,d=1/a−B`，直接正系数积分

\[
P_n=\sum_{k=0}^n {n\choose k}B^{n-k}d^k/(k+1),\qquad w=P_5/P_3.
\]

这里 `t^-5dt` 对应 `z³dz`，`t^-7dt` 对应 `z⁵dz`；公因子 d 消去。该路径没有使用相同差幂重排。B=0 的完整尾也直接合法。未来每个授权 band 必须两路准确 Fraction 相等；目前没有实际求一个新权或 q，也未重放 153 已完成的全尾/两段权重符号基线。

## 已核静态 schema 与尚未绑定项

只读代码实际确认：151 minimal jets 的 `values` 有 5 个区间，`values[0]` 为 H，随后四个是物理 H 导数；coeff 用 `jet_id` 和 `raw_record_sha256` 绑定原记录。151 的旧 numeric checker 说明 Delta/log-gradient 来源；本项只读源码，没有调用它或读取 jets 数值。

092 compact donor 是 `{key,point,endpoints,original_raw_sha256}` 投影，`endpoints[j]={t,H:{lower,upper},...}`；095/103 为保留完整 evidence 的同名端点。173 作者报告拟用 411 个092、5个095、39个103，共455源；该计数在本项仍标为**作者准备元数据，尚未独立读取455记录绑定**。092 应连接其 nodes 的 original raw 路径/SHA，并验证原始 endpoints 完全一致；095/103 应核 producer kernel/registration/fingerprint、bodytable、point 和同模型完整73time。不能因旧 raw 缺 model 字段就补造一个原字段，只能引用已冻结 typed 外层来源。

`binding_template.json` 明确 173 manifest/schema 尚为 UNBOUND，实际 case 注册数为0。未来上限为 455 q、22295 正 band products、至多35872 source/receiver pair各16角点；该模板不是实际运行授权。原168/158 baseline 真值/去重消费门仍须独审/主接纳，不从待审标签推断。

167 transport 只证明旧球载荷身份与无损搬运，不是新的 H 数值证书；455 donor 的 H 真值仍共享原已认证数值来源。同 cutoff 交叠也只是单位/来源的一致性诊断，不能证明球含真值。

一般 L1 与全尺度输入链没有改变：有限闭叶完整覆盖 F 后才可得相应 H/G floor；固定 cutoff H 不作齐次扩展，只通过 G 的−4齐次性和对称性扩到 cube。029 的有限正源、同一接收集、一般 L1、全部正共同半径接口仍适用。局部新 q、若干通过 pair 或旧 target115 pending 不构成完整 floor/I/C5。本项未消费全域证书或合成新常数。

唯一新数学工具运行是预登记后的纯符号分母残差：PID **8177**，exit0，0.014008458005264401秒，identity difference=0；真实 H 字段/q/bandproducts/corner 为0。源 SHA、代码 SHA、时间键 metadata 和结果见 registration/manifest；不改167/160/147冻件，不用MATLAB/Astra/下界，不提交或推送。主保留最终资格及证明责任。
