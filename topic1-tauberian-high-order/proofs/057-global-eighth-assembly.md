# R057 候选证明：同一新八层的精确八阶结算与全部重数传递

2026-11-05。**已供前提下的完整总装候选，供非作者严格复核；未发布定理、未改正文或仓库。** 只写本独占 scratch。R053、R055 以及二维 H56 已正式供货，本稿新增的二维补集重放和全输入总装仍待独立审核；不凭局部价格或标量核的通过直接宣布结论。

本稿真正新增的连接是：同一 N8 的正常模型与完整补集只减一次 J；实际删除质量单独保留 sharp 导数；源缺量和重复量在八层中正抵消；最后直接移实际 M，才计算全输入分母。全部多线性重数的精确八阶传递另给五阶参数截断，不能继续使用旧 quartic 截断的同阶误差。

## 0. 冻结依赖与供货状态

固定 n≥2，记
\[
 a=v^n=2^nt,\quad v=2t^{1/n},\quad A_v=(1+v)^n,\quad
 k=n^2/27,\quad B_n=\binom n2,\quad K_n=B_nk,
\]
\[
 c_{8,n}=-13/18-n(n-2)/3+400/243,\qquad\Gamma=(10/9)k.
\]

| 依赖 | 本稿所用接口 | 当前状态 |
|---|---|---|
| R048 同族前缀 | 同 N8 的第七前缀定量粗价、第六前缀领先上包及原八层公共账本；各层用本层实际对象重放 | 数学模块源 `f19499528e7fe5d322cf2a2d0333d7743c2fd6ffa3ac059eaa5798adc34b98cd`；完整正式审 `333e7ffcbfe200b0e8adbd0bb013663c91c6071c88f5ecad41facd4659d2f1a9`，已实核并读同族前缀范围 |
| R049 | 竞争输入的同 N8 稳定性、固定正常集合、raw、端移、面量 O(v^8S)、D_N=O(v^2S) | 已审供货，源 `5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4`；正式审 `e8018aac5929e2ac83ae926edafdefa98ccf580cb7f4f3e9ed03cf2cda1a128c` |
| R050 | n≥3 的同 N8 全部高支撑 gross=o(v^8S)，不再消费第二份 J；n=2 此交集为空 | 源 `22f848879b3c9ba839edb5d307a205ae3245e5ba3a89a99f7ec2ce2974010748`；正式审 `68e9532a8b80926806b0ed11219c7630a46fe7aa4a9d72c0a5e94a73b520377e`，已实核并读固定竞争族范围 |
| R051 | 整个有效源 chart、固定实际原目标帽与 root、同原见证高度必要性及真实侵蚀层 | 源 `37bbb0ec63fac1dce33fb9883e14499488baef490341df32a6fdfb0da427c979`；正式审 `535a2ccdf58676dac81bd8d8cd5862f0ab9cdbc5bde1e14beeedace0fd1382ba`；已审供货，n≥2 同层固定有界 packet |
| R052 | 同一实际帽算子的反单调第二反馈与一次强 CDF；packet 可含固定有界非正常目标 | 源 `d655c1e160c88914774248c81688c6bdde3f27614387c4eb7f0b9c1a0a0a5296`；正式审 `2ec88399d618522fc0a57332b91597121781e770d47359d1f54441f766df6985`；已审报核合 R051，历史 H51/H52 已解除 |
| R053 | **只用 (13M)**：未减任何实际 J 的正常模型 P_N≤q_N W_N+o(v^8S) | 源 `d83dc772fc2c1ca5698a39eb2bd15702395d5e8238b11209fd2b6e7765b1e642`；完整审 `a6ebc52f1274e1d5cb00eadc7aad1ce8dd791d8713006a880f878e2ebee8cbf0` 与第二审 `7ab3b0eec28e61df3f5324bbef8c0425277343d80379c59da13e84d0f846a0f0` 均已实际读取相应接口并核 SHA |
| R055 | 同 P_N 的整源共同高度，完整 J 一次，补集任意固定 δ v^8S+o_δ(v^8S)，保留槽欠量 | 源 `0684d01f3665a5a064b4786b6eb1d78ea5780ccdae64f96df3ad01647224a1a1`；完整审 `aede0585bc545d564484cc2fca82d8601c254794e4c3748ae73b744fe5dcbfae` 与第二审 `f46f1e0e5fef543c20eec7748cebdb3e11100aca147fad9ed66873d5ccaf9a0d` 均已读相应接口并核 SHA；正式范围 n≥3 |
| 二维 local 重放 | 本稿 §2 给出 R055 §§1–9 在 n=2 的原位证明，使用已审实际四角/四槽桥 | 新总装的一部分，**待本稿独审**；不是声称旧 R055 已审范围自动扩大 |
| R054 内部正预留 | 同 N8 的 P_I+P_∂=O(v^8S)，二维 boundary/env 晚于八阶 | 已实际读源 §2 与完整审对应段；源 `bc74ac7f20b63be1387bf1f465316fb61e00513b6567002246ac32095b4f00d1`，完整审 `fc43cc5d4b3361d60daaf080e728f36ad7843a4bae866b7e6e8b5b1855ae35c2`；本稿显式新增调用，不谎称 R055 曾使用它 |
| H56（已供应） | 二维 inside sharp 常数 2K_2，且实际 fresh 基线的删除系数为 1+O(v)，不破坏原 D/f 账本 | 新源 `18fa60b96a6b2a2be7092d06dc59b4741a049e4214941323a33f846dad12ffc0` 已实际全文读取；完整终审 `a87f991f05c0e79ac661cf89fd388be0f3ffaf51a9c98d738b92ae90e0cc1ab4` 全文读取，其 §8 明确供本总装；第二审 `2e5f13ea480534ea1d3ef440fc8fa18f267059de0da455a858b3f2e20e8785bd` 的 §8 已读并核新源字节绑定 |
| R041 下界 | n≥3 所用 J_n 的实际有限输入下界，对每个充分小实参数；不是上系数或安全模型 | 源 `bba23ac1dba9cfe6c091426595fa162bebe3c585d3afce6d9d2968f497d11845`；几何正式审 `7d0e7808f21d15ee35b6465631301940913299e0f8e314307359ee5a68d0832f`；代数正式审 `9c3c97f20bd94c65206971148f4b20c571799ab30797715fc40e7d23ebee5130`，均实核其实际有限族范围 |
| R054 下界 | 二维新两尺度 lower J_2+8/27，同样对每个充分小实参数 | 源 `bc74ac7f20b63be1387bf1f465316fb61e00513b6567002246ac32095b4f00d1`；完整审 `fc43cc5d4b3361d60daaf080e728f36ad7843a4bae866b7e6e8b5b1855ae35c2`；lower 第二审 `625f1a0301e39f34694d5a979c97549faaa79a4a03e20c16c12b4a7d88e3e49b`，已读范围与严格实际见证段 |
| 既有精确多线性约化 | 共享同一半径、公共交集正有限测度归一化；σ_m≤log2<1 的全 m 共同有限分支 | 既有 main 源 `b976cf01d89b73ff2ec540a19aebc5efbc44848cd2d2fc6bf271667eb4eb63f5`；正式 main 审 `7f42e40b738bcce1f389c8f730ed0566069a7ffaefa4156d919f5ef2b2929ec9`；本稿仅接用该既有约化 |
| 新代数核 | 实际分母、ordinary b8、八层中心化、五阶 σ、fresh 归一化 | Terra 已实跑并返回 WL 与独立 Fraction；本稿已实际读取代码/结果，详 §10；仅作标量证据 |

表列新增源与正式审报已实际读取相应接口并逐文件核对完整 SHA；R048 同族前缀、R050 高支撑、R051/R052 必要性与迭代、R041 实际下界以及既有精确约化均显式列出，不以 R053/R055 的引用代替供货版本。H56 正式供应后，二维 inside 不再是未完成前提；其常数来自已审实际引理，绝不是把旧 dyadic 的任意常数 C 改写成 1。**本稿的新二维补集重放和全局总装仍须独立审查，尚未发布精确八阶定理。**

## 1. 实际竞争序列和不可改动的同族结构

仅需考察任意固定 B<∞ 下的实际输入序列
\[
 v\downarrow0,\quad V=|E|>0,\quad
 K\subset\mathcal H_n(1-t;E)\setminus E\text{ 紧},\quad M=|K|,
\]
\[
 1+M/V\ge P_{7,n}(t)-Bav^8,
 \quad P_{7,n}=L_n+K_nav^6-(n+26/9)K_nav^7.
\tag{1}
\]
先按原 P/S 图固定 I、回收/边界价格和全部原见证；一次完成八层有限递归，统一 A_8 后取同一个强 CDF。各层源 T_(i,h)、F_(i,h)、目标 G_(i,h) 均来自本次族。写
\[
 \mathcal T=V+M,\quad S=\sum_i r_i^n,\quad
 R_{i,h}=(g_{i,h}/a)^{1/n},\quad
 m^0_{i,h}=|T_{i,h}|,\quad\mu_{i,h}=|F_{i,h}|,
\]
\[
 \delta_{i,h}=m^0_{i,h}-\mu_{i,h},\quad
 \mathcal D_h=\sum_i g_{i,h}-M',\quad
 \mathfrak d_h^0=M'-\sum_i m^0_{i,h}.
\]
外部源求和恒有
\[
 \sum_i(g_{i,h}-\mu_{i,h})
 =\mathcal D_h+\mathfrak d_h^0+\sum_i\delta_{i,h}.
\tag{2}
\]
原有限递归给
\[
 \mathcal D_h\ge f_h\ge0,\quad
 \mathfrak d_h^0\le f_{h-1}+aP_\partial,\quad
 f_0=A_0vM',\qquad1\le h\le8.
\tag{3}
\]
同层领先质量为
\[
 S/\mathcal T=1+O(v),\quad M/V=a[1+O(v)],\quad
 M/a\le S\le V/(1-a)\le\mathcal T/(1-a).
\tag{4}
\]
第八层正常集合、正常模型、旧正常实际层、槽分配和 R049 稳定常数先固定。后选补集 δ、有限 T、m、软权，只改变测试 packet 和消失误差常数；**不改变 I、源、目标、root、原见证、CDF、正常集合或正常模型**。每个辅助参数先固定再取 v→0，没有 T(v) 或 δ(v)。

## 2. 二维补集的原位重放（本稿新增，待总装独审）

R055 正式源仅 n≥3。本节在 n=2 重建其 local 接口，不能仅把维数改为 2。

在同一二维 N8 中，所有二轴边的轴对都是 {1,2}。global-min-face chart 的最浅轴必属于这对轴；不在 pair 的图完全为空。因此 R053/R052 的正常 B_N、Y_(2,N) 是同一实际二维图，三近面名义项与其退款核均为空。这不会删掉其他实际源图。

给固定 δ>0，沿 R055 次序先选 T>T_N、再选 0<m<c_N/2。正常与有界非正常共同 packet 的实际高度为 h_A=h_N+h_B≥0，均是同一见证中实际目标帽质量之和；帽的 transverse 乘积取 1，真实原始目标锚和有效源锚照旧区分。

本节的量从二维实际图原位定义，而不是沿用一个高维投影：λ_v=(1+v)v/a，
\[
 \mathsf M=\lambda_v\sum_{\text{二维 charts}}\int h_B^0(y)\,dy,
 \quad\mathsf P_N=\sum_NI_e-\lambda_v\sum_i|Y_{2,N,i}|,
 \quad W_N=\sum_NR_e^2,\quad W_B=\sum_{\mathcal B_{T,m}}R_e^2.
\]
所有定义保留实际完整横宽（二维为空积 1）、实际 R_i 分母和两个有向端，不把 M 换成带小权的角量。

取全 packet 的同一非负侵蚀 ε_i(y)，令 B_N^-、B_A^- 分别为 min(H_i,(h_N−ε_i)_+) 与 min(H_i,(h_N+h_B−ε_i)_+) 的子图。对全部 H,u,w,ε≥0，
\[
 0\le\min(H,(u+w-\epsilon)_+)-\min(H,(u-\epsilon)_+)\le w.
\]
故单一 raw 转移常数仍精确为一：
\[
 \mathcal R_N^-\le\mathcal R_A^-+(1+v)^{-1}\mathsf M_{T,m}
 \le\mathcal R_A^-+\mathsf M_{T,m}.
\tag{5}
\]
不要求正常/非正常层物理分离，也不允许在正常价格中先减一次全 J。

在非零二维目标项的支撑内，唯一底变量距有效近面为 O_T(v^3R_i)，纤维总高度为 O_T(v^4R_i)，而盒半宽是 cvR_i 的固定正倍。因此真实 fibre 顶端 H_i≥cvR_i 大于总高度。没有非临界自由底变量，R055 的 upper-clip 损失在该目标项上为零。M_3 亦为零；完整有界非正常名义成本只有两个实际排序的 M_2，仍保留两方向不等根和真实截幅。

旧正常实际层的锚定、R049 raw=O(v^8S)、目标核全族范数 C_Tv²/a 仍给
\[
 \Delta_A\le(1+v)\mathcal R_*+o(v^8S),\quad
 \mathfrak F_B\le C_{T,m}v\mathsf M_{T,m}+C_Tv\mathcal R_*.
\]
由增量恒等式及最终小 v，得到 M_(T,m)=O_(T,m)(v^8S)，且其归一化上极限≤固定 C_*，不随先选有限 T,m 增长。退款和 clip 剩余为 o_(T,m)(v^8S)。

因此实际增量和第二反馈满足
\[
 \lambda_v\sum_i|B_A(\mathscr D)_i|
 =\lambda_v\sum_i|Y_{2,N,i}|+\mathsf M+o_{T,m}(v^8S).
\]
对 B_A^- 用已供应的二维整源强 CDF **一次**，得
\[
 \sum_iJ_i\ge\lambda_v\sum_i|Y_{2,N,i}|+\mathsf M
       +c_R\mathcal R_A^-+c_E\mathfrak E-o_{T,m}(v^8S).
\tag{5J}
\]
这从实际层的高度可加性而来，不是 normal-only 扣全 J 后再加坏层成本。

全域非正常标量强制性在同一二维正部模型上保留 γ_ν 和 c_B。宽度误差依旧为 v^6Σξ_iR_i²=O(v^(17/2)S)=o(v^8S)，根比替换为 CvM=o(v^8S)。故 resolved 非正常净价有固定折扣，并保留一份完整 M，而非只有 b_eD 的加权成本。

具体使用 R045 已审全正部的强制式 G≥c_ν(b+D)，其中 G、D 取同一 x,y,ζ 归一化；二维 b=1。另对原始 p,q≥m 的完整两方向名义高度积分，正部凸性给固定正的 c_m v^6R_e² 下界。二者只在准确回代各自尺度后比较，不把 raw p,q 的高度成本与 G 中的归一化 D 混同。故可在正常窗、T,m 之前固定 γ_ν,c_B>0，使实际有界坏边满足
\[
 \sum_{\mathcal B_{T,m}}I_e-\mathsf M+c_B\mathsf M
 \le(1-\gamma_\nu)A_6(v)W_B+o_{T,m}(v^8S).
\tag{5B}
\]
它对极短边不使用，极短边由下面的 gross 软容量处理。所有后选窗误差只乘消失的 v 或全局宽度 L1，不留下 C_Tv^8S 的未知常数。

**二维容量来源单列。** 使用已审 R045/R048 的实际四角/四胞与质量根桥：每个正常方向槽唯一，实际 clean 占用给
\[
 W_N=\mathcal C-D_N+O(v^3S),\quad
 \mathcal C=\sum_iR_i²,
\]
全部软边仍满足 W_soft≤C+O_Csoft(v^3S)。这里不引用仅为 n≥3 证明的 transverse helper 或反用软容量。由此 near-zero 角 gross 的归一化上极限≤TmC_*。

长二轴条直接沿已审二维同见证轴帽证明：目标正常层在长帽内的比例≤C(v³/P+v²)，源长短条至少有一个轴的近面深度>K_Nv³R_i，因此避开二维正常二轴小管。条不必避开新整源大管；root、任意长短比、完整实际短边都保留。固定重数给 long gross≤(C_L/T)R_*+Cv(R_*+F+端移)，故归一化上极限≤C_LC_*/T。C_L 来自固定正常小管，不随新 T,m 改变。

边界相交端用实际缺量分解和有限度：Σ_Adj(∂I)ε_iR_i²=O(v^7S)，故其 gross=O(v^9S)。二维 q≥3 的高支撑交完全为空；不借 R050 的高维论证。零/一临界 pair 依原强面容量无正体积交。各分组覆盖全部外部 pair。

因此先 T 大、再 m 小，使 C_LC_*/T 与 TmC_* 之和小于预定 δ 的固定份额，最后取 v→0，同一个整源强 CDF 只减一次完整 J，给 R055 (6) 的二维版本：
\[
 H_8+c_Dv^6D_N+c_B\mathsf M+c_R\mathcal R_A^-+c_E\mathfrak E
 \le q_N\mathcal C+\mathsf E_N+\delta v^8S+o_\delta(v^8S).
\tag{6-2D}
\]
此式可直接由 (5J)、(5B) 和小角/长角 gross 相加取得：先有 H_8+正剩余≤P_N+(1−γ_ν)A_6W_B+δv^8S+o_δ8，再用 W_B≤D_N+O(v^3S)、W_N=C−D_N+O(v^3S)。因为 q_N−(1−γ_ν)A_6≥c_Dv^6，槽欠量的正价只保留一次。模型 excess E_N=P_N−q_NW_N 是同一个正常模型。
这只是 local 接口。它没有消除 inside 删除，也没有从 (6-2D) 推出 H56；二维 inside 恰是下一节保留的独立成本。

## 3. 同一个 P_N 接合，槽欠量不重复收取

取 R053 (13M) 的整体标量
\[
 q_N(v)=\frac{k}{(1+v)^2}v^6
       +\frac{\Gamma}{(1+v)^3}v^7+k c_{8,n}v^8.
\tag{7}
\]
它对小 v 为正且 q_N=A_6+O(v^7)。固定它后，R055 的模型 excess
\[
 \mathsf E_N=\mathsf P_N-q_NW_N\le o(v^8S).
\]
将 **这一条模型式** 代入 R055 (6)（n≥3）或 (6-2D)，只取得一份 c_Dv^6D_N：
\[
 H_8+c_Dv^6D_N+c_B\mathsf M+c_R\mathcal R_A^-+c_E\mathfrak E
 \le q_N\mathcal C+\delta v^8S+o_\delta(v^8S).
\tag{8}
\]
不再加 R053 (14) 的第二个 D_N，也不加 R053 (13) 中已减全 J 的 normal-only 式。左侧的 H_8 始终是全部 I_e 减完整 ΣJ_i 一次。

完整实际槽乘积给
\[
 \mathcal C=B_n(1+v)^{n-2}\sum_iR_i^n+O(v^{7/2}S),\quad
 \sum_iR_i^n=(M'+\mathcal D_8)/a
 \le\mathcal T/A_v+O(v^6S).
\]
误差乘 q_N 为 o(v^8S)，所以可取
\[
 Q_{8,N}(v)=K_n\left[
 \frac{v^6}{(1+v)^4}+\frac{10v^7}{9(1+v)^5}
 +\frac{c_{8,n}v^8}{(1+v)^2}\right],
\tag{9}
\]
得到 H_8≤Q_(8,N)T+δv^8S+o_δ(v^8S)。最后一项也可换成 c8v^8，差为 O(v^9)，不得改首退款的 (1+v)^−5。

## 4. 删除质量的 sharp 导数与共同基线

对源 i 令 f_i(m)=(R_i+m^(1/n))^n、g_i=aR_i^n。实际质量下界给
\[
 (1-A_8v)g_i\le\mu_i\le m_i^0\le g_i.
\]
写 u=m/g_i，则整个 [μ_i,g_i] 上精确有
\[
 \frac{f_i'(m)}{v/a}=(v+u^{-1/n})^{n-1}=1+O_{n,A_8}(v).
\tag{10}
\]
最后一步只在 u∈[1−A_8v,1] 使用固定光滑导数界，不需要逐盒额外近等号。因此
\[
 f_i(\mu_i)\ge A_vg_i/a
 -(1+O(v))(v/a)(g_i-m_i^0)
 -(1+O(v))(v/a)\delta_i.
\tag{11}
\]
这里将两段真实导数积分分开，**δ_(i,in) 的系数单独保留 1+O(v)**；只能对余下 𝔡_h^0 选择固定 b。不能把 δ_inside 纳入 b(v/a)(g−μ) 后仍宣称领先常数为一。

原实际乘积为 |X_i|≥f_i(μ_i)+J_i。Bonferroni 因而给 ∪X 的下界 Σf_i(μ_i)−H_8；将 (11) 插入后才扣实际重复质量及环境。由 (2)，D_h 的净系数为
\[
 [A_v-a-(1+O(v))v]/a\ge1/(2a),
\tag{12}
\]
而 𝔡_h^0≥0 仍以固定 b(v/a) 为负价格保留到望远镜。J 不重新添加，δ 不从 H 中再减一次。

当 n≥3，inside 的原 a 系数、其余删除和环境的绝对误差为 O(v^(n+6)S)=o(v^8S)，盒外 CvF=O(v^9S)。第八层无需用未知八阶常数购买这些费用。

二维使用现已正式供应的 **H56**：
\[
 \frac va\sum_{deep}\delta_{i,in}
 \le(2K_2+O(v))v^8S,\quad K_2=4/27,
\tag{H56-in}
\]
它来自实际全点原半径覆盖与正常层，而非 dyadic C=1 的改名。与 (11) 合用仍仅花 (2K_2+O(v))v^8S。R054 已审 PI8 正预留给 boundary/env=O(v^10S)，盒外=O(v^9S)，没有第二个隐藏的八阶地板。

记 κ_n=2K_2·1_(n=2)=8/27·1_(n=2)。用 S≤T/(1−a) 将这唯一正费用换成 κ_nv^8T+o(v^8S)，不要求 S/T 的未知一阶系数。因此取
\[
 Q_8(v)=Q_{8,N}(v)+\kappa_nv^8.
\tag{13}
\]
第八层误差上包为
\[
 E_8\le Q_8(v)\mathcal T+\delta v^8S+o_\delta(v^8S).
\tag{14}
\]
本式 n≥3 不需要 H56；n=2 的 H56 已由表列终审供应。恢复原免费面如有需要只允许一次，其非负收入可丢；本式不另取一份全 J 或正常禁层支付小源。

## 5. 同族前缀、八层正抵消及 δ 次序

第七前缀必须在本次八层族实际重放已审 R048 定量粗价，可取 E_7=Q_old(v)T+O(v^8S)，其中
\[
 Q_{old}=K_n[v^6(1+v)^{-4}+(10/9)v^7(1+v)^{-2}].
\]
它与 Q_8 的差为 O(v^8)，所以 E_7−Q_8T=O(v^8S)。第六前缀为 E_6=(K_n+o(1))v^6S；(4) 给 E_6−Q_8T=o(v^6S)。前五层仍只需 O(v^6S)。这些对象都不是旧第七层重命名。

各层加入同一标记半盒和内部突出回收。真实内部式保留 (1−a)，同一公共左侧为
\[
 \mathcal T\ge A_vM/a+c_{rec}P_I-C_bP_\partial-2nvS_I
 +c_{mass}\mathcal D_h/a-b(v/a)\mathfrak d_h^0-E_h,
 \quad c_{mass}=1/2.
\tag{15}
\]
δ 不改变 source/CDF 常数，所以 b 可在辅助 δ,T,m 之前固定。τ=bv/c_mass，Ω8=Σ_(j=0)^7τ^j。由 (3)，精确有
\[
 \sum_{h=1}^8\tau^{8-h}(c_{mass}f_h-bvf_{h-1})
 =c_{mass}f_8-bv\tau^7f_0.
\tag{16}
\]
初端除 aΩ8 后 O(v^9S)。先减同一个 Ω8Q_8T，而非只改末层：
\[
 \overline E-Q_8T=\Omega_8^{-1}\{
 [E_8-Q_8T]+\tau[E_7-Q_8T]+\tau^2[E_6-Q_8T]
 +\sum_{h=1}^5\tau^{8-h}[E_h-Q_8T]\}.
\tag{17}
\]
四类贡献依次为 δv^8S+o_δ(v^8S)、O(v^9S)、o(v^8S)、O(v^9S)。其中 τ²o6=o8，τO8=o8。统一图价格吸收固定边界/内部体积价，f8 和剩余正回收可丢，得到
\[
 \mathcal T\ge A_vM/a-Q_8\mathcal T
 -\delta v^8S-o_\delta(v^8S).
\tag{18}
\]
若原 δ 分配带一个固定系数，可在最初将 δ 除以这个与后选窗无关的系数。本稿仅需要 (18) 误差的 δ 系数统一有界且与新 T,m 无关。

外部指标空时所有外部源、模型和交量取零，不定义空源 CDF；原内部共享质量/突出回收账本直接给同方向。该分支不能通过给空 normal 集强加一个正成本来处理。

## 6. 移实际 M 的真实分母：候选精确系数由这里产生

令 x=M/V，S/V≤(1−a)^−1，将 (18) **直接移实际 M**，不以 aS 替换 M：
\[
 1+x\le\frac{A_v}{A_v-a(1+Q_8)}
       +C_n\delta av^8+o_\delta(av^8).
\tag{19}
\]
小 v 时分母正且有固定正下界。与 L_n=A_v/(A_v−a) 的差精确为
\[
 \frac{aA_vQ_8}{(A_v-a)(A_v-a-aQ_8)}.
\tag{20}
\]
Q_8=O(v^6)，所以把最后一个分母替为 A_v−a 的误差除 a 后为 O(aQ_8²)=o(v^8)。于是所需标量是
\[
 Q_8\frac{A_v}{(A_v-a)^2}.
\tag{21}
\]
容量项与 inside 项分别保留到这一步，不能先套用 R041 可删预算。

形式展开应为
\[
 Q_{8,N}=K_n\{v^6-\tfrac{26}{9}v^7+(40/9+c_{8,n})v^8\}+O(v^9),
\]
\[
 A_v/(A_v-a)^2=1-nv+
 [n(n+1)/2+2\mathbf1_{\{n=2\}}]v^2+o(v^2).
\tag{22}
\]
因此待全输入总装核验的上系数为
\[
 K_n[40/9+c_{8,n}+26n/9+n(n+1)/2+2\mathbf1_{\{n=2\}}]
 +\kappa_n.
\tag{23}
\]
其与
\[
 J_n=\frac{(n-1)n^3(2609+1971n+81n^2)}{26244}
       +\frac8{27}\mathbf1_{\{n=2\}},\quad
 \widehat J_n=J_n+\frac8{27}\mathbf1_{\{n=2\}}
\tag{24}
\]
的相等性现已由 Terra 的最终成功实跑 `r057-computation/settlement.wl/log` 核为残差零。**式 (23) 来自真实分母和一次已审 inside 费用，不是预先把 J_n 当上系数再机械平移。** 一般 n≥3 的系数核为 J_n；二维正常分母给 15694/6561，再加独立 inside 8/27 后给 17638/6561。程序只验证这里已经推得的标量，不供应集合几何。

通式核的合法分支如下：先令 a=v²w、0≤w≤1，归一化 (20) 后的前三系数为 K_n、−(n+26/9)K_n、J_(n,0)+2K_nw，其中 J_(n,0) 是 (24) 不含二维指示项的多项式。分母在 v=0 统一远离零，固定 n 的三阶 Taylor 余项对 w∈[0,1] 一致为 O_n(v³)。n≥3 时再代 w=v^(n−2)，新增 v²w=v^n=o(v²)；n=2 时 w=1。这样得到全体整数维数的解析结论，不是仅以 w→0 代替一条未知的非一致展开。

## 7. 全部输入上确界与匹配下界

若目标上极限≤widehat J_n 失败，只能选一个固定 η>0 和 v_j→0，使全输入归一化余量大于 widehat J_n+η；**不要求该序列趋于无穷**。按 C_n 的上确界定义取实际有限输入，再以 o(Vav^8) 精度取外部严格 halo 紧集。密度定理保证输入几乎处处在自身 halo；外部 halo 可测有限，内正则性适用。

得到同一实际序列满足 (1)，例如固定 B=1 即可。先固定这个序列与 R049 构造，再给定足够小的固定 δ，使 (19) 的 C_nδ<η/4；依次选 T,m、soft 常数，最后 v→0。逐序列 little-o 模量可不同，但对该选定序列确趋零。矛盾后才令 δ→0。没有把未证统一速率移出上确界，也没有重新选择 I 或 CDF。

已审 R041 对 n≥3 给实际有限输入 lower P7+J_n av^8+O(av^9)。已审 R054 二维两尺度输入给 lower P7+(J_2+8/27)av^8−O(av^9)。因此本稿的待审精确八阶候选定理为
\[
 C_n(1-t)=P_{7,n}(t)+\widehat J_nav^8+o_n(av^8).
\tag{25, candidate}
\]
现有前提均已供应，但 (25) 尚待本稿二维补集及全局总装的非作者严格复核。不得因已有前提和标量核通过，就跳过新总装审查或宣称已发布精确八阶。

## 8. ordinary b8：完整基线第八项不能省

在 (25) 成立的相应维数范围内，令 0≤j≤8，
\[
 \beta_{n,j}=2^{n+j}\left[
 \sum_{q=1}^{1+\lfloor j/n\rfloor}
 (-1)^{j-(q-1)n}\binom{n+j-1}{j-(q-1)n}
 +\mathbf1_{\{j=6\}}K_n
 -\mathbf1_{\{j=7\}}(n+26/9)K_n
 +\mathbf1_{\{j=8\}}\widehat J_n\right].
\tag{26}
\]
则 C_n(1−t)=1+Σ_(j=0)^8β_(n,j)t^(1+j/n)+o_n(t^(1+8/n))。解析 L_n 在 v^(n+8) 的 coefficient 仍包括 q≤1+floor(8/n) 的全部碰撞，不能把 ordinary β8 直接写成 2^(n+8)widehat J_n，也不预设 β8 的符号。j≤7 与旧 b 完全相同。Terra 的 WL 与独立 Fraction 新核均给二维 baseline 第八项为 256，故 β_(2,8)=1737988096/6561；这仍以 (25) 的总装成立为真实系数解释的前提。

第八项的 q 上限分别为：n=2 时 5，n=3,4 时 3，n=5,6,7,8 时 2，n≥9 时 1。因此包括 n=8 的整数幂碰撞；不能沿用第七项的“n≥8 只有一项”分界。

## 9. 全部 m 的精确八阶传递：五阶 σ 才足够

模型仍为共享同一半径、公共交集正有限测度归一化的精确约化。记
\[
 s=\sigma_m(\varepsilon)=m[1-(1-\varepsilon)^{1/m}],\quad
 \varepsilon\le s\le-\log(1-\varepsilon)\le2\varepsilon
 \quad(0<\varepsilon\le1/2).
\]
s≤log2<1，所有 m 同时在有限支路。

旧 A4 的 bracket 只到 ε³，导致 σ−εA4=O(ε^5)；二维 p8=5 时不够 o(ε^5)。定义五阶参数多项式的 bracket
\[
 A_{5,m}(\varepsilon)=A_{4,m}(\varepsilon)+d_m\varepsilon^4,
 \quad d_m=\frac{(m-1)(2m-1)(3m-1)(4m-1)}{120m^4}.
\tag{27}
\]
这里 εA5 本身到 ε^5，不是把 bracket 的阶数与 σ 的阶数混同。置 q=1/m∈(0,1]，第六导数是
\[
 \sigma_m^{(6)}(z)=\prod_{r=1}^5(r-q)(1-z)^{q-6}.
\]
在 0≤z≤ε≤1/2 上介于 0 与 5!·2^6；Taylor 给全 m 一致的
\[
 |\sigma_m(\varepsilon)-\varepsilon A_{5,m}(\varepsilon)|
 \le(32/3)\varepsilon^6.
\tag{28}
\]
m=1 高阶项及余项均为零。Taylor 余项非负也给 ε≤εA5≤σ≤2ε。

设 p8=1+8/n≤5。标量同一个 r_n(t)=o_n(t^p8) 的模量
ω_n(δ)=sup_(0<t≤δ)|r_n(t)|/t^p8→0，所以
\[
 \sup_{m\ge1}|r_n(\sigma_m(\varepsilon))|/\varepsilon^{p8}
 \le2^{p8}\omega_n(2\varepsilon)\to0.
\tag{29}
\]
对九个固定 p=1+j/n∈[1,5]，均值定理和 (28) 给
\[
 s^p-[\varepsilon A_{5,m}(\varepsilon)]^p
 =O_p(\varepsilon^{p+5}).
\]
有限求和误差 O_n(ε^6)=o_n(ε^p8)，最紧二维仍有一阶余量。因此 (25) 在其成立维数范围内推出全 m 一致的
\[
 C_n^{(m)}(1-\varepsilon)=1+\sum_{j=0}^8\beta_{n,j}
 [\varepsilon A_{5,m}(\varepsilon)]^{1+j/n}
 +o_n(\varepsilon^{1+8/n}).
\tag{30, conditional}
\]
自然参数则不作截断，a_s=2^ns、v_s=2s^(1/n)，直接由同一个标量模量得到
\[
 \sup_{m\ge1}\left|
 \frac{C_n^{(m)}(1-\varepsilon)-P_{7,n}(s)}{a_sv_s^8}
 -\widehat J_n\right|\to0.
\tag{31, conditional}
\]
这不是交换逐 m 的极限与 sup；有效阈值来自同一个标量小区间，均与 m 无关。固定维数常数不推广为维数一致。

## 10. 实际计算材料与边界

计算由已明确委托的 Terra 在本机实际执行，本作者实际读取最终代码和输出，没有重跑旧批次。必要材料保存在

`/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/r057-computation/`。

| 文件 | 实际 SHA256 |
|---|---|
| settlement.wl | `505dd6488197f3dd1a7fb463ad18c5306c8136c20caf0ede3384c99699be0f1d` |
| settlement.log | `520489ad350b5d580b261b8b2d405da2b72af539163f263e4e93c3f2546be8f4` |
| fraction_check.py | `0482fcf2efee7d08335ddf2623340294c566aa4333ace8a14350d6c616703841` |
| fraction_check.json | `793244cfb3221b56578b350e3bac2cfcc62947547acffd1b425ff408e7485263` |

最终 WL 日志为 Mathematica/Wolfram 15.0 ARM、PID 21045：真实比值差分解残差 0；一般 n≥3 三系数差 {0,0,0}；二维正常系数 {4/27,−176/243,15694/6561}；加 inside 后 17638/6561。第八层精确源缺量 telescope 残差 0，减 Ω8QT 的精确恒等式残差 0，四类误差系数次数 {8,9,8,9}。第五参数系数残差 0，六阶导数乘积超过 120 的违反集合 False，统一 Taylor 余项常数 32/3。fresh 导数的准确归一化代换残差 0。

日志中较早一行 `degree check ... {8,9,8,9,9}` 仅为次数算术，不是望远镜恒等式的机器证明；其后新增的两条残差 0 才是准确恒等式核验。第六前缀的 e6 在本稿解析论证中必须趋零，不能从数字次数 8 推断 little-o。

独立 Fraction 程序直接将二维实际比值展开至 v^10（这里 av^8=v^10），给 normal_J2=15694/6561、inside_J2=17638/6561、baseline_j8=256、ordinary_b28=1737988096/6561、all_exact=true。它不是从显示小数重构有理数。q=0、1/2、1 的 σ 第五项分别为 1/5、7/128、0，作为通式边界的准确抽查。

首次 WL 日志有未完全化简表达和错误的导数代换，首次 Fraction 日志误抽 v^8 而非 v^10；这些保留在 first-unresolved/first-invalid 文件中，不作为通过证据。本表绑定的最终成功代码/输出已修正并实跑，不能把首次错误日志的退出状态当成认证。

复现分别为 `wolframscript -file settlement.wl`、`python fraction_check.py`。这些程序只核标量代数、有限系数与参数阶数；不承担 R053 实际帽积分、R055/R056 的集合几何、§2 的二维补集、同族构造、预算不重用或全输入上确界的证明。

## 11. 候选范围及待审问题

本稿 §2 二维补集重放和 §§3–7 全局预算是本轮需要非作者完整审查的新内容；§§8–9 是 (25) 下的系数及全部 m 传递。现有 R053、R055、H56 等先验已按实际版本供应，本稿没有剩余的“等待某个常数一引理”条件，但新总装仍处于候选状态。

应重点核对：全部对象确为同 N8；完整 J 一次、D_N 一次、inside sharp 一次；δ 的选择不改源族与 CDF；二维容量/长条有自身实际桥；移实际 M 以后才得到 (23)，不把安全可删预算冒充全输入必要性。只有两份非作者总装复核通过，才按实际完成范围对外发布 (25)。本稿不改论文、不提交仓库，也不把机械展开单独记为已完成的一轮进展。
