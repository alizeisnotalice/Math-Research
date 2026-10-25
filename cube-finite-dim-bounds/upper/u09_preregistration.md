# U09预登记：体对角finite-frame迁移与真实势正部截断

2026-10-25，基于主已推第8轮4ea9d04f3f8b77a7149e53d06ab1463899fddfe3。已读主009登记、snapshot、frame接口和外部设计，准确输入bytes/hash在u09_inputs.json。只原clone upper/u09_*，不改旧源、bridge、总账，不提交推送扩招。不是唯一工作者，不撤销他人改动。主负责g/J成本，bridge负责p/G准确锥式；本批不重复。

固定V_n={(1,±1,...,±1)}, n3/4，L_V=−Σ|v·D|，真实保守全空间Cauchy和，Lebesgue、同一共同尺度r。先完成一个finite-spanning-frame迁移定理，再一个K_c=(G−c)_+命题；完整G由正密度积分定义，不假定锥和每项正或在frame lines有限。

证明预算：只此固定frame。迁移须逐项列transient上界、operator域H1、u L1能量、a.e饱和、质量cutoff、非局部接触、核原点cutoff/CND/无限源、同D全r和一般L1。截断须真实付G局部可积/weak-L^(n/(n−1))、正部L1、mollified正部的全空间分布极限及原点原子/可能无限线奇性。若失败给准确反例/条件；不能用off-line点值代替分布符号。

拟核严格路径：heat kernel min(L∞,Ct^-n)→G weak-Lq→(G−c)+ L1；G_ε=G*ρ_ε及凸近似的非局部Jensen→L(G−c)+≥−δ0，若成立再用大球cutoff计正源质量。line分类单独检验：one-large-jump下界，off all frame lines的two-large-jump上界；不预设log/dist精确渐近。

实际Python有限诊断：仅一个5节点对称正率path矩阵，准确Dirichlet Green零延拓及正部截断的Jensen/源质量（不模拟R^n）；SymPy只核n3/4 weak-Lq层积分、径向mass公式与c=g/n代数。MATLAB给同一有限矩阵和质量代数脚本，由主串行；无clear/path。不得加matrix、参数grid、成本最小搜索或alpha变化。

停止：一个迁移定理和一个截断命题或严格反例及上述有限诊断完成后停止。实际低于8/16必须另有主g/J准确证书并送原audit；本分支不据构核接口或浮点比值宣称新U。

推导中简化（仍同一命题/批预算）：真实potential有P_tG=G−∫_0^t p_s（a.e.有限点），正部1-Lipschitz/Jensen直接给P_tKc−Kc≥−∫_0^t p_s。主证明优先用此与L∞测试生成元极限，避免把mollifiedG的H1当额外假设；mollifier路线仅作等价核查。此项无新增计算或trial。线分类附于主结论后，不阻塞源合法性。
