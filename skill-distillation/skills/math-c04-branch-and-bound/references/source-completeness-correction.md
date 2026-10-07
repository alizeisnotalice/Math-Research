# 来源完整性更正与有效数学范围

原8页PDF物理完整。PDF第8页左栏末尾的“simply set”在右栏顶部续接；MD抽取将续文放在357–358行而非文件末尾。旧卡的“原附录缺页/截断”事实理由撤回，历史卡保留原字节。

原Lemma3在PDF第7页、MD342–343行严格要求a≤c<d≤b；证明跨PDF7–8页，MD344–352、404–422和右栏续文357–358行。两组独立审核已确认原符号、全部允许端点及局部恒等式；c=d反例仅排除错误扩域。

局部证书可独立核对：置z=x−a、C=c−a≥0、L=d−c>0、R=b−d≥0、B=C+L+R、s=√[CR(C+L)(L+R)]。取α=(BL+2CR+2s)/L²、β=α−1≥0。C=R=0时α=1、β=0、g=h；其余情形β>0，取γ_z=[B−α(2C+L)]/(2β)、γ=γ_z−a。g−αh关于z的判别式L²α²−2(BL+2CR)α+B²=0，故恒等于β(z+γ_z)²。此证明只使用L>0；不宣称L→0系数一致有界。

有效状态是“原PDF完整；Lemma3严格域已局部核验；完整SOS层/定理1及外引证明链仍待核”。Putinar完整假设、矩/SOS对偶、有效阶数和数值证书仍需独立证据。本补正不新增全篇阅读、不认证全SOS或中心立方体课题结论。qBnB(3)的独立操作门另要求L₃>0、r>0、可逆正则化Hessian及迭代域；对应P-922909d5e8e47a30 PDF9–11页、MD471–558行，作者的L₃≥0字面端点不能直接进入Newton除法/逆矩阵。

旧sidecar：`evidence/reviews/TEAM_CD/behavior/SOS-appendix-completeness-retraction.json`，SHA-256 `b265e1e9b817a25d3524d70432dbd32dc653aa1e6bdcc0eb404549e2efe598c1`；其中pending审核状态由以下正式receipt替代，只覆盖局部Lemma3。

局部补正卡：`evidence/reviews/TEAM_CD/evidence/SC-P-deb72842d1067b4f-appendix.json`，SHA-256 `f4be8f23e416f334d997266d34dee82ea32d523538ea73842ca50fdc7ac6d7ae`。

独立receipt：`evidence/reviews/AUDIT_CD/followups/P-deb72842d1067b4f-lemma3-correction.json`，SHA-256 `96e0c5930902520af50c1f3d1ed46f29766062699eb53fa94e5b7cf9dbc48410`。

独立receipt：`evidence/reviews/AUDIT_AB/round2-cd/P-deb-scoped-correction.json`，SHA-256 `77e9f8bd9edba61e4b1e606530f34770abad1ffcd8b00093e4d3348bb59b22fc`。

原PDF SHA-256：`deb72842d1067b4fa98ad54677019c6a61e469f3b205f7fb4839659b59e447d1`。
