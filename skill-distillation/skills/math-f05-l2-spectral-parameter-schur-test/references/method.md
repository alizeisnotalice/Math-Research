# L² 谱参数与 Schur 估计：方法与验收边界

1. 在有限支撑稠密子空间定义A(λ)，写明两Hilbert空间、底测度及λ的含义；矩阵行列Schur、逐项Schur乘子和谱分解是不同命题。

2. 取正权p_i,q_j，验证Σ_j|a_ij(λ)|q_j≤R p_i及Σ_i|a_ij(λ)|p_i≤C q_j，两侧同时成立且R,C对所需λ范围统一。

3. 已检查Cauchy–Schwarz构件：把每行Σ_j|a_ij||x_j|拆为√(|a_ij|q_j)与√(|a_ij|/q_j)|x_j|；平方后求i、交换非负和并用列界，得||Ax||²≤RC||x||²。

4. 积分核同理先检可测性、正权、Tonelli与有限支撑估计，再密度延拓；权可依赖λ但界及定义域必须可控。

5. P-3da1b59567d35c36的Luna全文卡现覆盖99页；第3节行列和界为√(RC)，第4节正PSD逐项乘子另为||H∘A||≤sup h_ii||A||，不能相互替换。全文覆盖不等于全部外引定理独立认证。

6. 若实际问题是逐项乘子且H=LM*，给D_L=sup_iΣ|l_ir|²与D_M=sup_jΣ|m_jr|²，才能用√(D_LD_M)界；Bennett逆向、Peller与Birman–Solomyak双谱积分是作者转述，外部原始证明未另核。

7. survey双谱积分入口另需正交谱测度E,F和核h(mu,lambda)=∫m(mu,x)l(lambda,x)dx，两因子L2 norm的esssup受控；仅矩阵行列界或有界h都不能替代该因子证书。Birman–Solomyak/Peller等原始外引证明仍未另读认证。

8. 参数端点、极点、零权及正交谱块分别检查；常数1的有限Weyl谱/奇异值幂和比较也不提供具体A(λ)一致界。范数统一界不证明核逐点收敛、Hilbert矩阵最优权或任何私人谱接口。

## 不可省略的限制

两侧正权行列界缺一不能调用Schur；survey全文卡已覆盖99页，但不含私有参数λ一致定理，且Bennett/Peller/Birman–Solomyak为外引依赖。独立Cauchy–Schwarz只认证该弱构件，不能代替真实核、谱重叠、端点或private接口。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。

局部证据补充见 [ROOT记录](audit-root.md)；其审读者与范围以记录为准。
