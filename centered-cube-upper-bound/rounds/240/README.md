# 第240轮：距离加权的全尺度来源预算

日期：2026-10-17。一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。

真实actual-soft重叠满足I0(y,z)≤a^-n exp(−E_z(y))O_R(z)，E_z=Σ(|yi−zi|/b−1)_+。同一来源测度上的FIRST与点态包络联合Hölder，将粗质量β(D)替换为显式m_w=∫S exp(−E/2)dβ，保持原fullν后验。全部j/M和外系数为84Λe⁵∫δ[R_w+.5log(1/δ)]³dν，R_w=1+log+ρ_w+log+K+logB。

一次远门F_z={E≥2log+ρ_D}满足ρ_w,F≤1，其费用为84Λe⁵[1+log+K+logB]³||f||₁（原K跨外资料共同固定时）。保留K：尚未认证原合法K(n)增长，不能直接称无条件根阶付款。远门可能为空；近门与一般补集仍未支付。辅助附件给出真实共同soft归一空间权，所有菜单、积分及n费用显式。

## 冻结归档与实际验证

[research.zip](research.zip) 保存11个必要成员：完整LaTeX、主路线/核工具两附件、独审、合同、报告、数值代码/结果/说明、复现driver和manifest。不编译PDF。

- ZIP：39,815 bytes；SHA256 `6faa047b672964969b93ff8532ece18a4e875fe6a649b2b25d18673824c8a5ba`。
- 先独筛命题后审证明，强O因子、端点、完整节点常数、一次远门及原K资格均在报告范围核验通过。
- 本机32个正式配置：13实际核积分、10共同Holder数组、9远门代数；另实际运行60项真实核/CDF/tensor检查。
- mp30/50核积分对比最大差2e−32，是普通高精证据，不是区间证书；全部数值不认证FULL历史或一般维度阶。
- 仓库外完整driver复跑成功，38,760-byte主JSON逐字节一致，SHA256 `99a0b9cac34154af4757aae67903571bde737af8cf16e581eaea2dfd301650d5`；10项manifest在复跑前后通过。

解压后以Python3.12及mpmath1.3.0执行 `python3 reproduce.py --output-dir <包外不存在的新目录>`；本机依赖路径和参数见numeric_notes.md。冻结目录 `/Users/zhengzhihao/Desktop/T/output/centered_cube_round240_20261017_final_checkpoint`；完整回执 `/Users/zhengzhihao/Desktop/T/output/centered_cube_round240_20261017_final_replay/receipt.json`。

本轮仅新增README和ZIP两个tracked文件，另更新项目总览；临时work在仓库外备份与逐成员核验后移除。下一瓶颈是近来源加权质量与共同soft/真实接收几何的联合预算，另保原K及所有分支、极限费用。
