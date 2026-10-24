# 第47轮运行前独立源码gate

全文静态审查结论GO，无阻断；未运行模型、点、LP或证书。

- `search.py` SHA `715678d73fa80e7eb108d9642a6b02fde242773be7feec8d44931359510731d2`。
- `solver.m` SHA `d37a5f3194e8fce144f7eef94a2c2f84de03579180a68769abecf1c161ddcd8c`。
- `run_once.py` SHA `2622068da96b9d24df5f42c38e03093504bf1bcf219313d66533150c176808a8`。

- 必要模型严格删去36720充分Bernstein行和所有halfslack；51574旧点原目标+P1/P2两A行+24abs+1结构，共51601×26。缓存只用于一次数值发现，未给float伪造误差界。
- Ψ/S1复用baseline Spec但分别调用**未乘c**的Psi/APsi_point或S/AS；12T独立Spec，192bit在构造常量之前设置，D4质量、h及固定12列来源正确。空aSpec已避开，无dummy系数。
- 内/外AK仅用pin core闭式；AK正坐标且r≠R，外区J活跃片断言0<y/(x+y)<R/(x+y)<1及x−(R−y)>0。非finite/不支持停，不跳行、不换支持，无新增积分路线；K点排除原点。
- 全局(种类,u,v)缓存14列，392上限含28初始列及top26支撑重评；调用前后账本/trace与源hash齐全，不重复同点列。结构/abs精确不花点评价。
- MATLAB一次20秒soft、第五输出lambda；从有限lambda.ineqlin准确二进制有理量化2^-40，非负截取top26（并列原行号）后固定，不重选。该y不假称准确stationarity，最终残差全额补偿。
- 支撑点行A端点由F区间反号交换，z列零；b为原目标负号，结构/abs准确。Ahat/radius、ei=ΣδMj、relaxed_b=b+ei、r=f_lower+Ahat^T y、L=−y·relaxed_b+Σmin(rj*lj,rj*uj)逐项方向正确。
- 保存40次根夹逼重核后，ψ目标系数取α/(2−α)乘root/scale，不含候选c；c盒非负使其为成本严格lower。S1/T用准确半质量且保留符号，原26维准确盒未换float。
- 300秒自有组hard含导入/LP/点/保存，首异常/超时不重试，无核候选或全域checker。独立保存值核验另30秒，仅重组端点/来源/λ/残差/根，不调用生成代码或数学原语。
- 结论范围严格为固定α/R、12方向、参数盒及c≥s*c1子族；即使L>17/5也不是C2下界、全核族不可达或新上界。L≤17/5/不闭合不证明存在3.4可行核。
