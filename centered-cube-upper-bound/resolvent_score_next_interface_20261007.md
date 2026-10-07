# 下一阶段候选：原预解核的 Gamma 来源 score

2026-10-07。仅登记解析推导候选，尚未独立审计或安排新原核数值；不登记为主费用。不新增或重复运行已有实验。

原 G=(I+cB)^−1=∫e^−t P_ct dt，G²=∫t e^−t P_ct dt。利用同一 Exp(1) 隐变量 T 和同一条件位移，可由条件 Jensen 推出 χ²(G²||G)≤E(T−1)²=1。允许条件位移含奇异部分，用测度 RN 导数；不把它误当实际后验独立。

p_r=(1−r)δ+rG，q_r=G*p_r=(1−r)G+rG²。可在带标记的潜变量空间定义基准p的 holding标签概率1−r、Exp标签概率r：目标q的RN密度为 holding标签0、Exp标签 (1−r)/r+T。其均值1，方差1/r+r−1。把标签送到实际位移再用条件Jensen，得到实际χ²(q_r||p_r)不超过此值。坐标乘积下中心化score独立，预期给

||S K_r||TV ≤ sqrt[n(1/r+r−1)], 0<r≤1.

同一固定 H_r 的实际Poisson跳跃表示，总数J~Poisson(nr)，SH_r 的latent signed score n−J/r；因此 ||SH_r||TV≤sqrt(n/r)。S是nI−ΣGi。

这些仅是单参数卷积核TV候选。取sup_r以后没有免费的同阶L1或弱型结论；分别三角估计还丢失D=K−H的取消，尤其小r。下一阶段应先独立核代数及原核守卫，再决定是否值得研究差项score的联合耦合，不能据本候选宣称geom已付。
