# 杀死过程、Duhamel 与首次跳跃标签分解：方法与验收边界

1. 明确杀死是独立指数时钟、状态依赖杀死率还是吸收跳转；写出含墓地状态的生成元与存活半群。

2. 对时间齐次 càdlàg 强 Markov 过程定义首跳时 T、跳后位置 Y 和次概率核 F_x(ds,dy)=P_x(T∈ds,Y∈dy,T<ζ)。有界可测 f 满足
   E_x[f(X_t);t<ζ]=E_x[f(X_t);T>t,t<ζ]+∫_[0,t]F_x(ds,dy)P_{t−s}f(y),
   其中 P_t f(y)=E_y[f(X_t);t<ζ]；t 处可能有原子，故闭区间不可随意改成开区间。只有另外建立 hazard 结构后，才能把 F 拆成生存因子、强度和跳核。

3. 对共同定义域的 C₀ 半群生成元 A,B，若 V=A−B 延拓为有界算子，则
   S^A_t−S^B_t=∫₀ᵗ S^A_{t−s}V S^B_s ds.
   在共同稠密定义域上对 S^A_{t−s}S^B_s 求导，随后用强连续性、V 有界和稠密性延拓。不同定义域或无界生成元差不在本结论范围。

4. 若互斥可测事件 A₁,…,Aₙ 覆盖事件 A，且 Z 可积，则 E[Z1_A]=ΣᵢE[Z1_{Aᵢ}]；把 Aᶜ 的零跳/余项单独保留。可数划分另需绝对可积或可交换求和的条件。

5. 若引用 P-c580aef7365f8786 Corollary 5.5，仅有每个有限λ的域外占用惩罚表示；λ→∞到首次退出杀死不是本篇已证明的极限。另核验边界正则、路径占用与首次退出事件是否等价，以及所需极限交换。

6. Miles–Keener P-ee636b8a2c2329ad 的式(3)/Theorem 1针对一维扩散的 absorbing pre-jump 密度 p：∂ₜp=Lp−λp（L是 forward operator），且 p_ℓ(x)=∫₀∞λ(x)p(x,t)dt、p_τ(t)=∫_Rλ(x)p(x,t)dx。原文假设足够正则、无爆炸且强度有限，并在证明中用 p→0 于无穷远及积分交换。迁移时另核验边界通量、域和可积性；只有 P(τ<∞)=1 才归一化为概率。若 reset kernel 为 J，landing law 是 pre-jump measure 的 J-推前，须分开报告。

## 首跳 hazard 来源和 reset 记号

Miles–Keener 的实际符号、事件密度、条件和源文页码见 J02 claims。对于已定义的 reset Markov kernel J，pre-jump measure μ 与 landing measure ν 的关系 ν(B)=∫μ(dx)J(x,B) 在 [J02 step 6 local derivation](../../../audit_current/jk/derivations-j02-step06.md) 中按测度定义核对。

## 受限恒等式的证明记录

首跳分解、半群 Duhamel 恒等式和事件标签分割的具体假设与逐步证明见 [J02 steps 2–4 local derivation](../../../audit_current/jk/derivations-j02-steps-02-to-04.md)。这些结果仅在各自列出的条件下成立，不补足目标模型缺失的强 Markov 性、停止时、生成元共同定义域或可积性。

## 不可省略的限制

杀死、越界跳跃和普通转移可能是不同机制；遗漏墓地质量或把边界杀死当一次存活跳跃会破坏质量守恒。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。

局部证据补充见 [SOL_HN记录](audit-sol_hn.md)；其审读者与范围以记录为准。
