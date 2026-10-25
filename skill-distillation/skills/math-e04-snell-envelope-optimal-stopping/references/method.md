# 最优停止与 Snell 包络：方法与验收边界

1. 离散有限时域且X_t适应可积时，从S_T=X_T逆推S_t=max(X_t,E[S_{t+1}|F_t])；以条件期望检超鞅与支配。

2. 对任意支配X的可积超鞅Y逆推得Y_t≥S_t，证明最小性；τ*=inf{t:S_t=X_t}为停止时，有限停止后的S保持鞅给E X_{τ*}=E S_0。

3. 随机化停止另给H适应右连续非降G、G(0)、终端总质量和未停止收益；质量≤1是subprobability，不能默认与强制τ≤T同一问题。

4. P-faa1738580ee2154 Theorem3.1按印刷GH仅G(T)≤1。有限T,k≡−1时普通停止值−1而G≡0值0，故全符号收益等价有反例；须修正为总质量1并处理t=0，或显式加入不停止零收益选项后另证。

5. 广义逆α(r)=inf{s:G(s)≥r}的停止事件≤/≥在右连续非降G下有效；打印G(α(r))=r a.a.在jump失败。Stieltjes换元单独由α#Leb_(0,G(T)] =dG证，保初始/终点质量及可积性；该local repair不修负收益mass gap或逐分位最优对应。另，PDF p.9 Theorem 4.2(a)打印的G*(t)=1_{t≥τ*>0}+1_{τ*=0}若P(τ*=0)>0，则G*(0)=1，违反Problem 2.3的G(0)=0；值逼近/相同不等于精确optimizer map可行。

6. Theorem 4.2(d)的singular-control jump inverse必须明确dG=e^ξdξ及Lebesgue–Stieltjes integrand的前/后跳约定。在按Problem 2.4字面post-jump e^{-ξ(s)}解释时，单位G原子需要1=ΔG=z e^{-(a+z)}，有限z>0不可能（z e^{-z}≤1/e）；但原文未明确此跳值约定，所以仅保留该literal-interpretation条件结论，并将jump optimizer correspondence标为convention-dependent unresolved。

7. 连续时间Snell需右连续、class-D等适用假设；一般信息流论文没有证明经典Snell最小超鞅定理，不能充当该定理来源。
## 不可省略的限制

有限离散Snell是已检查构件。一般信息流来源允许随机化剩余质量，有限T负奖励给其值等价反例；缺终端质量/未停收益约定时不引用该等价。连续时间和无界停止需独立正则性及可选抽样条件。

# 有限离散Snell包络：已检查基础推导

有限时域0,…,T，适应奖励Z_t可积。令Y_T=Z_T，
Y_t=max(Z_t,E[Y_{t+1}|F_t])。逆向归纳得Y适应、可积且支配Z，
并满足Y_t≥E[Y_{t+1}|F_t]，所以是超鞅。
任何支配Z的可积超鞅U，由U_T≥Z_T和逆向归纳得到U_t≥Y_t。
令τ*=inf{t:Y_t=Z_t}（T必可停止）。在τ*前Y_t=E[Y_{t+1}|F_t]，
停止的Y是鞅，故E Z_{τ*}=E Y_0；任意τ≤T有E Z_τ≤E Y_0。
这里不覆盖无限时域、连续时间或非可积奖励。
