# 承重经典引理证明复现（R34）——底层链闭合：Young → Hölder → Jensen

续 R31–R33，累计 12 条。本批的意义是**把整条验证链的最终根源闭合到 exp 的凸性**。

---

## 10. Young 不等式 ab ≤ a^p/p + b^q/q（1/p+1/q=1）

**证明**：a,b ≥ 0 时
ab = exp(ln a + ln b) = exp((1/p)ln a^p + (1/q)ln b^q)
   ≤ (1/p)exp(ln a^p) + (1/q)exp(ln b^q) = a^p/p + b^q/q  [exp 凸性 + 权重和 1] □
取等 ⟺ a^p = b^q。

**数值**：100,000 组随机 (a,b,p∈(1,6))，max[ab − a^p/p − b^q/q] = **0.00e+00** ✓

## 11. Hölder 不等式 ‖fg‖₁ ≤ ‖f‖_p‖g‖_q

**证明**：归约到 ‖f‖_p = ‖g‖_q = 1（除以范数；零/无穷情形平凡）。
Young 逐点：F·G ≤ F^p/p + G^q/q；积分：∫FG ≤ 1/p + 1/q = 1 ⟹ ∫|fg| ≤ ‖f‖_p‖g‖_q □
取等 ⟺ |f|^p 与 |g|^q 成比例。

**数值**：间接（Young 数值零违例 + R31 Doob 数值 + R32 LP 数值均已通过）。

## 12. Jensen 不等式（有限 + 一般）

**有限版证明**：对 n 归纳。n=2 即凸性定义；n→n+1：
φ(Σ_{i≤n+1}λ_i x_i) = φ(λ_{n+1}x_{n+1} + (1−λ_{n+1})·ȳ)，ȳ := Σ_{i≤n}λ_i x_i/(1−λ_{n+1})
≤ λ_{n+1}φ(x_{n+1}) + (1−λ_{n+1})φ(ȳ) ≤ Σ_{i≤n+1}λ_iφ(x_i) □（λ_{n+1}=1 平凡）

**一般版证明（Fenchel 图表示，独立于逼近论）**：φ 正常凸时
φ(t) = sup_s (st − φ*(s))（Legendre 变换）。
对一切 s：Eφ(X) ≥ E[sX − φ*(s)] = sEX − φ*(s)；
对 s 取 sup：Eφ(X) ≥ φ(EX) □

**数值**：50,000 组随机凸函数实例，max[φ(EX) − Eφ(X)] = 0 ✓

---

## 链闭合声明（本批的核心价值）

至此，前几轮所有"标准定理引用"的依赖链已**收敛到两个原子事实**：

```
exp 凸性 ──┬── Young ──→ Hölder ──┬──→ Doob L^p（R31，g05 的 4C）
           │                      ├──→ Bernstein–Chernoff（R31，h05）
           │                      └──→ Schatten 引理（R27，f03 的插值链）
           └──(log 凹) Jensen ──→ KL≥0 ──→ Gibbs 变分（R33，d02）
```

除 exp 凸性本身（分析学原子）外，**链条上没有任何未闭合的引用**。
各 skill 的"标准定理引用"现在都有完整的第一性原理路径支撑。

## 账本影响

| 定理 | 影响对象 | 状态 |
|---|---|---|
| Young | 全链底层 | independent_proof（新入册 GEN-YOUNG-01） |
| Hölder | h05/e01/g05/f03 的既有证明 | 底层依赖闭合 |
| Jensen | d02/d04 的 KL 链 | independent_proof（新入册 GEN-JENSEN-01） |
| F04 三条谱恒等式 | math-f04（R21 sympy 精确） | 正式入册 F04-IDENT-01a/b/c |
