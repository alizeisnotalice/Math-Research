# 卡↔PDF 保真性验证（R47）——P-2f81500e5dfa70c6 通过（含 2 处条件表述偏差已修订）；重依赖卡覆盖完成

---

## 1. P-2f81500e5dfa70c6 / Oertel《On Jump Measures of Optional Processes with Regulated Trajectories》（arXiv:1508.01973v1，11 页）

依赖度最高的未核验卡：被 3 个 skill 依赖（i02 / j01 / j05）。

**逐字比对结果**（PDF p4–p5、p9–p10）：

- **Proposition 1（PDF p4）**：迭代枚举 T^A_1=inf{t>0:ΔX_t∈A}、T^A_n=inf{t>T^A_{n−1}:…}、
  evanescent set 上严格增停时、disjoint graphs ·[[S^A_n]] 并集恰枚举 {ΔX∈A}——
  与卡陈述**逐字一致** ✓
- **Theorem 4（PDF p9）**：optional + regulated + ΔX_0=0 ⇒ ΔX optional；每条轨道
  jump set 非有限 ⇒ 存在严格增停时序列，disjoint graphs 枚举 {ΔX≠0}（thin set）——
  与卡**逐字一致** ✓；证明四步（X− predictable/X+ adapted → D_m 有限分解排序 →
  début of optional set + Theorem 3 → relabeling）与卡 proof_steps 吻合 ✓
- **Theorem 5 + Corollary 1（PDF p10）**：j_X 为 integer-valued random measure、
  评价 = 停时上可数和——与卡**逐字一致** ✓
- 卡的 gap 记录（"pathwise infinite-jump condition excludes paths with only finitely
  many jumps…not intensities or compensation identities"）与原文 Problem 1/2 的
  开放问题分区正确 ✓

## 2. 发现并修订：2 处条件表述夹带非原文内容（与 E03 同类，条件层）

1. **Prop 1 assumptions**："0∉A, hence is bounded away from zero"——PDF 仅假设 0∉A；
   "hence" 推断不成立（A={1/n}）。原文证明 lim∈Ā⇒ΔX_{t*}≠0 一步需要 0∉closure(A)，
   系**原文证明 gap**；卡隐式修补但未标注。
2. **定义卡 jump measure**："A bounded away from zero" 的有限性条件（Lemma 2）——
   原文仅写 A∈B*；0∉A 不足以保证 N^A_X(t)<∞（A={1/n} 在 accumulating jumps 下
   可数无穷多 A-jumps）。卡的修补数学上必要，但原文 p5 "A⊆R\(−ε,ε) for all ε>0"
   量词系原文疏漏（应为"存在 ε"），卡应标注为卡的澄清而非原文。

**处置**（append-only）：三份卡副本（i02/j01/j05，sha256[:16]=402d5374d5e04090 一致）
原文保留至 `*_original_r47` 字段，新表述含【修订 R47】内联标注；三个 skill 的
verification-log 各加更正条目；SKILL.md 正文未转录这些假设，无需修订。

---

## 卡↔PDF 累计（R36–R47 终版）

```
十二轮：26 张卡核验 —— 陈述主体 26/26 无虚构定理；
条件/常数层偏差 3 处（E03 势公式、本卡 2 处），均已修订并留痕
重依赖卡（in-degree≥3 全部 + in-degree≥2 抽样）覆盖完成
```
