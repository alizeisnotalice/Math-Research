# 独立复核记录（append-only，2026-10-06）

本文件由外部审核批次写入，不修改本 skill 任何既有内容。
分级含义：A 完整证明 / B 附加条件下成立 / C 仅数值实例检查 / R 文献确认 / G 证据不足。

## 复核分级：C（案例级）→ 引用忠实性已核（R14）

- **G02-CASE-01**（instance_check_only，numeric_experiment）：outer Lp size 实例
  （supavg ≤ L2·√(N/10) 的粗粒度检验通过；注意这是**实例级**，不是定理验证）。
- **G02-REF-01**（instance_check_only→引用核对，source_comparison，R14）：
  SKILL 步骤 4/5/7/9 对 Theorem 4.1（强型 1<p≤∞）、5.1（strong p>2、p=2 仅 weak、
  参数 0<|α|≤1、|β|≤0.9、b≤2^{-8}）、Theorem 1.3（σ-finite、p>a sup / p≤a sup-inf、
  μ(A\B)≤μ(A)/2、C(a,p,q,r) 非统一端点）、Theorem 2.1（p 有限且 p>a、p∞ 需
  Remark 4.3 size–mass）的引用与源卡（5/5 全文阅读）**逐条一致且更细**。

## 更正条目（R14）：档位提升依据

待证据 → **部分可用**（引用忠实性已核 + 案例实例检查通过）。
仍待证据：源定理证明本身（阅读≠重证）；中心立方体接口为推断（skill 已标注）。


---

## 更正条目（R37）：卡↔PDF 保真性验证扩展——G02/D01 主依赖卡通过

- **G02 Theorem 5.1**（PDF p30）：0<|α|≤1、|β|≤0.9、b≤2⁻⁸、强型 p>2、p=2 仅
  L^{2,∞}、常数依赖 α,β,b,φ,p——与源卡/skill 逐字一致 ✓
- **D01 Theorems 1/2/4/5**（PDF p5–6）：对偶空间识别与 KR 对偶、Theorem 4 的
  AC 边缘前提、Theorem 5 的 m>1 反例——与 skill 步骤 3 及 R19 对照表逐条一致 ✓
- 卡↔PDF 累计 5/5 卡，0 偏差；本 skill 的引用链三层闭合（skill↔卡 R19↔PDF R37）。
