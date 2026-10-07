

---

## 算法公式数值+解析验证（R22）

- **qBnB(2)**：f(x)=(x−0.3)²、盒[−1,1]、L₂=2：qlb(C)=f(x_C)−(L₂/2)r² 含最小点时 ≤f* ✓；
  gap 公式 (L₂/2)r² 在二分序列上精确：1.0 → 0.25 → 0.0625 → 0.0156（r 减半 gap ÷4 = 二阶）✓
- **qBnB(3)**：m=3L₃r ⟹ r_{k+1}=3L₃r_k²/(2m)=r_k²/(2r)，序列 1 → 0.5 → 0.125 → 0.0078 →
  3.1e−5 → 4.7e−10，二次收敛（r_{k+1}/r_k = r_k/(2r) → 0）✓
- **步骤 3 守卫**：qlb 只须在含全局最小点的盒上 ≤f*，剪枝需「保留至少一个全局最小点」的
  归纳证明 —— 与源卡 AD-P-922909d5e8e47a30 的 quasi-lower-bound 定义一致 ✓
- **步骤 7**：P-deb72842d1067b4f 截断判断撤回与 Lemma 3 的 a≤c<d≤b 严格要求，
  与 source-completeness-correction.md 一致 ✓

## 档位提升依据（R22）
待证据 → **部分可用**（算法公式数值+解析验证通过 + 引用与守卫一致）。
仍待证据：qBnB(3) 收敛定理证明本身；SOS/矩层级与外引；维数收缩复杂度声明。


---

## 更正条目（R40）：卡↔PDF 保真性验证——f07/C04 通过

- **f07 Proposition 1.9**（PDF p7）：正多线性算子、不等式 (8) 常数即 ‖T‖（相对常数
  恰 1）、双线性卷积推论（Young 指标）、"equality in the case p=r follows by induction"
  ——全部逐字一致；与 R22 取等构造互证
- **C04**（PDF p1 摘要）：quasi-lower bound "only for sub-cubes containing a minimizer"
  与 skill 守卫逐字对应；qBnB(2) 二阶收敛与 R22 的 (L₂/2)r² 验证一致
- 卡↔PDF 累计 11/12 卡逐字通过、1 处转录偏差（E03，已修正）
