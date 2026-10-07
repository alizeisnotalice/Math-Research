

---
## 核心断言验证（R24）

- **步骤 5 softmax 界 = A 级**：max q ≤ β^{-1}logΣe^{βq_j} ≤ max q + log(m)/β
  解析一行（e^{βmax}≤Σ≤m·e^{βmax}）+ 2000 组数值偏差 ≤2.2e−16 ✓
- **步骤 3 指标核**：F(k)(ξ)=∏sin(ξ_j)/ξ_j 可取负（sinc(4)=−0.189）⟹ 立方体指标核
  非 Fourier 正 —— skill 正确区分「正定」与「逐点非负」✓
- 步骤 4：P-8205499d（S¹×S¹ 频谱支集刻画）、P-2ebc469（stretched Gaussian f̂>0）
  均不定义任意立方体核/softmax —— 守卫与源卡一致 ✓

## 档位提升依据（R24）
待证据 → **部分可用**（softmax 界完整证明 + 指标核符号验证 + 引用一致）。
仍待证据：无限指标尾控制；课题 soft-square-supremum 接口。