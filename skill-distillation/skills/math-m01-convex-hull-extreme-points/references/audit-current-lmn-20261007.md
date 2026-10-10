# M01 当前审核记录（2026-10-07）

本记录只更新 M01 的有限维支撑界、锐性边界和相关论文范围。它不把本组整体或 Krein–Milman 部分标作完整验收。

## 断言与证明

- `M01-BND-01a`：若 `x∈conv(S)⊂R^d`，则 `x` 有至多 `d+1` 个 `S` 中点的凸组合表示。完整仿射相关消元见 [`method.md`](method.md)“有限维Carathéodory消元”及“上界与锐性的独立证明”；假设仅有限维 `R^d` 与有限凸组合，不需 `S` 紧或闭。
- `M01-BND-01b`：对 `d≥1` 的标准单纯形，重心具有全部 `d+1` 个顶点的唯一重心坐标表示，故统一上界不能减到 `d`。证明见 [`method.md`](method.md)“上界与锐性的独立证明”。有限三角形案例另见 [`examples/positive.md`](../examples/positive.md)。随机抽样不是锐性证据。
- 弱星闭凸包与有限原子凸组合不同：黎曼和给出 Lebesgue 测度为有限 Dirac 凸组合的弱星极限；Lebesgue 测度本身没有点原子，故不是有限 Dirac 凸组合。证明见 [`method.md`](method.md)“弱星闭包与有限原子表示”，负例见 [`examples/negative.md`](../examples/negative.md)。

## 论文范围

本轮全文阅读了包内来源 `evidence/papers/P-3d2c7d539b055402/paper.pdf`，PDF共17页，阅读页1–17，未发现附录，参考文献页也已覆盖。若需打开PDF，在保留 bundle root 的环境中设置 `EVIDENCE_ROOT` 和本 Skill 的 `SKILL_DIR`，然后运行 `python3 "$EVIDENCE_ROOT/audit_current/delivery/resolve_evidence.py" --skill-dir "$SKILL_DIR" --paper-id P-3d2c7d539b055402`；resolver 会按来源索引核验SHA。其 Theorem 2（PDF p.3）是在额外可度量、Banach函数空间、连续限制及分离条件下的 exposed-point 闭凸包结论，证明在 pp.9–10 并依赖外引。该文不是有限维 `d+1` 上界或锐性的来源，也不替代一般 Krein–Milman 定理。

相应页码/条件、版本SHA和全文覆盖记录保存在包内 `audit_current/ln/paper-reviews/P-3d2c7d539b055402.json`；本组相关候选论文的完整未读清单保存在 `audit_current/ln/papers-read.json`。单独安装 Skill 不包含这些共享账本；需保留 bundle root，并按 [`portable-audit-access.md`](portable-audit-access.md) 使用 `--artifact-path` 解析。该总账的 `full_read` 仅登记本轮逐页覆盖的论文。
