# N04 当前审核记录（2026-10-07）

本轮完整证明仅覆盖给定数列递推的条件结论，不代表某个具体熵、能量或立方体极大算子应用已满足前提。

设对所有 `k≥0`，`H_{k+1}≤H_k−d_k+ε_k`，`d_k≥0`，且 `H_k≥H_inf`，其中 `H_inf` 是统一有限下界。对任意 `m≥0` 求和：

`Σ_{k=0}^m d_k≤H_0−H_{m+1}+Σ_{k=0}^m ε_k≤H_0−H_inf+S_m`，其中 `S_m=Σ_{k=0}^m ε_k`。

若 `sup_m S_m<∞`，则 `d_k` 的部分和单调有界，故 `Σd_k<∞`。特别地，逐项非负且可和的误差满足此条件；绝对可和或条件收敛误差的部分和也有界。逐项非负误差不是必要条件。

如果不限制 `S_m` 的上界，取 `d_k=1/(k+1)`、`ε_k=2/(k+1)`、`H_0=0`、`H_{k+1}=H_k+d_k`，递推取等且 `H_k≥0`，但 `Σd_k` 由调和级数比较发散。因此“每步有下降”本身不能给出总下降有限。发散由解析比较证明，有限截断计算不作此证明。

本轮全文读取了两篇仅有方法词汇邻接的论文：Leonetti (P-909d1441688f2461) 的来源PDF位于包内 `evidence/papers/P-909d1441688f2461/paper.pdf`，14页全读；Hou–Wang (P-cf46d6211b76b15e) 的来源PDF位于包内 `evidence/papers/P-cf46d6211b76b15e/paper.pdf`，18页全读并覆盖 pp.15–18 附录。若需打开任一PDF，在保留 bundle root 的环境中对相应 `P-…` 使用统一包级 resolver：`python3 "$EVIDENCE_ROOT/audit_current/delivery/resolve_evidence.py" --skill-dir "$SKILL_DIR" --paper-id P-…`；resolver 会按来源索引核验SHA。前者研究矩阵/算子族的统一可和性及理想收敛，后者研究有理函数的代数移位差分可和性。二者均无上述 `H_k,d_k,ε_k` 下降预算条件，因此不作为N04证明来源。页码和版本SHA分别见包内 `audit_current/ln/paper-reviews/P-909d1441688f2461.json` 与 `audit_current/ln/paper-reviews/P-cf46d6211b76b15e.json`，或索引 `audit_current/ln/papers-read.json`；共享账本可用统一 resolver 的 `--artifact-path audit_current/ln/paper-reviews/P-909d1441688f2461.json`、`--artifact-path audit_current/ln/paper-reviews/P-cf46d6211b76b15e.json` 或 `--artifact-path audit_current/ln/papers-read.json` 定位。单独安装 Skill 不含 resolver、共享账本或PDF语料；需保留 bundle root 并按 [`portable-audit-access.md`](portable-audit-access.md) 通过 `EVIDENCE_ROOT` 或本机 root registry 解析。
