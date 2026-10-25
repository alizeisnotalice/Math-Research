---
name: math-e02-nonlocal-hjb-supersolution
description: "用于非局部 HJB 超解的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前审核状态（2026-10-07）**：主源 P-9873072d1f9301d3（27 页）及其引用的 Sun–Xie–Xie arXiv:2004.02595v2（24 页，含附录）均已逐页核读；后者的 Theorem 2.3 是无控制自治 SDE 结果，向受控/路径依赖终端项 (4.48) 的迁移仍未验证。Theorem 2.2 的速率陈述与 §4.2 末式 (4.49) 存在证明支撑缺口，不能据此称定理为假。页覆盖、定位与分级见 `audit_current/e/source-reading/` 和本组 `claims.json`；其他候选源仍待相关性筛查。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


> **Portable evidence access.** A complete unpacked bundle has sibling `workspace_revised/`, `audit_current/`, and `evidence/` directories; the shared resolver is `<package-root>/audit_current/delivery/resolve_evidence.py`. From an installed Skill, obtain the absolute package root from the unique `math64-20261007` entry in `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` (configure the example registry with the actual absolute root first; never leave `${MATH64_SOURCE_PACKAGE_ROOT}` unresolved). These commands derive paths from that registry, not the current directory; each resolver call clears any inherited `EVIDENCE_ROOT` so it cannot override the selected registry root. `set -e` makes a missing claims file, cases file, or evidence PDF stop the block immediately. If the registry/root/bundle/index is missing or ambiguous, report the configuration gap instead of guessing. Resolver `PASS` checks paths and current bytes, not reading or proof correctness.
>
> ```sh
> set -e
> REGISTRY="${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json"
> ROOT_ID="math64-20261007"
> PACKAGE_ROOT="$(python3 -c 'import json,os,re,sys; from pathlib import Path; d=json.loads(Path(sys.argv[1]).expanduser().read_text(encoding="utf-8")); r=[x for x in d.get("roots",[]) if isinstance(x,dict) and x.get("id")=="math64-20261007"]; assert d.get("schema")=="math-skill-evidence-root-registry-v1" and len(r)==1; v=os.path.expandvars(os.path.expanduser(str(r[0].get("package_root","")))); assert v and Path(v).is_absolute() and not re.search(r"\$[A-Za-z_{]",v); print(v)' "$REGISTRY")" || exit 2
> RESOLVER="$PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
> SKILL_DIR="${CODEX_HOME:-$HOME/.codex}/skills/math-e02-nonlocal-hjb-supersolution"
> test -d "$PACKAGE_ROOT/workspace_revised" && test -d "$PACKAGE_ROOT/audit_current" && test -d "$PACKAGE_ROOT/evidence" && test -f "$RESOLVER" && test -f "$SKILL_DIR/references/handoff-evidence-index.csv" || { echo "configuration gap: complete bundle, resolver, or installed Skill index missing" >&2; exit 2; }
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/claims.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/cases/cases.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --skill-dir "$SKILL_DIR" --evidence-id E0567 --verify-sha256
> ```


# E02 · 非局部 HJB 超解

当前证据状态：**已核对特定论文与方法边界；不构成对所有受控非局部 HJB 模型的通用定理背书。**

## 输入与产出

受控状态过程、漂移/扩散、跳跃核与补偿约定、运行奖励、终端及边界数据、候选值函数。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

入场检查应记录时域、状态/控制域、控制的可测或变差条件、生成元定义域、Lévy 截断和补偿、奖励符号、终端/边界条件及候选函数正则性。缺任何会影响生成元或 Itô 项的约定时，不给出超解认证；把可直接验证的残差与依赖文献迁移的部分分开报告，并定位每个失败条件。

## 执行步骤

1. 从给定受控过程写生成元，列漂移、扩散、Lévy核、截断及补偿约定；先确定最大化奖励或最小化成本及时间方向。
2. 局部接触的非局部测试将小跳作用于测试函数、大跳作用于候选解；不得只对局部PDE检残差，或任意删补偿项。
3. 光滑候选逐项验证残差、终端及边界支配，再用适用的含跳Itô与停止论证验证界；非光滑候选需完整viscosity定义和比较条件。
4. P-9873072d1f9301d3 的模型限于原文假设：对称加性 α-stable 噪声、α₁,α₂∈(1,2)、快漂移强耗散及规定的正则/增长条件。有效 Hamiltonian (2.6) 是对原控制集 U 先取上确界再积分；一般不能与先对快状态取平均再取 Bellman 上确界交换。原文 (2.8) 使用扩张控制集 U^ex 才有 Bellman 表示，迁移时须保留该区别。
5. Theorem 2.2（p5）**陈述** `|uε−ū|≤C(P)ε^p`，其中 `1<p<2α₁α₂/(α₁+2α₂)`。该指数区间非空当且仅当 `α₂>α₁/[2(α₁−1)]`；结合 `α₂<2` 必有 `α₁>4/3`。这是陈述域的代数条件，不等于证明已支持该速率。
6. §4.2 证明末端 p20–23 的 (4.49) 写为 `|uε−ū|≤C(P)(ε^{p(1−1/α₂)}∨ε)`；它与 Theorem 2.2 的指数表达不一致，记为 `source_proof_gap`。目前只报告“陈述与展示的证明末式不匹配”；不得称定理已被反驳，也不得把中间估计导出的候选弱指数升级为已验证定理。
7. 对 (4.36)，一般控制路径若进入 corrector 链式法则，`∂vΦ·dv/dr` 需要绝对连续/有限变差及相应链式法则条件；在原文精确的 Ab,L 快慢解耦假设下，Poisson corrector Φ 不依赖 v，此项才消失。单独这一处不足以推出任意控制版本失败。终端估计 (4.48) 引用无控制的 Sun–Xie–Xie Theorem 2.3；向受控、可能路径依赖的模型迁移所需条件与推导未建立，记为独立的 `source_transfer_gap`。
8. 原文半极限→Liouville→corrector→perturbed test→terminal layer→comparison 是作者给出的证明结构；本卡覆盖论文页码，但不把每个引理都称为已独立复现。常数边界强最大值原理的引用范围限于原文条件，不外推到全空间比较定理。
9. 不外推至 α<1、乘性噪声或任意跳核；模型正则性、增长、可积性或极限顺序缺失时只报告条件结论或候选超解，并逐项列出缺项。

## 证据与失败处理

超解方向依赖奖励/成本和时间约定；小跳/大跳及补偿不能省略。原文特定多尺度stable模型的局部一致收敛无通用统一速率或中心立方体结论。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
