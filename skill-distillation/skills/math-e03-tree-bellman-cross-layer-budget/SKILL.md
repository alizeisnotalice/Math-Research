---
name: math-e03-tree-bellman-cross-layer-budget
description: "用于树 Bellman 跨层预算势的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前审核状态（2026-10-07）**：主源 P-ee5e7f1d2f3397ce（13 页、无附录）已全文核读，视觉核回 p2 势定义；准确公式是 `B=4(F−f²/(v+A))`。其余候选源和跨专题迁移仍按共享阅读队列筛查。逐页覆盖与断言见 `audit_current/e/source-reading/`、`claims.json`。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


> **Portable evidence access.** A complete unpacked bundle has sibling `workspace_revised/`, `audit_current/`, and `evidence/` directories; the shared resolver is `<package-root>/audit_current/delivery/resolve_evidence.py`. From an installed Skill, obtain the absolute package root from the unique `math64-20261007` entry in `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` (configure the example registry with the actual absolute root first; never leave `${MATH64_SOURCE_PACKAGE_ROOT}` unresolved). These commands derive paths from that registry, not the current directory; each resolver call clears any inherited `EVIDENCE_ROOT` so it cannot override the selected registry root. `set -e` makes a missing claims file, cases file, or evidence PDF stop the block immediately. If the registry/root/bundle/index is missing or ambiguous, report the configuration gap instead of guessing. Resolver `PASS` checks paths and current bytes, not reading or proof correctness.
>
> ```sh
> set -e
> REGISTRY="${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json"
> ROOT_ID="math64-20261007"
> PACKAGE_ROOT="$(python3 -c 'import json,os,re,sys; from pathlib import Path; d=json.loads(Path(sys.argv[1]).expanduser().read_text(encoding="utf-8")); r=[x for x in d.get("roots",[]) if isinstance(x,dict) and x.get("id")=="math64-20261007"]; assert d.get("schema")=="math-skill-evidence-root-registry-v1" and len(r)==1; v=os.path.expandvars(os.path.expanduser(str(r[0].get("package_root","")))); assert v and Path(v).is_absolute() and not re.search(r"\$[A-Za-z_{]",v); print(v)' "$REGISTRY")" || exit 2
> RESOLVER="$PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
> SKILL_DIR="${CODEX_HOME:-$HOME/.codex}/skills/math-e03-tree-bellman-cross-layer-budget"
> test -d "$PACKAGE_ROOT/workspace_revised" && test -d "$PACKAGE_ROOT/audit_current" && test -d "$PACKAGE_ROOT/evidence" && test -f "$RESOLVER" && test -f "$SKILL_DIR/references/handoff-evidence-index.csv" || { echo "configuration gap: complete bundle, resolver, or installed Skill index missing" >&2; exit 2; }
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/claims.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/cases/cases.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --skill-dir "$SKILL_DIR" --evidence-id E1485 --verify-sha256
> ```


# E03 · 树 Bellman 跨层预算势

当前证据状态：**已核对指定树 Bellman 源定理与常数；迁移时仍须验证树测度、增量和终端条件。**

## 输入与产出

根树 T、层级、节点费用与质量/转移概率、允许的子节点，以及候选预算势 Φ。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

入场时校验树有限/无限深约定、节点质量是否由给定转移概率诱导、预算为逐节点平均还是逐路径、终端收益符号与极限条件。若节点质量、预算类型或终端项缺失，暂停常数引用和望远镜结论；只给出已由输入支持的有限层代数恒等式，并列出失败项。

## 执行步骤

1. 先定义根树、到达质量m(v)、子概率p(w|v)与终端约定，说明逐路径还是平均预算。
2. 平均预算逐节点检Φ(v)≥c(v)+Σ_wp(w|v)Φ(w)；逐路径预算则用sup_w，二者不互换。
3. 乘以m(v)并有限层求和，逐条用m(w)=m(v)p(w|v)抵消内部边；保留终端势项及其符号。
4. P-ee5e7f1d2f3397ce 的 PDF p2 定义 `B(F,f,A,v)=4(F−f²/(v+A))`，4 乘整个括号。适用域至少包括 `f²≤Fv`、`0≤A≤v`；其定义及增量结构不可删改。若把式子误读成 `4F−f²/(v+A)`，主不等式望远镜求和的归一化项不再按原文抵消。每次实际求值另须 `v+A>0`；在 `0≤A≤v` 下以 `v>0` 作明确守卫。若 `v=A=0`，则 `f=0` 而原式是未定义的 `0/0`；不得赋值为 0 或宣称连续延拓，须查明来源是否给出边界扩展或另设零质量分支，否则报告缺口并停止该点认证。
5. 在该域上 `0≤B≤4F`；原文 Theorem 1.3 的常数为 4，Theorem 2.1 的最大函数常数为 32（分别按原文归一化与停时设置）。不可互换两个常数，也不外推到不同树测度或 bi-tree 问题。
6. 复用原文树测度时，重写 λ、Λ(I)、α 及子树 Carleson 条件，核对零质量分母；bi-tree 仅有盒条件不自动满足单树证明条件。
7. 令深度趋于无穷前证明终端项收敛或受控。共享输入、非嵌套中心立方体与终端误差不能由普通 Bellman 望远镜自动编码。

## 证据与失败处理

平均预算与逐路径预算不同。节点质量、终端余项或无限深极限任一未经核验，都会破坏望远镜求和。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
