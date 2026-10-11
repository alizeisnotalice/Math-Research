---
name: math-e04-snell-envelope-optimal-stopping
description: "用于最优停止与 Snell 包络的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前审核状态（2026-10-07）**：P-faa1738580ee2154（15 页）已全文核读，识别出随机化停止质量不足以代表强制停止的边界、τ*=0 时打印优化器映射的端点问题，以及 singular-control 跳点逆映射的约定缺口；普通有限离散 Snell 递推另作独立验证。其他候选文献仍待相关性筛查。页覆盖和原文定位见 `audit_current/e/source-reading/` 与 `claims.json`。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


> **Portable evidence access.** A complete unpacked bundle has sibling `workspace_revised/`, `audit_current/`, and `evidence/` directories; the shared resolver is `<package-root>/audit_current/delivery/resolve_evidence.py`. From an installed Skill, obtain the absolute package root from the unique `math64-20261007` entry in `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` (configure the example registry with the actual absolute root first; never leave `${MATH64_SOURCE_PACKAGE_ROOT}` unresolved). These commands derive paths from that registry, not the current directory; each resolver call clears any inherited `EVIDENCE_ROOT` so it cannot override the selected registry root. `set -e` makes a missing claims file, cases file, or evidence PDF stop the block immediately. If the registry/root/bundle/index is missing or ambiguous, report the configuration gap instead of guessing. Resolver `PASS` checks paths and current bytes, not reading or proof correctness.
>
> ```sh
> set -e
> REGISTRY="${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json"
> ROOT_ID="math64-20261007"
> PACKAGE_ROOT="$(python3 -c 'import json,os,re,sys; from pathlib import Path; d=json.loads(Path(sys.argv[1]).expanduser().read_text(encoding="utf-8")); r=[x for x in d.get("roots",[]) if isinstance(x,dict) and x.get("id")=="math64-20261007"]; assert d.get("schema")=="math-skill-evidence-root-registry-v1" and len(r)==1; v=os.path.expandvars(os.path.expanduser(str(r[0].get("package_root","")))); assert v and Path(v).is_absolute() and not re.search(r"\$[A-Za-z_{]",v); print(v)' "$REGISTRY")" || exit 2
> RESOLVER="$PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
> SKILL_DIR="${CODEX_HOME:-$HOME/.codex}/skills/math-e04-snell-envelope-optimal-stopping"
> test -d "$PACKAGE_ROOT/workspace_revised" && test -d "$PACKAGE_ROOT/audit_current" && test -d "$PACKAGE_ROOT/evidence" && test -f "$RESOLVER" && test -f "$SKILL_DIR/references/handoff-evidence-index.csv" || { echo "configuration gap: complete bundle, resolver, or installed Skill index missing" >&2; exit 2; }
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/claims.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/cases/cases.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --skill-dir "$SKILL_DIR" --evidence-id E0133 --verify-sha256
> ```


# E04 · 最优停止与 Snell 包络

当前证据状态：**有限离散 Snell 包络与特定随机化停止源已核对；连续时间扩展和源定理的全符号收益等价均需额外条件。**

## 输入与产出

滤过概率空间、适应奖励过程 (Xₜ)、停止时域、可积性类别和允许的停止规则。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

入场时记录有限离散/连续时间、奖励过程适应性与可积性、停止时是否允许随机化、未停止时收益、终端时刻是否强制停止，以及连续时间所需的路径正则性。条件不全时，仅在明确的有限离散条件下使用递推；不得据子概率停止结果推断强制停止值，也不得在缺少右连续或 class-D 等条件时宣称连续时间 Snell 包络结论。

## 执行步骤

1. 离散有限时域且X_t适应可积时，从S_T=X_T逆推S_t=max(X_t,E[S_{t+1}|F_t])；以条件期望检超鞅与支配。
2. 对任意支配X的可积超鞅Y逆推得Y_t≥S_t，证明最小性；τ*=inf{t:S_t=X_t}为停止时，有限停止后的S保持鞅给E X_{τ*}=E S_0。
3. 随机化停止另给 H 适应右连续非降 G、初始值、终端总质量和未停止收益；质量≤1 是 subprobability，只有明确零收益/未停止选项及终端约定后，才能与强制 τ≤T 比较。
4. P-faa1738580ee2154 Theorem 3.1 的随机化停止约束只要求 G(T)≤1。有限时域、k≡−1 时，强制停止值为 −1，而 G≡0 给随机化问题值 0；因此全符号收益下不能宣称二者等价。要改成总质量 1 并处理 t=0，或加入明确的不停止零收益选项后另行证明。
5. 对广义逆 `α(r)=inf{s:G(s)≥r}`，核实事件不等式方向与右连续约定；跳跃处一般不能写 `G(α(r))=r`。Stieltjes 换元应以推前恒等式 `α#Leb_(0,G(T)] = dG` 证明，并注明端点质量与可积性。该恒等式不修复负收益下的总质量差，也不自动给出逐分位最优停止对应。另核对精确优化器映射：Theorem 4.2(a) 打印的 `G*(t)=1_{t≥τ*>0}+1_{τ*=0}` 在 `P(τ*=0)>0` 时违反 Problem 2.3 的 `G(0)=0`；值可逼近或由另一随机化策略达到，不等于该映射本身可行。
6. singular-control 跳点逆映射须明确 `dG=e^ξdξ` 和 Stieltjes 积分使用的前/后跳值。按源 Problem 2.4 的字面 post-jump 被积函数 `e^{-ξ(s)}` 解释时，单位 `G` 原子不能由有限 `ξ` 跳满足该关系；原文没有明确这一跳值约定，故将 Theorem 4.2(d) 的 jump optimizer correspondence 保持为 convention-dependent unresolved，不把条件性计算泛化为对所有解释的反例。

7. 连续时间Snell需右连续、class-D等适用假设；一般信息流论文没有证明经典Snell最小超鞅定理，不能充当该定理来源。
## 证据与失败处理

有限离散Snell是已检查构件。一般信息流来源允许随机化剩余质量，有限T负奖励给其值等价反例；缺终端质量/未停收益约定时不引用该等价。连续时间和无界停止需独立正则性及可选抽样条件。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
