---
name: math-e01-occupation-measure-flow-lp
description: "用于占用测度 LP 与流守恒的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


> **Portable evidence access.** A complete unpacked bundle has sibling `workspace_revised/`, `audit_current/`, and `evidence/` directories; the shared resolver is `<package-root>/audit_current/delivery/resolve_evidence.py`. From an installed Skill, obtain the absolute package root from the unique `math64-20261007` entry in `${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json` (configure the example registry with the actual absolute root first; never leave `${MATH64_SOURCE_PACKAGE_ROOT}` unresolved). These commands derive paths from that registry, not the current directory; each resolver call clears any inherited `EVIDENCE_ROOT` so it cannot override the selected registry root. `set -e` makes a missing claims file, cases file, or evidence PDF stop the block immediately. If the registry/root/bundle/index is missing or ambiguous, report the configuration gap instead of guessing. Resolver `PASS` checks paths and current bytes, not reading or proof correctness.
>
> ```sh
> set -e
> REGISTRY="${CODEX_HOME:-$HOME/.codex}/math-skill-evidence-roots.json"
> ROOT_ID="math64-20261007"
> PACKAGE_ROOT="$(python3 -c 'import json,os,re,sys; from pathlib import Path; d=json.loads(Path(sys.argv[1]).expanduser().read_text(encoding="utf-8")); r=[x for x in d.get("roots",[]) if isinstance(x,dict) and x.get("id")=="math64-20261007"]; assert d.get("schema")=="math-skill-evidence-root-registry-v1" and len(r)==1; v=os.path.expandvars(os.path.expanduser(str(r[0].get("package_root","")))); assert v and Path(v).is_absolute() and not re.search(r"\$[A-Za-z_{]",v); print(v)' "$REGISTRY")" || exit 2
> RESOLVER="$PACKAGE_ROOT/audit_current/delivery/resolve_evidence.py"
> SKILL_DIR="${CODEX_HOME:-$HOME/.codex}/skills/math-e01-occupation-measure-flow-lp"
> test -d "$PACKAGE_ROOT/workspace_revised" && test -d "$PACKAGE_ROOT/audit_current" && test -d "$PACKAGE_ROOT/evidence" && test -f "$RESOLVER" && test -f "$SKILL_DIR/references/handoff-evidence-index.csv" || { echo "configuration gap: complete bundle, resolver, or installed Skill index missing" >&2; exit 2; }
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/claims.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --artifact-path audit_current/e/cases/cases.json
> env -u EVIDENCE_ROOT python3 "$RESOLVER" --root-id "$ROOT_ID" --registry "$REGISTRY" --skill-dir "$SKILL_DIR" --evidence-id E1970 --verify-sha256
> ```


# E01 · 占用测度 LP 与流守恒

**本轮审核状态（2026-10-07）**：核心来源 P-610d8ca4e35b8a23（29 页，含附录）已逐页核对并在原PDF复核关键公式；有限MDP、Dynkin流、Cantelli/VP 分支及 RU/measure ES 可在显式假设下使用。Theorem 5.5 与固定标量终止时刻问题之间存在同域反例；原子分位处印刷式 (3) 与后文 RU/LP 标准 ES 不一致。完整证据卡见上级审核工作区的 `audit_current/e/`。其余候选直接/邻近来源及共享阅读队列仍需按独立SHA逐篇验收；本状态不代表 E01 全部引用已完整验证。

## 输入与产出

状态空间 S、动作集 A、初始分布 μ₀、转移核 P（或生成元）、折扣因子 γ/时域、奖励或成本 r(s,a)。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

入场时先明确目标是轨迹控制值、弱占用测度 LP、峰值尾概率还是 ES；核对 μ₀ 为概率、时间质量约定、测试函数域和生成元可积性。ES 还需给出尾概率 ε 与分位点原子处理。缺初始律、终止规则或测试域时，停止等价与无间隙结论标为未验证；不得用 LP 可行性代替轨迹可实现性。输出同时列明输入检查、所用来源页码、结果性质（定理/推导/数值）及失败原因。

## 执行步骤

1. 先选有限γ折扣MDP、定时域确定性控制、扩散/跳控制、停止风险、均场相互作用或变分域模型；时间/折扣/停止occupation的质量及界方向不同，均场二次目标不默认纯LP。
2. 有限MDPγ∈[0,1)、μ₀概率、P为归一化转移核：x(s,a)=EΣγᵗ1_(sₜ,aₜ)≥0，检Σ_ax(s,a)−γΣ_(s′,a)P(s|s′,a)x(s′,a)=μ₀(s)，总质量1/(1−γ)。正占用处恢复π=x/Σx，零占用状态另定。该有限MDP结论由审核包中的独立有限维推导支持，不由 P-610 的连续时间 Dynkin 方程推出。
3. 确定性OCP对C¹测试q保初末Liouville等式。扩散P-3b44b82c06eccd85需C^{1,2}及Dynkin可积性，Aq=∂_tq+f·∇q+(1/2)tr(ggᵀ∇²q)；终端质量1、运行occupation质量T。跳过程改完整跳生成元，不能漏跳项。
4. 最小成本弱measure LP可能严格低于原轨迹控制最小值；多项式紧域球约束下moment-SOS层级渐近趋该弱LP，不认证任意有限阶或真实控制恢复。
5. 均场来源P-078f65ff156d77ea与P-22e8303b48a27276用同一时刻人口边缘的二次交互。P-078的PSD条件使未缩放交互泛函F凸；原目标含λF时还须λ≥0（λ=0合法），W PSD本身不能推出λ<0时凸。P-22也须核对矩阵核在所有有符号测度上的PSD及交互权符号。冻结梯度后的LMO线性分离不证明原目标凸；exact LMO、可测经典最优选择、状态约束内点和smoothness/曲率条件才调用FW速率。详见[分模型来源与定位](references/frozen-source-locators.md)。
6. 变分源P-b8164dc1ce4dcf5d的控制无gap结论除OC1–OC3外，还须保留Section 5.1基线：Ω有界连通、边界分片C¹；U,U∂完备可分，约束数据Borel可测、成本可测且局部有界，状态属W¹,p、控制可测；p<∞时测度有相应p阶矩，p=∞时支撑紧。一般仅M_aff≤M_r≤M。reachability体积、危险区期望停留时间和最小成本的目标/界方向分开；有限矩匹配不能恢复真实控制轨迹，详见补充审计。
7. 局部分片来源使用有限确定性时空分区；不要把受控过程/跨片停时的适应性误写成分区适配。局部primal聚合须固定初始节点为全局初始律在各片上的限制、相邻时间片共享同一节点测度，并保反对称signed界面通量以验证抵消；扩散空间界面要连续，确定性单向放宽不适用扩散。作者有限阶更紧比较另需闭片的严格更丰富多项式表示，非爆炸/停时及通量可定义性另核。
8. P-610 的 peak measure/SOC 给上界：Cantelli 使用 `r=√(1/ε−1)`；unimodal VP 简化分支另要求 `0<ε≤1/6`。这对应 VP 单侧不等式分支条件 `r²≥5σ²/3`，不声称全局最大有效 ε。Problem 3.1/5.1 的 `t*` 是确定性标量，而 (29) 不施加共同确定终止时刻。P-610 Theorem 5.5 的字面无间隙结论与其固定标量终止时刻问题 (27) 及随机停止测度 LP (29) 的可行域不一致。`audit_current/e/independent-checks/E01-theorem5-5-terminal-time-gap.md` 给出紧区间停止布朗运动见证：固定时刻 ES <0.91994，而随机停止/LP 值为 1；该 SDE 类的 A5–A7 依据原文 p7 对紧域 SDE 的明确分类（外引 [20]，本地未独立取得该文）。同一文件的两态纯跳例固定时刻 ES ≤0.19802、随机停止值 1，仅作辅助数值见证，不单独声称已证明有限时域时间空间图核心 A5–A7。结论只针对固定标量 `t*` 与 LP (29) 的等价，不削弱 Theorem 5.4 的上界方向、LP 本身或改为随机停止后的命题。
9. 标准 ES/RU（(4)）与 measure LP（(5)、(28)、(29)）在 `0<ε≤1`、`ψ=p#μτ` 为概率测度时允许拆分分位点原子：`εν+νhat=ψ` 且 `ν` 质量为 1。P-610 的 Definition 2.1 式 (3) 对整个事件 `{V≥VaRε(V)}` 积分后除以 ε；阈值为原子时，该量不等于 RU/LP ES。它是印刷的另一尾事件统计量，不能与后续标准 ES 等同或静默替换。例 `V=3(.2),1(.3),−2(.5), ε=.4`：式 (3)=2.25，RU/LP拆分阈值原子给2.0。
10. 作者 tail-SOC 强对偶附录 B 的 (60c)–(61) 用 `E[p]≤Π₁` 控 `E[p]²`，对可能为负的 `p` 不成立。概率测度下 Jensen 给 `(E p)²≤E[p²]≤Π₂`，可把 trace 界局部修为 `6(1+2Π₂+Π₂²)`；这不替代对抽象对偶、停止表示与层级收敛外引条件的核验。
11. 闭环控制、稀疏分区实验、数值重构与论文理论分开；无限离散跳态的半代数外包、所有共同输入和中心立方体流接口另证。
12. P-0297391ab7ac5c7b finite-mode切换：μ_j为u_j加权投影；逆向以Σμ_j及δ_ej合并，weak LP等价不等于原switching无gap。用v=t检Σ质量T，time-only tests检Σπ_tμ_j=dt；joint reconstruction须保共同time marginal并真实积分验证。Archimedean polynomial hierarchy仅渐近给weak LP lower bounds；source p24乘性dynamic与linear chattering计算矛盾，Cor2 support无权积分暂不作认证。

## 证据与失败处理

occupation质量与界方向取决于模型；弱LP未必无gap或可恢复轨迹，有限阶数值亦不认证。随机局部分片需合法通量/界面条件；peak来源的确定性共同终止时刻问题与未限制共同终止时刻的停止测度LP不等价、ES分位原子错误及tail-SOC符号界缺口保留，不能当一般无条件occupation定理。 人口交互一般为二次测度泛函，核PSD或权重符号缺失不认证凸性；精确oracle速率不覆盖Adam/网格近似，无gap变分定理不泛化到任意动态控制。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
