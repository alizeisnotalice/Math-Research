# 第75轮：软层覆盖反例、远尾统一支付与定价证书

2026-10-11。**一般 sqrt(n)n^{o(1)} 上界仍未证明**；保留全一般 O(n log n) 基准。

- 固定代表的价格多面体有精确贪心刻画；优化代表后的列系数既不普遍超模也不普遍次模。一般分离来源块的 T、B_t、O_t 精确相加，全部输入可无损归约到按自身质量半径连通的块，无递归深度收费。
- 固定 n=2 的真实完整低背景分支严格否定未截断物理门的任何一致正覆盖率。这是中间覆盖命题的反例，不是弱型反例。
- 对任意源、任意有限尺度族 subset [1,2]，半径 R_n(a)=1+2 log[2n(1+n log2)]+2a 外，所有软层远尾一起强 L1 付费至多 2 exp(-a)||f||_1。近区覆盖及全尺度回代仍待完成。
- 费用归一化的逐源核支配门修复该反例；bbox 包络与混合门支付式成立。全部截断核心可网格覆盖，但费用仍为 O(n log n)。真实近区子事件排除廉价 bbox-only 普遍覆盖，不否定完整分支或混合门路线。
- q=2 固定代表最小割分离全部来源子集；目标指数 q=1+1/k 的有理凸包/辅助节点定价保留 epsilon*T 误差预算，不额外乘原子数。

数值：10个较大有限训练模型、2,560训练点及10,240独立验证点，全部有理上下界通过审计；其中9个量化模型原/对偶等目标，64点平面稠密源保留约0.0363273间隙。目标指数测试覆盖 n=1,4,9,16,64，含150次全子集定价对照、30个有理LP证书及并列赢家/粗舍入回归。所有有限模型与浮点诊断均明确限定范围，不推断连续覆盖或渐近阶。

父代理独立重放 q2 主脚本4.587秒、目标指数脚本0.682秒，两个JSON逐字节一致；独立综合审计7.035秒。几何、实际软层、核包络及独立软层守卫全部通过。精确stdout、源码哈希、数值误差和复现范围见总报告及inventory。

[research.zip](research.zip) 恰含7项：proof.md、proof.tex.txt、pricing_probe.py、pricing_results.json、parent_pl.py、parent_pl_results.json、inventory.json。LaTeX以txt保存，未编译、未生成PDF。仅依赖Python标准库，实际运行Python3.12.14。解压到临时目录并保留JSON备份后运行两份主脚本。五份守卫内嵌于proof.md，提取如下：

```sh
python3 - <<'PYGUARD'
from pathlib import Path
import re
blocks = re.findall(r"```python\n(.*?)\n```", Path("proof.md").read_text(), re.S)
names = ["coverage_guard.py", "soft_guard.py", "envelope_guard.py", "audit_guard.py", "soft_audit_guard.py"]
assert len(blocks) == len(names)
for name, code in zip(names, blocks):
    Path(name).write_text(code + "\n")
PYGUARD
python3 -B pricing_probe.py
python3 -B parent_pl.py
python3 -B coverage_guard.py
python3 -B soft_guard.py
python3 -B envelope_guard.py
python3 -B audit_guard.py
python3 -B soft_audit_guard.py
```

覆盖、软层及核包络使用Astra medium；数值与独审使用6.1 Sol high。64项canonical技能入口及437处直接引用已重新核验；其可用性不等于全部文献或定理已验收。仅新增本说明和压缩包，并更新项目总览；草稿留在仓库外，历史依赖保留，未宣称全量历史去重已完成。

归档：283,340字节，SHA-256 `8482f5e6edaaf86600eb874285d476139896d7ca5e07335e5e2340fc78cb89fa`。

已复核[GitHub官方限制](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)：100 MB单对象、2 GB推送为强制限制，3,000条目录宽度为建议。网页diff的300文件显示上限不是仓库总文件数上限。本轮仅新增两个Git文件，归档亦低于1 MB单对象建议值。

提交前工作区规模：含本轮共 2804 个文件，项目 1291 个；最宽目录 821 项，最大文件 16,150,642 字节。Git对象包261.91 MiB、松散对象12.76 MiB；这是本轮整合远端前的快照，未把其他项目后续提交计作本轮新增。
