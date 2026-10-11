# 第83轮：可调平滑似然与共同容量分块

2026-10-11。**目标未完成：一般 sqrt(n)n^{o(1)} 尚未证明；保留全一般基准为 O(n log n)。**

本轮在任意输入的同一真实 hard 选择器上，证明可数个可调 E2 宽度的极低似然分支并集可按 O(sqrt n) 支付。宽度 eta<=1/n 时，可支付截止值达到 n^{-o(1)}；费用不含每个窄平滑极大算子。固定物理尺度下，全噪声宽度变化也有 O(sqrt n) 预算，但自适应物理尺度 L(x) 的对应积分仍未支付。

共同 Hall 主分块具有严格下降容量、饱和公共流、所有凸负载费用最优性及精确二次回流恒等式。低容量块可支付；高容量块尚无一般根号阶几何预算。真实输入反例否定“小容量差自动带来小捕获损失”，以及对所有最大捕获 Hall 割只支付并集外边界源质量的路线；它们不否定主目标。

本轮使用更新64技能的适用边界。覆盖/软层/几何由 Astra medium 执行；数值与独审由 Sol high 执行。主代理补核硬端点一致收敛、复核论证并在临时目录实际重放。未编译 LaTeX、未生成 PDF。

## 精简归档

[research.zip](research.zip) 内恰有5个必要成员：

- proof.md：合并推导、范围、反例、数值笔记及4个可抽取核验代码块。
- proof.tex.txt：详细 LaTeX 源文本；仅作静态结构检查。
- principal_probe.py：独立标准库完整空间程序。
- principal_results.json：实际来源、完整分区、逐行来源损失、分块与有理流证书，供独立算法审核。
- inventory.json：哈希、技能版本、主代理真实重放输出与证据范围。

解压后执行 `python3 principal_probe.py`。完整便携审计可用以下方式运行：

```python
from pathlib import Path
import re, subprocess, sys
blocks = re.findall(r'```python\n(.*?)\n```', Path('proof.md').read_text(), re.S)
assert len(blocks) == 4
Path('audit_guard.py').write_text(blocks[-1] + '\n')
subprocess.run([sys.executable, '-B', 'audit_guard.py'], check=True)
```

数值程序三阶段覆盖20个实例、123724个完整原空间单元、34个主分块流（含2个零块）、363行、999个空间模式；有理压力搜索360次、种子830612，基础输入种子820611。维数1--3，至多8个来源标签、9个窗口尺度。三个组件检查还包括60个有理图、12组真实L1参数、连续E2几何门槛的全区间有理包络，以及维数16/64/256对应的E2标量诊断。有限网格原子证书、真实连续L1构造和浮点积分分别标记，不互相替代。

本轮只新增此README与research.zip，更新项目总览；本地组件草稿、重复数据和打包工具不提交。历史依赖继续保留，无未经核验的删除。目录宽度与文件规模按 [GitHub官方限制](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits) 核对；总文件数不等同于网页diff显示上限。

归档：250591字节；SHA-256 `d0cf4a31803c17b9febb7504ee0f28d12f8a5bc279c67a8d9a37437d538bf0f2`；恰有5个成员。

规模快照：仓库共2820个跟踪文件，项目1307个；最宽目录821项，最大文件16150642字节。本轮仅新增2个跟踪文件，未触及100MB单文件和2GB单次推送限制；网页diff显示限制不作为仓库文件数量上限。
