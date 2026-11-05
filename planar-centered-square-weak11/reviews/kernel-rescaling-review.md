# portable-v1 最终独立审查

2026-11-05。结论：**PASS，可按当前声明范围交付；无须修订数学核验。** 本次审查没有修改 `scratch/kernel-public-package`，只在本审查目录保存副本和运行证据。

## 审查对象身份

冻结公开包共有22个文件、893402字节。实际重新计算的两个外部锚点与所给值相同：

- `package-manifest.json`：`9455aec41dd26d08d1812430cc2946cb34e3387d2fe371ee1f245f5bd64d9cb3`。
- `release-receipt.json`：`a9d1f89c60302fc897d8241314a473af4a428f1dd6ba4a528cc0ac133f5998d6`。

receipt列出的其余21个文件逐项核对哈希通过，集合等于实际交付文件集合减receipt，无遗漏。静态manifest绑定16项输入，其余是manifest本身与4项历史包验证输出。没有符号链接。vendored `hl_check.py`、`hl_core.py`、`LICENSE`与原公开clone逐字节相同，wrapper记载的原审查脚本SHA也与本目录原脚本相同。公开包的运行输出、receipt不伪装为额外数学前提。

## 原脚本到便携脚本的数学检查保留情况

已全文检查冻结版 `verify.py`，与原审查脚本及本次设计要求对照。

| 项目 | 最终状态 |
|---|---|
| 精确缩放 | 保留：c和全部1181类a以Fraction逐项比较，类键集合相同，全部其他数学字段/键集合相同，新claim精确有理相同。 |
| C1/C2新计算 | 保留：旧/新核分别实际跑`check_c12`，192-bit，无失败，cell计数，C2直接取`lower().fmpq()`与κ/1比较。 |
| C3独立积分 | 保留：不用`check_c3`，展开全部整数索引轨道求和并直接评价积分，输出精确有理端点，比较新质量上下简单有理界及claim；还增核旧质量小于原claim。 |
| 旧C4完整覆盖 | 保留：`args.only`全C1-C4、`args.cells`空；expected集合、报告集合与1368计数同时相等，排除遗漏或重复替代；全部worker无失败，叶数和443754。 |
| 旧C4身份 | 保留：旧证书SHA、旧checker/core源码SHA与包内实际内容绑定，所有worker的spec SHA及历史provenance摘要一致。 |
| 新C4 | 保留：新旧支持cell集合相同；数学证明使用精确正比例的分布线性传递，没有冒称重跑新C4。 |
| 旧13项冻结身份 | 明确改为只核实际复制的5项子集，原13项manifest本身哈希完整保留。遗漏的大进程采样、旧runner、环境和其他输出不参与数学证明。 |
| 主机依赖 | 删除旧provenance必须等于本机provenance的断言；分别记录历史与实际provenance，checker/core源码内容仍必须相同。输入全部相对`__file__`定位。 |

没有把输出float或Arb十进制字符串提升为证明前提。`-O`模式明确拒绝，防止Python关闭数学断言。`requirements.txt`固定python-flint 0.9.0；包没有绑定旧主机的FLINT二进制文件或旧绝对路径。

## 独立实际移位执行

本次另复制整个冻结包到 `scratch/kernel-rescaling-review/portable-v1-independent-copy`，使用已有数学环境，从`cwd=/tmp`以`-I -B`运行副本的`verify.py`，输出写在本审查目录。返回码0，用时7.45秒，状态PASS，精度192位。实际加载的checker/core路径确实位于该副本；新运行provenance与历史provenance不同。包内五项历史子集、16项静态输入、全缩放及新旧C1/C2/C3精确比较均通过。

新质量仍严格满足

\[
3.590867754338446537
<\tfrac12\int K_{\rm new}
<3.590867754338446538
<3.590867754339.
\]

运行前后重新核公开包receipt及其21文件哈希，原包保持冻结。没有重跑C4，也没有安装任何包。

运行证据见 `portable-v1-audit-result.json`、`portable-v1-independent-verification.json`及同名前缀stdout/stderr日志。新验证输出SHA256为 `2cb3f3fff310cbe29f331f2f6f7c183165064e2995cb2a410aea80bbfb98b7fe`。本次核验程序为 `audit_portable.py`。

## 范围说明和数学结论

README明确说明旧C4报告是完整运行的逐knot-cell汇总，不是逐叶Bernstein/余项证书；默认核验检查覆盖和历史provenance，没有数值重算每一叶。PROOF同样明确以旧完整C4运行记录为前提，正线性传递新C4。README给出原核及新核的可选完整重跑命令，命令引用的输入、checker和依赖都在包内或依赖列表中，支持换主机重新执行。

PROOF正确保留Theorem2.2/Lemma3.1依赖、原点外分布条件、二维原点可积性、整体支持外尾部及r>R样条部分；没有误加AK原点附近可积或有限测度假设。全部正尺度和一般L1量词、菱形到中心轴平行正方形的线性变换、严格/非严格水平集的全局常数口径均清楚。

因此包的声明与实际工作相符：在原审报所核的Lin–Tomar比较原理和旧全C4运行基础上，给出一般L1、全部正尺度的新可认证严格上界 \(C_2<3.590867754339\)。没有宣称新的Lean构建、核族最优、锐常数或极值达到性。包的便携化没有削弱这项数学结论。
