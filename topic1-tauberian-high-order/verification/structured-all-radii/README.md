# 第18轮：结构化类全部实际半径的最小复算包

本目录只保存最小计算复现材料，不裁决冻结稿之外的数学结论。三份证明稿与三份审查记录按源字节复制；验证目录收录三组 WolframScript 与两组 Fraction 脚本及现有日志。

## 文件映射

| 归档目标 | 冻结源 | SHA-256 |
|---|---|---|
| `proofs/018-grid-all-radii.md` | `r018-beyond-safe/actual_radius_grid.md` | `7224c44ed69231b46a7b3b8f6bc71f6cc2f332ecca6fbeba4f31153184f23dd9` |
| `proofs/018-staggered-rows.md` | `r018-beyond-safe/staggered_rows_2d_kfix.md` | `e2c7938fa347c6957631ebd4b804aa99df0cc7837d7c06797e5288ebce4ee872` |
| `proofs/018-macro-shear.md` | `r018-beyond-safe/macroscopic_shear_nd.md` | `619cbeb62a47d95c9849277f32e71e1fa11b832a4bb84880a3a1d6cffc89593d` |
| `reviews/018-audit.md` | `r018-review/audit.md` | `1307c75c7caebbffc46cb6d122375f943a4bb0c157b4492d00cf69da15d49e3c` |
| `reviews/018-staggered.md` | `r018-review/staggered.md` | `145611f997e6ed4f413df6931a2896a831df7298acc7065306650b799ce9b5ad` |
| `reviews/018-macroscopic.md` | `r018-review/macroscopic.md` | `304ab3d1d56d171c364e94cce2a3bc34d6d9f805ad0d57b4108c0308a26391a3` |

`verification/structured-all-radii/` holds five source/output pairs:

- `grid_checks.wl` + `grid_checks.log`
- `fraction_grid.py` + `fraction_grid.log`
- `row_checks.wl` + `row_checks.log`
- `fraction_rows.py` + `fraction_rows.log`
- `shear_nd.wl` + `shear_nd.log`

## 内容范围

- **近整格、全部实际半径（各维）**：全局整数标号的非周期网格；逐棱法向步长 `1+O_K(u^6)`，切向步长 `O_K(u^3)`，逐格异孔/质量根在统一界内。审查记录给出该结构类的 `binom(n,2)n^2/27` 系数上界、有限物理轴平行盒截断误差及固定且足够大 `K` 时的取等例。不是一般输入的结构约化，也不覆盖任意退化凸体。
- **二维任意行相位错列**：`n=2`，行相位可变，竖棱横移为相位差加 `O_K(u^3)`，其他棱满足相应近单位条件。固定 `K` 下按小相位格面比例给 `4/27` 上界；宏观错列的系数为零；允许相位且 `K` 足够大时可达 `4/27`。
- **高维固定宏观单剪切**：`n>=3`，指定单一仿射剪切参数 `s in [s0,1-s0]`，逐格异孔。比值至少为周期基线的族被迫满足反向剪切 `O(u^6)`，并有报告所列的严格系数损失；固定单剪切的精确最优值仍未知。

三项结论都限于各稿的明示几何类、共同的逐格外包和质量根界，以及物理轴平行有限截断条件 `R u^6 -> infinity`。它们共同补充 R15 下界的特定构造/类上界信息，但不能宣称一般 Lebesgue 输入的精确系数已经确定，不能把模型内最优等同于 `C_n` 全输入最优。

## 审查记录的版本对应

`018-staggered.md` 的主体审查对象是修订前 `staggered_rows_2d.md`；其末尾明确记录了唯一的第 24 行取等量词修订，并核验最终 `staggered_rows_2d_kfix.md` SHA。归档证明文件对应这个最终 kfix 版，未复制旧版。其他两份审查记录直接对应各自归档证明稿哈希。审查中记载的 `proof-valid` 结论仅针对其精确范围。

## 实际运行与清理检查

现有 `grid_checks.log`、`row_checks.log`、`shear_nd.log` 均为 Wolfram Language 15.0.0 日志；各项输出为 `True` 并以 `OK` 收尾。两份 Fraction 日志来自脚本所列有限精确参数样本：`fraction_grid.log` 核对 12 格 Jacobian 积分和跨格有符号平均；每格混合乘积差可以显著非零，报告的平均小负值是其有符号抵消结果。`fraction_rows.log` 的 `captured=False` 是脚本显式断言的预期未捕获样本，不是失败；18 个样本均按预设的真假状态完成。

本次对待复制文件检查未发现私有绝对路径、空文件或失败日志。日志中的 `FAILED` 仅作为 Wolfram 检查脚本的失败分支字符串，实际输出没有该标记。没有发现需要主代理修改或复算的日志项；所有数值证据仍是有限代数/有理样本，不能代替其报告对应的几何证明。

## 复现命令

在本目录运行 `wolframscript -file grid_checks.wl`、`wolframscript -file row_checks.wl`、`wolframscript -file shear_nd.wl`、`python3 fraction_grid.py`、`python3 fraction_rows.py`。全部输出打印到标准输出，按需重定向至对应.log；无额外Python库、私有工作路径或自动时限。三份证明稿头保留原冻结“候选/待审”状态以保持已签SHA，最终独审及二维量词修订状态见项目reviews/018-*.md。本轮主代理未重复这些未改代码的批次。
