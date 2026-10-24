# 第36轮最小积分契约补丁审查

**最终静态 gate：GO。** 仅允许固定B1核的一次MKA补验，Q8、Adepth7、Kdepth6、6工人、checker运行180秒进程组硬限。尚未执行数学；本审计没有导入flint、调用积分、评价核、生成模型或运行LP/checker。GO不预言补验通过，不产生新常数。

## 冻结来源与独立静态核对

已全文审`make_patch.py`、精确unified diff、两处修正函数、PLAN及`run_once.py`。另实际用标准库读取文件、核hash、重建diff、解析AST；未执行补丁生成器或数学模块。`approved_inputs.json`内5项全部与实际文件一致：

| 文件 | SHA256 |
|---|---|
| fixed_checker/hlcert.py | `4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05` |
| fixed_checker/cert_run.py | `079906db354febdd57d537f0d747428bc90bd37fa88297449cd1d46c7446310d` |
| make_patch.py | `ca5302ef0ad515fab7503306882c8dee9db3c8440b8efee1c71963ea4a070974` |
| patch_manifest.json | `3ac50400cd7603476c18369acf2daf8e7678cd6757fdf2862188d8123a82b622` |
| run_once.py | `9a618bbf9920e0c786721f22da1ee78087fd855e1cf57c014b850a7fa00a42e8` |

原`hlcert.py`仍为`4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5`；修正driver与原driver逐字节一致。实际diff准确等于保存的`hlcert.patch`，其SHA为`b3ac6a012c3ac942e87c74a6a7d9187d256b7f00db62eaf3a31c844e8f11ee01`。原核B1固定SHA为`21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`。

独立AST核对准确限定新增`analytic`关键字、末位`an`实参及局部pw的`an`形参；擦去这些转发后，完整AST与原件一致。必须按末参名称限定，不能误删原`S_power`内另一个同名pw的正常三个参数。

## 补丁覆盖与契约

`Apsi_outside_point`的三个复log均变为`log(analytic=analytic)`。`_Hint`的局部pw增加`an`参数，唯一log变为`log(analytic=an)`，H、Ha、Hb、Haa、Hab、Hbb_int六种回调各3个pw调用，共18处全部传自身`an`。总增量为4个log关键字、18个pw实参、1个pw形参，无遗漏回调。

此修改满足[python-flint积分回调的解析性要求](https://python-flint.readthedocs.io/en/latest/acb.html#flint.acb.integral)：复杂球触及log分支时不再默认忽略解析性检查。原积分区间、常数Ra、幂次、符号、2倍归一化、Hbb边界项及Hessian包围式全部保持；未换成新核机制或另一种生成元。标准库补丁生成/AST检查只说明补丁范围，不重新证明旧数学公式或新运行结果。

## 参数域及非finite传播

外区正坐标且`r>R`时，方向积分上限`L=R-b`。`L<=0`保持原空支撑积分0；`L>0,b>0`时实路径满足`s+b>=b>0`、`a-s>=r-R>0`、`a+s>0`，没有实路径内部或端点分支问题。该事实不能替代复杂球解析性标志，补丁同时保留两者。

外区Taylor仅在`u0+v0>R`等原条件下使用。若近角坐标b=0，导数积分可能非finite，不能把它当成有限实积分、擅自设0或截断。现补丁保持原`_outside_taylor`对余项及全部B1系数的`is_finite`检查；非finite返回`(None,None)`，`sing_parts`退回原单调常数路线。若有效下界仍无法形成，NaN-safe的`lo_min/min_lo`不丢弃该不确定项，最终`lb.lower()>=0`不能成功，递归或保存失败盒。没有为追求PASS添加任何有限替代。

这只审原有路径保持不变和已定位契约得到修复，不声称重新独立证明全部单调/Hessian公式。内区I0/I1/I2的2F1表达没有改动。

## 唯一运行接口及证据边界

wrapper运行前核对5项approved hash和准确B1 hash；以修正driver目录加载同目录修正hlcert，固定`--only MKA --Q 8 --maxdepth 7 --Kdepth 6 --procs 6`。原96bit精度不变。已有meta/log/report即拒绝，独占新session进程组；超时杀该组并等待，正常终态也清理该组。运行后记录源码/核hash、日志和存在时的report hash。

180秒从短manifest/hash准备结束后、Popen之前起算，覆盖checker启动、导入、全部计算及退出，不应称整个turn或准备阶段总硬限。没有失败重试、换核、加深或搜索接口。

旧first历史结果和源码hash不改写；只有本次修正程序实际闭合后，才能把它计入已修复的first全域证据。失败下界仍是充分包围不闭合，不自动成为真实负点；超时也不能借旧PASS补足。独立闭式路线及旧准确成本没有被本补丁更改，B1的中央说明沿第35轮更正口径保留。
