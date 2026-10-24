# 第36轮：修正first checker解析性契约后，仅复核已证B1

状态：静态准备，未执行数学。主侧下一turn审读冻结后负责唯一180秒hard全域运行；本代理只制作最小修正副本及可审diff，不改upstream，不改历史结果、不运行候选或旧失败核。

目标是补齐旧first checker外区acb.integral回调没有传播analytic标志的接口缺口。复数log在analytic=True时必须遵守解析性契约；实路径为正且避开奇点也不能免除复杂球的此要求。旧历史材料保留并由中央ledger区分证据强度，本轮不追改其SHA或宣称旧运行自动获得修正。

## 唯一文本修正

固定原hlcert.py SHA256 `4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5`，原cert_run.py SHA256 `079906db354febdd57d537f0d747428bc90bd37fa88297449cd1d46c7446310d`。只在round17/fixed_checker创建两文件：

1. `Apsi_outside_point`回调里的三个`.log()`均改为`.log(analytic=analytic)`。
2. `_Hint`的局部pw(z,e)增加an参数，唯一log改为`.log(analytic=an)`；H、Ha、Hb、Haa、Hab、Hbb_int六种回调的全部18个pw调用均传自身an。
3. cert_run.py逐字节复制不改。其原精度96bit、默认设置、各公式、阈值、舍入/非finite处理、积分限及作者归一化全部不动。

make_patch.py只用标准库读写、hash、字符串替换和AST，不导入flint/checker、不求模型或函数值。它pin原SHA、限定两段替换范围、保存精确unified diff为hlcert.patch和原新SHA/计数为patch_manifest.json。另将新AST只擦去新增analytic关键字及an参数后与原AST准确比较，确认除此之外结构相同。固定计数：4个log关键字、18个pw末参数、1个pw形参。该静态AST检查不证明原数学公式，只审补丁范围。

已有副本/manifest拒绝覆盖，静态准备与将来的数学运行分开记录。root拥有run_once.py/approved_inputs.json，audit拥有独立审读文件，不由本代理修改。

## 预定唯一数学运行

唯一核为第一轮准确B1：`../round01/scaled_certificate.json`，SHA256 `21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`。主侧将用修正目录的cert_run.py及其同目录hlcert.py，固定Q8、A maxdepth7、Kdepth6、procs6、onlyMKA，180秒进程组hard，包含导入/全部核验/保存。示意命令（实际以主侧冻结wrapper为准）：

`python -I fixed_checker/cert_run.py ../round01/scaled_certificate.json --Q 8 --maxdepth 7 --Kdepth 6 --procs 6 --only MKA --out corrected_first_report.json`

只能运行一次，不加入其他核、不调精度/深度、无失败重试或额外点诊断。记录原新源码、核、命令、实际耗时、退出码、日志/report SHA及所有原有计数。失败下界不是严格负点；超时或不闭合如实报告，不以“旧已过”替代新结果。非finite维持原fail-closed逻辑，不截断成有限数或0。

若本轮完整通过，只表示修正后的first实现对B1已闭合，结合此前独立second及既有解析传递接口加强该已证上界的复核链；不产生新成本或更小C2。若失败，则该缺口仍未补齐，不能换核/放宽设置来追求通过。用户目标仍3.4，本轮只处理验证工具的严格性。
