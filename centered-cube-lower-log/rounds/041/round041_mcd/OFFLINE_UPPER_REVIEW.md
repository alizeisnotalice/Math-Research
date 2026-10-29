# accepted replay的独立离线上包复核

结论：accepted receipt与全量独立公式核对通过。没有导入upper_math/source/driver，没有新样本或核计算，也未重新计算统计CI。普通FP与原kernel/BE/r0前提保持。

root通知RECONSTRUCTION_RECEIPT accepted之后，才实际运行 `offline_upper_check.py`。脚本独立采用log-odds/logistic Cantelli公式及直接McD指数，读取保存的logmean upper、logvariance、S代价及旧40payload，覆盖4R×3模型×256样本×233尺度×21阈值，共15031296行。

|复算对象|与保存贡献数组的最大绝对差|
|:---|---:|
|旧upper|2.047e−12|
|variance/3-only upper|1.365e−12|
|per-scale min旧/variance3/McD后union upper|1.365e−12|

对所有保存行，old lower≤φlower≤new upper≤variance-only upper≤old upper及原贡献帽均无正违反。独立恢复conditional B最大约.9653，未越[0,1]；h与logD和旧40数组逐位相同。零S、strict等号的保守返回1约定也与冻结实现一致。

不仅信任receipt的布尔值，脚本还直接读取old40和base_replay，核对每R的count-stream、sampled-anchor、source-statistics三hash相同；所有旧NPZ键及数组逐位相同。upper-input文件hash与accepted receipt一致，旧RESULTS/receipt、replay plan/observer/upper_math/MCD证明hash链均通过。记录的最大source r约2.512e−7，低于预设r0=10^−6。

实际独立核验耗时.593秒，Python3.15/NumPy2.5.4；accepted重建为root已有的同数据协议重放，不是独立MC。本审不改变原alpha或任何结果。详细可复现记录在 `OFFLINE_UPPER_AUDIT.json`。
