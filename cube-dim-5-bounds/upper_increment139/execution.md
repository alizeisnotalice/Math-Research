# 139唯一授权积分批：实际执行附录

主全文审阅冻结证明、117行producer、71行checker、来源绑定脚本与固定输入后，授权唯一128面板/257根原子/192-bit批。原授权 SHA：

`82fd882a4e1c8bdb6a16157aa50b8c8a41ddae9746068fedf272d86248811bd6`

路径为 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/round016-authorizations/primary139.json`。冻结的准备期源码、proof、plan、inputs和原manifest均未修改；本附录更新其“尚未运行”的历史状态。

## 实际结果

producer PID **58186**、exit **0**、约 **0.274568秒**；随后已锁portable checker只运行一次，PID **58187**、exit **0**、状态 `PASS_portable_auxiliary_integral_certificate`。单进程，pool=null，完成128面板与257根区间，无续算、无加密或新case。

准确根端点与全部正权面板保存在原bundle中。以下小数只供显示，证明使用原准确有理端点：

- 辅助增量预算夹逼约 `[1.0128196266016944e−6, 1.0128305986937574e−6]`，区间宽度约 `1.0972092062791304e−11`。
- 新完整I候选上端约 `4.538418122604361e−5`。
- 准确Fraction比较确认：新J上端严格低于旧保存135cap；新I上端严格低于第015轮已接纳I上端，准确预算下降约 `1.8562940840435963e−8`。

上述下端仅是固定保守参数的辅助表达式下端，不是实际R/J或完整I下界。正式新界仍待140独立审计和主接纳；本分支 `global_Cn_claim=false`。

## 冻结bundle与哈希

唯一原bundle：

`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_increment139/run001/`

|文件|SHA-256|
|---|---|
|inputs.json|`fe6dac94e6de5c7fd19d1a1da2b32279945b54b6d9aae428f1c8bd42796315b7`|
|plan.json|`90c6c97da2b1d965517492a69222db2d244dc1807fdc355a67bda24ac6b1ae92`|
|registration.json|`fcf777093b794b319ef4918d4cc60cb90e74a287ded12a46136bb22405f15171`|
|certificate.json|`8ded2236301d10b829a854107fba84affa76f4c9b269da9f3b5580ee28351c19`|
|results.json|`ce244231cfc77e6cee22ce12ebc3f92d79420ed0e90ab322a210b9ee8f505534`|

外置相邻 `dispatch001/` 保存授权绑定、actual PID、退出码和原stdout/stderr。全部输入/源码/来源执行前后字节锁一致，两份锁表 SHA均为
`de712b06c16c45014f1024fbc4545cacbceac4154fcc2a13cde77ca596d97b24`。

生产和checker都没有调用H、R点值、梯度、log、exp、Gamma、范数或M优化，没有重跑135/136或旧toy，没有调用共享MATLAB。主15核H/PF作业不受这次低资源单进程短批控制。

## 复核与归档边界

140可直接读取上述冻结bundle，按独立路径核正权、根区间、完整I与源绑定；不需要重新生产这些257个根。原portable checker已经实际通过一次，后续不由本分支再跑。

`execution_receipt.json` 保存完整准确辅助J区间、完整I上端、准确成本四次幂、准确I预算下降及原文件哈希。`execution_manifest.json` 只新增本次执行metadata白名单并引用准备期manifest，不改写准备期锁。完整certificate继续保留在唯一原bundle中，供主按必要材料归档。

不再追加面板或优化，不修改冻结128/133/135/execution135，不启动下界或Astra，不提交、不推送。
