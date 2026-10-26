# audit_transition062：冻结 052→061 有限树转换独立审计

结论：**指定的 061 batch0001 转换通过本轮独立核验**。2048 个新父的 165888 个点均有严格正的准确有理余量；旧、新前沿在全部 52360 个 root 上分别 prefix-free 且 Kraft 和准确等于 1。退休集合恰是被新父覆盖的旧记录，全部缓存和全部保留记录逐字段不变。

这是有限冻结状态转换结论。queue 仍有 61478 项；本轮不推出全域 g 或 C5 结论，不授权或启动新数值链，终审由主任务负责。

## 冻结输入与执行

输入均在运行前后核 SHA-256 相等；没有读取后续活动 state。

|输入|SHA-256|
|---|---|
|052 `batch0043/state.json`（743364738 bytes）|`065382625da1645e267025ef7c26d292a29a62b65cb7ab97f49dfa9c06995dbf`|
|061 `batch0001/state.json`（735202818 bytes）|`d0d62d99d57b3f15d3acb1f4a635b188e196b0de6e155b461460ccbd64a72d31`|
|061 `certificates.json`|`a3331af02c0a3df5548977971b85e62f2dac3b39f591658289569ab4f402b151`|
|061 `retired.json`|`15a1196e32d8714d05211da7ff95697509bcf72524c98977ff09e96503ee1e64`|
|061 `registration.json`|`7058a5e4b981fdc08bc9393aa7b5d40e8240a5af81d4938b4d4b64666254681f`|
|061 `manifest.json`|`2e1c6a69d9fd55ff7c09d65cd356c88bf726151f07744c6b7311c4d43cda3504`|

052 绝对路径为 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover052/batch0043/state.json`，061 输入同位于 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover061/batch0001/`。

最终实际核验进程 PID **59356**，单进程、单线程，exit 0。总时间 **7.113759416 秒**，双流结构与对象比较 5.493187666 秒，2048 父整数证书循环 **0.530967667 秒**，峰值 RSS **148488192 bytes**。核验循环限 2048 父/120 秒，支持已验证前缀恢复；本批未触及时间限。

本检查器新 H、gradient、kernel、Green、log、exp、Gamma、Arb 调用及 MATLAB 调用均为 **0**。旧、新状态的 `total_kernel_attempts_completed` 和 `total_new_anchor_starts` 均准确保持 178722；这只证明保存计数不变，不另充当外部进程全历史的调用跟踪。

## 直接字段比较及追溯

旧、新 JSON 同时流式解析，避免把两个约 0.74 GB 状态双重整载。比较使用解析对象的直接相等性，**不是只比计数或只比整表哈希**：

- 全部 178722 条 cache 的 key 与完整记录逐项相等。
- 80496 条保留 passed、91662 条保留 splits、61478 条保留 queue 的完整对象与相对顺序逐项相等。
- 除 passed/splits/queue 以外，旧状态其余 38 个顶层值全部直接相等；新增键恰为 `task061/fingerprint061/status061/conversion061`。
- 新 passed 的 2048 条证书与 certificates 账本直接相等；每条旧失败证书通过准确源 split 索引、node 和完整证书 SHA 锁定。
- `conversion061` 中旧 state 绝对路径、SHA、旧 fingerprint、证书/退休/注册/缓存 SHA 和父计数均逐键核对；新链启动标志为 false。

独立重建的 cache 规范 SHA 是
`f3b174e2bab3fcaed2c291090b378ddac7414393a6fb7a7bceffedca0ffaccac`。
原 `fingerprint/status` 被保留为源元数据；新状态的 `status061` 为 `finite_zero_kernel_conversion_complete`，`fingerprint061` 为 `9dd811ebb3a7682e792bf16acf09028dced820c819fdd020eb992ed22ae7be1b`。本轮核了新增 fingerprint 的保存值及旧状态链接，没有复算作者新 fingerprint 的生成配方。

旧、新 `active=null`、`inflight={}`、`resume_jobs={}`。冻结 pool 收据有 close 和 join 返回、八个 worker 身份相合且全部 exitcode 0、is_alive false。本轮不把保留的旧 pool 收据解释为独立 OS 进程监视；旧 controller 已退出的外部终止核对由主任务先前完成。

## 独立物理盒与整数树覆盖

root 列表独立按 `0≤l≤k≤j≤i≤31` 枚举并与两个状态比较。对 node `(r,p)`，第 j 个坐标只读取 path 的第 j、j+4、j+8、… 位；若该子串长度为 k、二进制整数为 b，则物理坐标区间独立构造为

\[
[(r_j2^k+b)/(32\,2^k),\ (r_j2^k+b+1)/(32\,2^k)].
\]

没有导入或调用作者 box、partition、certificate 函数。每条旧 split 另核 axis 等于 path 长度 mod 4、children 恰为追加 0 与 1。新父物理盒逐端点与 certificates 账本比较，anchor 由缓存 key 回读，全部严格在重建盒内。

对每个 root，把前沿 path 映为共同深度 D 的整数 dyadic 区间 `[b·2^(D−d),(b+1)·2^(D−d)]`，排序后要求端点从 0 连续准确接到 `2^D`。该检查同时排除前缀重叠与覆盖缺口，而非仅核浮点 Kraft 和。两树最大 path 深度均为 3；52360 个 root 全部准确铺满。

另外，活 splits 集合分别准确等于前沿全部 path 的真前缀集合。新父自身及其 split 后代全部退休，真祖先继续保留；新父构成前缀无关集。每个移除的旧前沿节点有唯一新父，且是该父的严格后代。退休账本的每一项 source_index/node/unique_new_parent 与由旧状态独立产生的列表完全相等。

|记录|旧值|退休|新父增加|新值|
|---|---:|---:|---:|---:|
|passed|84802|4306|2048|82544|
|splits|93920|2258|0|91662|
|queue|61478|0|0|61478|
|前沿叶总数|146280|4306|2048|144022|

两树分别满足 `leaves=roots+active_splits`。本块全部新父是 coarse root，因此无新父的真 split 祖先；其余 root 的旧 split、叶与队列记录保持不变。

## 新父的准确 81 点证书

固定 epsilon=1/64、target=1/10000、模型 alpha1、x1=1。采用已经 056/060 审定的分闭象限下界

\[
\Phi(v)=\Delta_L+\sum_j\min(g_{j,L}v_j,g_{j,U}v_j)
 -4096\left(\|v\|_2^2/4-\|v\|_1^2/36\right).
\]

每个新父的 Delta、四个物理 log-gradient 区间与 anchor 字符串必须和对应旧 failed_certificate **原样一致**。缓存五元区间有序，H 下端严格大于 target；四个导数区间和 H 区间的所有端点组合均准确核 `g_L h≤d≤g_U h`，共 **32768 次交叉比较**，从而核安全区间除法的符号方向。Delta 的既有有向 log 下端资格沿冻结旧链继承，本轮没有重新计算 log。

对每盒检查 `{lo_j−a_j,0,hi_j−a_j}^4` 的全部 81 点。全部 anchor 严格在内，因此每盒恰有 81 个不同点。所有位移清到分母 V，Delta 与梯度清到分母 C。整数内循环计算

\[
36\Delta_iV^2+36V\sum_j\min(g_{j,L,i}z_j,g_{j,U,i}z_j)
 -4096C\left(9\sum_jz_j^2-(\sum_j|z_j|)^2\right),
\]

共同正分母为 `36 C V²`。每盒最小分子严格为正，准确最小值与保存 margin 相等；保存 worst_offset 确是某个最小点。内循环不作浮点运算，也不调用 Fraction 四则。

全部 2048 个新父通过，整体最小准确余量为

\[
\frac{2902335149548405616429378352696791399946902478226218534973017079219999399069111}
{29642774844752946028434172162224104410437116074403984394101141506025761187823616}
>\frac{97910}{1000000}>0.
\]

81 点推出整个闭盒成立的理由继承 056/060：每个闭象限内下界为凹二次式，分片矩形的全部顶点均在三值坐标集合中，公共边界也被包含。本轮核验其输入、端点、整数余量与树接入，未重新证明正 Cauchy 混合的 Hessian 定理，也未重跑旧 162 父试点。

## 文件与复现

仅保留本轮四件：registration.json、check.py、result.json（含全部 2048 条独立结果）、本报告。未修改或复制全量 state/cache；旧 audit_cutoff044 只读，无提交或推送。

|文件|SHA-256|
|---|---|
|registration.json|`1f89d0f030af64f17c204c72d320d5d51b9f5924e96385cc5b21c92f98c75cc9`|
|check.py|`566d48d2f1483e5398e9c912ce871b23c6fb5b45465a4bc5d3c10534db2fcec6`|
|result.json|`eb4dc9f228f7fb47da2f5e5c7f364e36db73ecf10300ea0ccd596645c81d62fc`|

```sh
cd /Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-audit-cont
python3 -B cube-dim-5-bounds/audit_transition062/check.py --run
```

无 `--run` 只显示登记信息；`--run --resume` 可恢复相同代码与冻结输入的已验证证书前缀。使用 Python 标准库即可，未调用作者环境、作者核函数或 MATLAB。
