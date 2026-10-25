# 第007轮精简交付：固定旋转格的完整 level1 检查

精确C间距、旋转cos=4/5、sin=3/5时，完整全尺度超水平集比值被准确包围为

    [1.73995838416010979469265955500564837669736363224817715476779,
     1.73995838418300349717663336813054885976344568823264611332333]。

该固定输入在level1不能超过完整C_cert。这里没有一般算子上界、其他角度、
其他间距或其他threshold的排除。007交付不推广正式下界。

文件职责：

- `model_and_completeness.md`：超胞、去重、平移覆盖、8/25归一化及全尺度完整性。
- `rotation_atoms.json`：唯一225原子邻胞输入。
- `geometry_bundle.json`：26个family的准确整数网格partition与固定原子ID。
- `case_evidence.json`：256个原始K/n案的准确结果、零案证据与真实执行来源。
- `aggregate_histograms.json`：四个必要聚合面积/体积直方图；68个可重建分块不归档。
- `result.json`：完整有理端点、参考C哈希及核验来源。
- `profile_setup_v4.json`、`radius_certificate.json`：22个半径区间的准确排除与真实PID/RSS。
- `replay.py`：可携证书算术、单块重放及可恢复组合入口。
- `profile.py`：每次仅重放一个半径区间。
- `capture_case.py`：按需重建一个原始K/n案，统一处理全部原始shape。

所有计算入口强制nice>=10、一个Python线程及六个数值库线程环境变量为1。
无NumPy、MATLAB或额外依赖。下面前缀用于每次单独调用：

    env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 \
      nice -n 10 python3 replay.py verify

该命令只重放已保存证书算术，不重新计算旧完整union。更深入检查按有限块调用：

    ... nice -n 10 python3 replay.py block plane 0 0
    ... nice -n 10 python3 replay.py block unrot 0 0
    ... nice -n 10 python3 profile.py 0
    ... nice -n 10 python3 capture_case.py 4 1

`...`表示相同六个线程环境变量。每次仅一个命令、一个进程、一个有限块。
输出写入`replay_checkpoints/`或新案例文件，不覆盖已锁定证书。每个plane块最多
36个x段，每个unrot块最多16段；完成全部块后可用`replay.py rollup-blocks`
核对组合结果。索引缓存、全部阶段报告及重复零partition保留在本地研究目录，
不进入本精简归档。

先前K1/K8、K2 append检查点及其manifest均未改写。源程序仅作路径与入口整合；
`case_evidence.json`明确保留真实执行时的来源hash，不将整合源冒充原运行源。
