# 第225轮：移动尺度硬边界责任与共同来源预算

一般 n^(1/2+o(1)) 上界仍未闭合，本轮没有新增根阶付款。依据新建议并结合第224轮已有进展，本轮聚焦实际共同来源几何，未重复支付固定J分配。

新估计：完整混合核满足尺度导数 P(T−H)dx+B，其中中心化平滑score有二阶矩≤nsIφ，Iφ≤805/27；hard责任与正边界通量均为n(1−s)。固定s时，尺度列价≤1+log(b/a)[n(1−s)+sqrt(nsIφ)/2]。任意可测实际s(x)、A(x)下保留准确的选中硬壳层减体积补偿；固定s中心化不可免费用于自适应s。

原两移动腿保全部Π、FIRST、η和唯一V/j，得到F1+F2≤sqrt(2dXB_d)，原分配价≤εX+[512Λ²K/(27ε)]d e^d B_d。仍需证明实际共同占用的B_d≤(Cd/K)log^A(n+2)Wν；该式明确未证。删除来源关系的完整回退仅给n²log²n级来源费，不登记为预算推进。更深clock也未付。

## 实际核验

- 240行LaTeX文本、两个作者稿与独立审查归档；无PDF编译。
- 本机Python/mpmath三阶段真实φ核：窗口24/48/96、25/35/50位、GL48/96/192并细分边界尖点；I≈0.99016508953918365。普通求积误差与解析尾界分开，不作区间认证。
- 512/1024/4096路径、n至128的有限模拟；不拟合维数阶。另252项精确硬面selector恒等式、7个C_n有理值及常数检查通过。
- 独审先筛命题后读证明，核查所有费用和未证标签；精确脚本独立运行完整JSON一致。
- 最终包外迁移重跑：numeric_probe约24.957秒、exact_checks约0.038秒，两份完整JSON逐字节一致，12个manifest项前后通过。回执在本机 `output/centered_cube_round225_20261016_driver_replay/receipt.json`。

## 精简归档

[research.zip](research.zip)：13成员、48,924字节，SHA256 `8fe30f2ef30ebecbe266507b48e9c6ff8c8f545645753ac099feb124dce9a3ad`。本轮新增2个仓库文件，另更新项目README。

展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round225_20261016_final_checkpoint`；独立ZIP：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round225_20261016_research.zip`。与本轮work逐成员核对后才移除work及派生缓存；原材料未删除。

复现：解压后 `python3 reproduce.py --verify-only`；完整复现 `python3 reproduce.py --output-dir <新的包外目录>`，可用 `--mpmath-path` 指定已有mpmath1.3.0。不需要MATLAB或新安装依赖。
