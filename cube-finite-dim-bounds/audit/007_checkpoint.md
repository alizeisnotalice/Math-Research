# 第5轮首checkpoint：独立准确审计通过

2026-10-25。输入SHA256 `5ebe44bd14b85cb8483ac648a690f9226a07da35081eda848cd0ca38785f8cd7`；P1准确改为1.630361663339249，其他参数的有理值及170个key逐项保留。字符串尾部补零未当成数学参数变化。

准确目标是前轮T3完整上端，来自audit/006-satellite_n3_results.json，结果文件SHA256 `51c9a979ea5bfdf9e9e110cbd8ebed107c305859e97a2a75332015941dfcc47d`，不是显示小数或旧R3。本次复用先前独立Newton开根+slab-bitmask程序，但重新枚举新坐标全部3009个原子实例、核170个内外盒。未导入作者root/volume，未重算旧并集。

独立准确下端L5满足

    L5 > 1.7939726890816552，
    L5 − T3_upper > 0.00000008467763。

完整有理上下端及正差在007_results.json；与作者准确区间相交。9个半径临界等式、b=2、索引最小(-1,-2,-1)、最大(2,2,2)、质量M=500280708211/500000000000全部保持。每盒最小边长仍严格正，故006的统一收缩未裁剪盒、半径rcrit−τ/2、同质量平滑η<τ/2，以及有限索引块N→∞、τ→0传递原样适用；无需重验旧证书或把临界等式误作严格阈值。

首次运行仅因源JSON对未变参数补零，字符串相等断言失败；改为Fraction有理值比较后通过。最终一轮约0.052秒，30个小整数网格自检全部通过。Checker SHA256 `7153b6afaa69d7715064e9a203a26f87a516899464ddb1b9a6cd09f65423694f`；同006共享独立实现和Python整数/Fraction基础，明确不声称新的无共同代码认证。

MATLAB实际另以完整联合170盒的轴向布尔覆盖矩阵求体积，得到体积1.7949798548099192，质量1.000561416422，比值1.7939726890816492。此为浮点复核，不承担准确改进比较；变量aud7前缀，未clear、未改path或工作目录。

只读复核（从已整合输入的仓库根目录）：

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/007_check.py --author-result cube-finite-dim-bounds/lower/n3/round005_certificate_result.json

默认输入为项目lower/n3/round005_certificate_input.json，比较目标默认同目录audit/006-satellite_n3_results.json；可用--input、--seed-input、--target-certificate显式指定。默认stdout-only，加--save才写本批结果。

这是首checkpoint，不更新账本。主代理已授权作者在新预登记下继续有限refinement；后续不同不可变输入须另记子批、重新准确比较。当前核验未宣称参数族最优或完整周期超水平集体积。
