# 第44轮复现

本轮只有候选，没有新的已证一般常数。候选成本约3.590862654958541，须后续全域认证。

临时根按既有scratch结构提供round01基线/成本算术、round10/enclose_cost.py、round14已冻结12方向及旧rows/模型、round24完整cover/独立报告及其历史来源。程序放round25；public independent源码版本/哈希见build源。依赖Python3.12/python-flint0.9.0、NumPy/SciPy，MATLAB本机R2026a的linprog。

空输出副本依次`python -B run_stage.py build`、build完整后`python -B run_stage.py solve`。固定300/180秒hard，只有一次LP、一次既定收缩和量化。原始证据目录拒绝重跑。

candidate_and_lp.tar.xz含candidate.json、structured_kernel.json与实际model.local.mat；archive_manifest逐件SHA及字节比对。矩阵保留了实际MATLAB输入，便于单独审阅有限LP。两份核文件用于不同核验接口；系数未独立舍入。

约221MB模型/region/rows可由build再生，不重复上传；model_manifest.json及实际build_result绑定源、计数与哈希。要审保存模型，先在另一个空目录复现build，再将匹配哈希的数据用于保存证据副本的verify_build_saved.py。该核验限90秒，不重新计算模型。

verify_candidate_saved.py限20秒，只审参数、系数、根夹逼及保存最坏行；全行和全域证明的适用范围见报告。运行必须在无对应审核输出的副本；不把复现当成原运行再次发生。
