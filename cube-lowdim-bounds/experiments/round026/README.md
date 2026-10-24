# 第26轮：缓存复用和六次MATLAB约束生成

**没有新证书或新弱型上界。** 六次LP达到事前上限，末解仍有两个AK数值违反点；准确已证基准仍为B1<3.59091，尚未证明C₂<3.4。不能把数值成本3.5903930148360588当作一般L¹界。

一次实际运行44.889370333秒。18个缓存文件hash通过，0次新B/L/奇异闭式求值；种子与六解共361018个完整池残差，新细方向4951104个条目。单独MATLAB batch六次linprog均exitflag1，216→288变量，累计最多96细类；没有basis热启动声明。末解K0个违反，AK2个违反，最坏−0.0003370508029512109。原活跃约束与完整池复查分别记录，有限全池不等于连续域认证。

六次求解均完成后，独占MATLAB退出等待超5秒，协调器TERM其进程，外包装清理同属本轮的残余进程组；不是在途LP超时。worker_done记录6/6，主侧之后signal0确认本轮PGID不存在，未操作共享MATLAB GUI。

## 最小复现与依赖

Python3.12/numpy2.5.4/scipy1.18.1/highspy1.15.1；MATLAB R2026a Update4+Optimization Toolbox，已实际路径/Applications/MATLAB_R2026a.app/bin/matlab。新建同结构scratch目录round01/round04/round05/round06/round07：前几轮基线、诊断点、细化恒等式、公式模块按各轮README放置；复现第25轮Python步骤会产生本轮所需cache.local/point_manifest/result/初始行清单；MATLAB种子必须使用第25轮归档的原始matlab_result.json，其固定hash不可由重新求解后的带新计时JSON替代。缓存必须与第25轮已记录hash吻合，缺失则本轮停止，不自动重建。

round07放search.py/batch_solver.m/run_once.py，使用`python -B run_once.py`。驱动pin第25轮源码与MATLAB种子hash，并核全部缓存、点序、列序。已有真实运行不需重复。本轮enclose_cost.py只作预备：它从第22轮独立整数/Fraction算法把目标字段改为17/5，未有候选所以未执行成本或全域核验。

## 保存范围

result.json保留每次解、完整残差摘要、输入/结果SHA、变量/行数与计数；worker起止、batch日志和cleanup观察用于区分求解完成与进程退出清理。local_artifact_manifest.json是终态后生成的本地数据hash清单，不是运行前证明或核证书。所有矩阵、全点值、系数快照和MAT约582MB留本地，不上传；下一轮复用最终288方向缓存，不能擅自续用本轮剩余预算增加第七次LP。既有公开公式与已证基线不重复归档，私人原件保持本地。
