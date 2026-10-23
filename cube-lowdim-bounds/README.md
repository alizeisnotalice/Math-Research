# 低维中心方体弱型常数

研究Lebesgue测度、中心、轴平行、所有坐标共用同一半径的方体极大算子。按有限轮次持续研究；不以精确最优常数作为阶段验收要求。

第1—2轮已完成。二维公开下界Phi约1.685509993355525；内部三、四维下界分别为R3约1.79393322779477与R4约1.90174261110095（既有构造，非本轮新结果）。公开基准二维上界准确值2787954611/718800000；第2轮参数证书比此减少1/20000000。三、四维使用一般上界8、16，未宣称已穷尽最新文献。

证明见 proofs/，每轮报告见 rounds/，标量与几何证书见 results/。解析上界依赖Lin的通用障碍定理；新参数未做Lean构建。内部原始附件与完整日志仅本地保存。

复现（项目根目录）：
```sh
python3 baseline/verify_bounds.py
python3 baseline/algebraic_bounds.py
python3 src/verify_upper.py 1
python3 src/audit_stage.py
# 数值发现需要先安装 requirements.txt
python3 src/search_stage.py 1
```
已有审计只在改动或新疑点时重跑。搜索输出是数值证据，不代替有理证书。
