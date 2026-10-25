# F008 α=3/2产品与确定尾预算独立审计

固定N=32、α=3/2，不重跑下界或MC/kernel。只写本目录。先核C007的固定有理点参数尾证明，再核主删一product、商空间coercivity/ellipsoid、正face质量、时间指数与全时间可积性。

有限数值对象：①C007一组ell=1/32、R=10^18，r=64、p=7/10；②coercivity常数3及8/3的两种已提出解析封装；③ε∈{1/10000,1/1000,1/100,1/50}四个短时预算。总计不超过8标量情形；固定几何距离和Gamma常数复核不启动参数搜索。仅Python/Fraction/SymPy/Arb256，无MATLAB、无代理、无push。

验收必须指出：L1正象限中的32因子；coarea第一坐标消元Jacobian；ellipsoid的√6是否真正来自Fourier换元；whole-face对奇异点/交叠是否有遗漏；粗At质量的全时积分是否发散；新预算是否已经形成竞争性g/I。任何解析/数值改善均不自动登记C6界。
