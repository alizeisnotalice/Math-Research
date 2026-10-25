# F007：混合误差与时间band空间外包

首批20分钟，零新kernel，最多三个解析测试：①独立精确D12及双轴预算门槛；②正函数的纯mixed Gauss反例；③准确Cauchy平移比与稀疏gauge包络。只写audit/round007，不改006，不调用MATLAB、不扩招、不commit/push/stage。

输入：B clone upper/round007/mixed_time_error/report.md；A本轮发来的band/J2/Taylor精确接口。仅使用Python/Fraction/SymPy/Arb256。A实际kernel和梯度由A计算，F只核证明与读取结果。新增空间包络不新增核点。

验收：全张量恒等式与总阶D_m正确；完整31项缺项如实登记；A局部半径必须由全box J2和严格中心梯度支持；新平移不等式提供完整证明及独立标量核验。不因局部盒改善或算法接口产生C6晋级。
