# 第032轮：一般覆盖的联合接口与适用范围障碍

本轮三条理论均经独立审核ACCEPT，但没有新的统一Cn数值上界。原算子始终为Lebesgue、中心轴平行立方体、共同半径、一般L1及全部正尺度；没有七维专项或下界计算。

multiscale_theory/REPORT.md给出外壳盒A=(-s,-r)^q×(-h,h)^(n−q)及公共前序区的一维准确Stieltjes/壳面积积分、有限a半径优先与捕获guard。h≤1,q<n时质量≤1/[2(n−1)]；h>1时公共前序仅q个负面，正预算要求q>2τ/logt，扣减≤nlogt·2^(−2τ/logt)。这些障碍只在本单盒族及明列目标尺度下成立，不能排除R033典型非盒外壳。多cut必须针对同一个缺额测度，非负γ满足逐点Σγ(1_A/a+1_P/p)≤1才可扣Σγ；共享前序区不可免费重复支付。

joint_cover_theory/proof.md给出原有限见证V_i=aQi的准确质量LP及packing dual。K_na=sup|并Qi|/M(V)满足Cn≤K_na≤(1+1/a)^n Cn，故inf_a K_na=Cn；这是准确等价接口，不是新数值界。另构造可随真实见证质量取预算的覆盖flow，把所选及拒绝体积合入同一负载；原实例显式映射可行，有限LP与knapsack dual方向均已证明。所有见证margin可同时趋于1，不能仅凭正margin获得统一改善。两个最小n1/n2例的有理primal/dual/flow已经实际检查，非高维扫描。

anisotropic_cover/proof.md只将辅助资格区改为坐标乘积，原见证仍为共同半径立方体。完整有限a exact-sup公式合法；c_i>1时准确角点公式及固定δ_min>0再令a∞可积收敛均已证明，禁止无控制交换无穷尾。包括c_i≤1的有限a模型预算不优于同δ极限。此族最小值=nlogn−nloglogn+O(n)，leading1障碍只限此族。尺度重权若用全尺度统一比值K±支付，公式变为(sqrt(AK+K−)+sqrt B)^2，不能改善常权；非零幂权某比值无穷，需新证明。无害端点δ_i=1,rho=1的面数应a.e理解，不影响积分。

独立审核汇总../../audit/round032/general_cover/README.md。主逐段核对原证明、界方向与全尺度桥接；实际代数/Fraction/SymPy/最小LP输出保留各verification文件。原件及来源hash在MANIFEST.json。没有将有限样本视为全域证明。公开比较沿已核原始来源，本轮为内部接口与障碍，不称发表或真实Cn下界。

本轮还保存后续q=1/20阈值的独立检查q005_threshold.md，仅与R034定理结合时可推Cn<T−n/(40sqrt L)，L≥2^22；本轮自身不提前宣称该界。下一轮转向变量P(w)和典型外壳。下界暂停，七维专项停止，goal active，旧论文未重编。
