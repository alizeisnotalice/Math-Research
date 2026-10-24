# node663：保存值的唯一16行Fraction诊断

事前范围：只读本轮node663、来源cover及B1，最多16行各一次准确Fraction分解，总60秒hard，先保存本计划再唯一运行diagnose_saved_once.py。已有failure_breakdown.json则拒绝重复。0新AS/Psi/Hint/积分/点/细分/T/LP，不改原运行文件。

核来源SHA，检查local/physical几何与跨R情况；读取已存F0/Fu/Fv准确球、两个Hab球及绝对上端，准确核Muv、du/dv和单侧D=Muv*du*dv。只报告这些保存球的宽度，不重新评价曲率或中心函数。

每项保存T1和AS Bernstein上下端/宽度、扣D与ES前的联合下端、c1D、ES及最终准确下界，核对原16行一致。汇总最坏索引、负行数、空间系数跨度和数值球宽，区分单侧解析扣项与数值误差。float仅显示近似，判号用Fraction。

本次不做4象限限制或64行试验；它们只能作为之后独立预算的实质模型路线。负充分下界不等于真实AK负值。输出诊断JSON及机制说明，不启动T方向。
