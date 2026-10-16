# 第205轮：共同交叠面与捕获恢复的尖锐几何

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮独审通过三项可用于原有限变化的新估计。

1. 任意维度、固定来源对、全部软硬参数 s∈[0,1]、同比例 target 尺度 γL≥L，真实双密度相对硬面核满足 sum_faces ∫ |(v−y)·normal| P_s,L ΔP_s,γL dH ≤ (1+16e/9)L^(-n)。两立方体交叠的尖锐极值与其它坐标的平均交叠衰减共同消除逐面求和的 n 因子。尺度归一化、接触迹和全部原权重保留。
2. 原截帽分母有限跳有尖锐包络：F J²[(q+J)^(-2)−(q+J+δ)^(-2)] ≤ BP[1−(P/(P+Bδ))²] ≤ min(BP,2B²δ)。这对有限跳成立并消去 J，但原 B=B0/Wg 与真实非线性面权 K_B 仍需收费。负分母下降必须保实际值，不能错误地用该正上界替换负补偿。
3. 实际同 α 的共同捕获增益 G 给 rho+²−rho−² ≤ d+² G/(H−+H+)。持续 target 保实际分母负项；仅增益路径的系数和以半个储备对数比控制。储备补充、零后重入、资格新入口和几何交换完整保留。

另证 clamp 面在明确乘 L^n 的共同尺度积分下，n≥2 的费用≤2，n=1≤log2；无此权时仍可能有 n 费用，未将归一估计冒充原付款。所有定理均有完整证明及独立审查，不依赖数值拟合。

主要缺口：真实共同面上的 K_B 加权占用，与捕获新入口、储备补充和 current 几何交换，尚未在原 βΠ²、j、time、window 及 receiver selector 下统一支付。源依赖倒 Wg、无权 clamp/continuous 项和原 n 阶能量的回代限制未消失。固定参数的全 s 核工具不等于一般弱型闭合。

最终主代理实际运行 driver，耗时3.802秒：13文件哈希全部通过，两份完整JSON逐字节重放一致。检验含2079个generic有理cube配置、30个高维构造配置、2274个明确排除的接触配置；恢复/共同来源/分母跳/二项式精确组；6个实际mixed-product核的25/40位分段积分；另9个有理clamp积分与3个n=1非区间诊断。最大mixed跨精度绝对差<9.02e-24。有限证据不代替一般证明，解析尾界不包围求积舍入，不认证FULL物理输入。860行详细LaTeX全文终审及最终静态结构通过，未编译PDF。

[research.zip](research.zip)：84,325字节，14成员，SHA256 `2f25614892bf0c5e999b5abc8dde2838c913290ec466a9ea3b54deab5eeb44fb`。解压后 `python3 reproduce.py --verify-only`；完整复现 `python3 reproduce.py --output-dir <新目录> --numeric-python <具有mpmath的Python>`。使用现有Python Fraction与mpmath1.3.0；本机PYTHONPATH为 `/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor`，更多范围与依赖见numeric_notes.md。跨环境JSON差异须核原因。

14成员及仓库外展开备份已逐字节核验，才移除本轮自身work。仅新增说明和归档，更新项目索引；无PDF、缓存或历史重写。复用1 Astra high、2 Sol high、1 Sol medium，沿用B03数学工具的逐点支配纪律；无逐模型实际扣额数据。归档不足85KB；提交前实查项目根目录821项、rounds下204项，Git loose对象112.19MiB、pack261.91MiB。本轮只改3个仓库文件，不改写历史。
