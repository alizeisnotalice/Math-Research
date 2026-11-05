# R073：二维真实宏微九阶下界

[完整实际构造](../../proofs/073-planar-macro-micro-ninth.md)和[独立重建](../../reviews/073-planar-macro-micro-ninth.md)给指定有限族的绝对九阶展开，不确定C₂九阶极限或最优微孔。每个充分小实u取有限U4、dyadic微网h≤u²⁰、N=ceil(u⁻¹²)，所有原长短见证及严格裕量在证明中明确。

实际本机算核：`wolframscript -file ninth_budget.wl`（PID81523）及独立标准库`python3 ninth_fraction.py`（PID81866），均从同一个真实质量比的有限有理式展开到u¹¹。WL另核任意高阶剪切系数消失，低阶和Jhat回归残差为零。`host_leading.wl/log`是独审另行实际运行的首帽矩/密度/预算核（PID82033），不承担连续几何证明。

纯宏四项面积所需的原精确12矩积分已在[第67轮I2包](../067-iterated-protection/README.md)归档；本轮不复制或重跑它。所有帽的一致微格误差、有限偶层正部、共同host/总M、外盒全部半径及O(av^(19/2))由解析推导承担。

独审原工具snippet只存JSON的问题已补齐：`micro_increment_fraction_recheck.py`导出相同Fraction核心并实际新跑（PID86148、exit0），完整payload等于原JSON；旧输出和冻结终审未动。详`reproducibility_supplement.md`。复现此核请将该py和`micro_increment_fraction.json`复制至新临时目录运行，它用exclusive creation保护现有输出。其新JSON与原stdout逐字相同，只归档一份。失败的首次WL回归因Coefficient与标识符写法问题保存在本地，不作为成功证据；最终重跑文件和原日志在此。

所有脚本依赖仅Wolfram15或Python标准库，不依赖缺失私有计算模块。`manifest.json`绑定实际程序/输出与冻结源。数值/代数核不能证明所有输入的九阶上界。
