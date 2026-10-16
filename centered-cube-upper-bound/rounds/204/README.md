# 第204轮：共同捕获容量与加权余权极值

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮独审通过：原捕获量的有界原函数与有限变化公式、同一 α 的共享来源容量、真实权重下平方余权的精确二矩包络。

令 d=1−k/a、ζ=αJ/a、ρ=(d−ζ)+，原函数 Ψ=d+²−ρ²。normalization 与原 capture/amplitude 导数联合结算后，α/Hcap 导数精确消去；实际 k/a、硬面、门及势端态费用保留。同一来源满足 ∫Uν+(1/2)∫aΨν≤Dq，原 Π 权转换仍含 η mathcalB/(ea)，不能直接用无权容量支付有权交通。

在同一真实 H_gate Π² 测度下，S=∫mathcalB d+²，T=∫_(d>0)mathcalB dζ，U=∫_(d>0)mathcalB ζ²，W=∫mathcalB ρ²。有限矩时严格有：U≤T 则 W≤S−2T+U；U≥T>0 则 W≤S−T²/U。抽象非负比例的二矩类内界精确达到，证明使用逐点多项式支配与匹配分布。真实 active 门若已限 ρ>0，只能使用低支；扩展零余权槽须同时重算全部矩。固定 τ≥1 的 Rτ=min(ζ/d,τ) 保持 W 不变，只需 S 有限便可使用有限截后矩。近极值结论只约束截后比例，不能反推原远尾位置。

尚缺：真实加权势负向净流、current k/a 界面、资格入口与重入及 selector 的联合全局预算。实际矩和权重层分布未获允许来源费用；所有 β、j、时间、窗口、源依赖帽及原 n 阶能量保留。抽象极值分布未证明可由原中心立方体实现；静态容量不控制路径正变差。

主代理实际运行最终 reproduce.py：13 文件哈希通过；22,320 个 Fraction 有理分布及 66,960 个截断配置通过，含三阶段扩展、精确亏损、稳定性与无效 d=0 槽回代；另重放 9 个有理路径、3 个同来源配置、2 个冻结系数真实 φ 积分例（25/40 位）。两份完整 JSON 均逐字节一致。数值枚举不是一般证明，高精度求积不是区间证书，冻结系数例不认证原 FULL 几何。923 行详细 LaTeX 全文独审及静态结构检查通过，未编译 PDF。

[research.zip](research.zip)：65,705 字节，14 成员，SHA256 `3bde5df51212331867a39031490db60d5473a8e04b0ac76f062e078237b9b90a`。解压后 `python3 reproduce.py --verify-only`；完整复现 `python3 reproduce.py --output-dir <新目录> --numeric-python <具有mpmath的Python>`。本机实际使用 Python fractions 与 mpmath1.3.0；PYTHONPATH 为 `/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor`。完整环境与有限检验范围见 numeric_notes.md。

展开材料已在仓库外 final_checkpoint 备份并逐字节核验，14 个归档成员也全部核对。只新增本说明及归档，更新项目索引，无 PDF 或缓存。复用 1 Astra high、2 Sol high、1 Sol medium，使用 B03 与 N04 数学工具技能；缺少逐模型实际扣额数据。归档不足 66 KB；不改写历史。
