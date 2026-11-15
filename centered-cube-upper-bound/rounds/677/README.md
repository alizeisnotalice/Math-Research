# 第677轮：历史合法共享见证与全局组合费用

2026-11-15。**一般纯 C_n≤K√n 及原 n^(1/2+o(1)) 上界仍未证明。** 本轮使用1个Astra、3个GPT-6.1 Sol；未开MATLAB，未编译PDF。

## 已证构造与完整总账

保持原q、真实owner和receiver不变，在原历史次序内为全部loss重新选择共享饱和见证。每组完整正来源率为τ，后续组的全部target parents避开此前见证。重新分类并一次支付far/upper，不能沿用旧residual类别。

每个非空新剩余组存在共同low core K′，率>3τ/8−(32n)^−n；显式删除seed self后，所有来源严格早于每个目标，原admission给出整个接收盒上的来源load≤1。同一core覆盖整组，无需等对数尺度细分。局部短采样界仍保留，但跨组来源复用未付。

完整账目仍为 V₀≤C_adm[U+(K_far+K_up)U/τ+V_rem,new]、P₀≥U/τ，保留原族全部packing变量。

## 当前唯一核心费用

每个含重复标签的multi-index组合只交给最早eligible组，用原seed半径形成扩张交集H_m。确定性multinomial覆盖给 T_min^k V_rem,new/k!≤Σ(q^m/m!)|H_m|，无需Poisson成功概率损失。

**右侧全局费用尚未获得足够上界。** core首次eligible不一定是空间首次进入；改为完整P首次进入得到更大halo，其费用必须重新保留，不能混用两种版本。原q-load只控制未扩张交集，单组估计相加仍重复花费来源质量。来源标记恒等式与纯重复标签预算不替代混合key总预算。

## 实际验证

本机Python3.12.14、Fraction精确有理算术：10个真实有限cube fixture，n=1,2,3,4,6、N≤18，223个全部可实现正体积mask；每policy覆盖69个loss。Canonical得到69组，最大未分配receiver体积policy得到19组。另复用输入完成k=1,2,3的264个组合恒等式、276个去重keys（180含重复label）、806个group/key关系和1318个实际闭接收盒包含检查，全部通过。

这些实验使用低τ=κ(32n)^−n，额外维度的中心坐标复用零；key实验用full P，而非高阈值K′。它们只核验有限路由与系数，不认证R673、hard regime或根阶预算。独审先独立检查命题，再对照最终证明，八件文件哈希已绑定；JSON审查限于元数据/源算法，不称逐项人工检查或重新执行。主代理实际运行完整性及来源绑定校验；manifest与ZIP CRC通过。

## 精简材料与下一步

[research.zip](research.zip)：13个必要成员，50,127字节，SHA256 `c425c41cd57a0fb785639604444d1137e90513d8dd44caac0d859212fdad1c0a`。包含详细LaTeX、独立证明与审核、复现脚本和必要输入结果；可再生的2.2MB逐key输出仅保留本地，归档只存摘要。REPORT SHA256 `e085caeb2617be5f74882bfd0566971f00ada616c2e8d0880683c7147b43884b`。

解压后用 `/Users/zzh/.local/share/math-research-tools/venv/bin/python verify.py` 核对归档；`WITNESS_CHECK.py`、`KEY_CHECK.py`按次序复现。只有三个仓库路径改动：本摘要、压缩归档和项目入口。下一轮集中首遇扩张交集的真实联合几何费用，不重复局部阶乘或Poisson常数调优。
