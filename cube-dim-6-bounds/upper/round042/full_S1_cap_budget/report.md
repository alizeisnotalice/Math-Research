# C042 完整正时间 S1 帽预算

本批闭合原body32 Cauchy模型在整个正face S1帽
`C={x=(1,y),y∈[0,1]^5:dist(x,S1)<1/64}`
上的**完整Green范数上预算**：

    ||G||_(L6/5(C)) <2118205441/10^15
                       =2.118205441×10^-6.

它同时支付选定大参数分支P的短时／完整长尾，以及补事件B的整个正时间轴。没有帽外区域、全face I、g、Hnorm或C6新证书。新定理和有限计算仍由主协调者保留最终证明及集成责任。

## 新增内容与完整账目

A041/F041已证明真实正分解p_t=P_t+B_t：P为同一Cauchy参数`|a_ones|≥1/2`，B为其补集，其他31个参数在两分支中都不受限制。A041原标量只支付`0<t≤.001`，本批利用其已审ALL-t点态工具补完整时间，而没有复跑原短时预算。

独立推导见 `theorem.md`：

1. P短时由真实L1切片与单个Q5峰插值得到，不能假定P自身有反射对称性。
2. P完整长尾的原函数准确为
   `R(ε)=4/(3ε³)−16/ε+32 atan(1/(2ε))`。
3. 只选一次解析松式最优值`ε*=(Cbig πV)^(1/5)`，其中`Cbig=R1/(8π6)`，`V=(4√6π²/45)64^-5`。未优化准确尾式或扫描参数。
4. B分支全时积分以`a_i=c_i/√A_i`化为三个二次因子积分，其准确值为
   `π(a+b+c)/[2abc(a+b)(a+c)(b+c)ΠA_i]`。
   符号由二因子部分分式独立推导；正双重／三重根通过微分身份和统一可积支配连续延拓，不需除根差。

参数事件P/B disjoint且exhaustive；P的两个时间段disjoint且exhaustive；B的一次完整积分包含两个时间段。两事件产生的空间函数并不互斥，Minkowski正向支付它们的和。S1帽内所有S2交叠均已由p=P+B支付，不另乘31或32。以后合并其他空间预算，须用正union bound或显式减去帽；本轮没有宣布空间管互斥。

## 门控与实际运行

先独立写证明和门控代码、登记哈希，等待时不执行任何标量批。预登记input SHA字段为
`6b106814e735e9e75a7f2ff20f40630d0b07de61ef344acf46797f50a131a06e`。
收到主代理正式release后，核对：

- A041 manifest SHA：`a1c195da4acad86e1ed9804458b1c27f35afb2cdf788f525c3434af659fff21d`。
- F041报告 SHA：`ef114d32c209b783dae1572d4087689733fcc3e37984c37bda29dc1855547339`。
- F041 manifest SHA：`a93aeaee9352f420c43a315de02fe804c4678ea42a5d4119825f281cf8b9c7ef`。

F041接受ALL-t两分支及pointwise_coarea_note；其input_manifest锁定与本批完全同一A041冻结文件。release后才运行唯一批次。等待准备474.5969235897064秒；正式批内部0.20584408400463872秒，release后进程墙钟0.2847858330060262秒，退出0。

严格六配置：三因子积分形式身份；双／三重根形式身份；长尾及解析optimizer形式身份；唯一帽常数和ε；P短时加完整长尾；B完整时间和总帽。前三配置的全部SymPy残差准确为0。后三配置使用Arb640；只有rho=1/64与唯一解析ε值，零扫描、零核调用。没有重新枚举旧frame/geometry/triples，也没有调用037或026计算器。

## 已支付的数值上预算

以下为便于阅读的近似显示；准确证书是 `results.json` 的完整dyadic mantissa/exponent。所有区间向外取整，不从短radius展示重造。

|项|范数上预算，展示值|
|---|---:|
|唯一ε*|0.0004422980819045475423|
|P在(0,ε*)|1.5035329609206304e-6|
|P在[ε*,∞)|5.8470589186795801e-7|
|B在(0,∞)|2.9966587932972498e-8|
|完整帽总预算|2.1182054407215609e-6|

`budget_ledger.json`还用Fraction把三个保存分支区间的**上端**准确相加，再向上取整为`2118205441/10^15`。这个上端不依赖显示小数或相关变量误差抵消。因此帽积分贡献满足

    ∫_C G^(6/5) <(2118205441/10^15)^(6/5).

它仅约占条件性whole-face阈值`4/48750`的2.58%，但其他空间区域完全未由此支付；不把小帽预算当作全I改善。原R_c路线的独立`g≥1/65000`和`Hnorm≤1/14000`仍未证，C6已认证上端保持64。

## 复核与封存

`evaluate.py`运行前要求`dependency_release.json`且逐一检查冻结输入SHA；已有结果时拒绝重复运行。`execution_receipt.json`保存唯一命令、源hash、等待／计算时间和退出。`verify_saved.py`纯只读核对来源、门控、六配置收据、准确正端点和完整事件／时间账，不重复SymPy或Arb配置，不调用核、不写文件：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round042/full_S1_cap_budget/verify_saved.py
```

使用已有Python3.15.0、python-flint0.9.0、SymPy1.14.0、FLINT3.6.0；数学运算为准确Fraction/SymPy或Arb640。库及最小源码哈希锁入manifest，无环境安装或系统路径修改。白名单不含缓存和旧binary副本。只写本新目录；不改A041/037，不做MATLAB、扩招、stage/commit/push。
