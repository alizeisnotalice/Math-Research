# 第231轮：局部密度帽与存活来源的联合几何预算

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮在原终端权w_g≥exp(-D)的实际分支上，将线性来源拥挤费降为二次对数；覆盖全部原实际L/R/s及三测试腿，不限定特殊输入模型。

## 新估计

ρ_V=||ν_V||/(qa^n)，A_D=(16/3)γ^n exp(D/2)，r_V=ρ_V exp(-U)/A_D，Ψp(r)=min(1,r)(1+log_+r)^p。原正交通切片满足

P_D,U ≤ 2Λ√(K/3)e^4 N_n A_D EΣ_V ||β_V|| min{204 Ψ3(r_V),264(D+1/2) Ψ2(r_V)}，N_n≤8+4√n log(n−1)。

关键是同一列同时使用真实接收域上的局部硬密度帽及QG列界，取较小值；再保原终端good的持久来源，早clock以配分函数log比支付≤D，末端短clock≤1/2。联合全部j与节点计数给常数11的二次对数价，先前带clock求和常数34给三次对数补充。没有重复领取硬帽或删除负信用。

D,U预先固定，原a,q与来源同V内不随j变化；其他对象须保未合并式。D,K固定、ρ≤(n+1)^B且B≥0固定的来源子族，其原far分支（U=0）实际费用≤C√n log³(n+1)Wβ。高deficit互补交通、一般加权来源对数矩、√K及其它补集仍未支付，不能视作全局闭合。

另证receiver区域分配与截断列公式：盒内用局部帽，盒外仅收真实未覆盖来源尾；硬窗口完全进入盒的来源退出尾账。矩形截断保四项乘积差和边界来源比例。尚未证明一般低价盒/责任权选择。

## 独审与实际复现

主估计、原持久资格、量词/常数、配套截断列及父稿均独审通过。210组节点配置（n至10^6、87空窗）以40/70dps复核，小n独立双循环对照；Fraction核验15组几何矩、6组三次身份、6组计数身份、38组比例转移。另6组原权形状的存活/删除模型经40/70dps分段求积复核。它们是有限计算证据，不是FULL或一般证明；二次阶由解析格点计数及存活来源恒等式证明。

独审还复跑配套稿内321项Fraction检查。最终父驱动在冻结包外实际运行两个主脚本，完整JSON均逐字节一致，12项manifest哈希运行前后通过。耗时约0.179秒（节点）与0.868秒（存活模型）；Python3.14、mpmath1.3.0。回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round231_20261017_final_replay/receipt.json`。无可归属本轮的实际额度扣减数据。

## 精简归档与下一步

[research.zip](research.zip)：13成员，80,129字节，SHA256 `6305898441578af6ddb42f81b694140c1e1ad9ced13122a1a8ad35b2b657b54a`。仅新增本README和ZIP两文件，另更新项目索引。保留详细LaTeX及必要复现材料，不编译PDF。

检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round231_20261017_final_checkpoint`；ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round231_20261017_research.zip`。逐成员核验work、检查点、ZIP及备份后才移除本轮work。

解压后 `python3 reproduce.py --verify-only`；完整包外复现：

```sh
PYTHONPATH=/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor python3 reproduce.py --output-dir <新建包外目录>
```

下一最小步骤：在原双后验积分里联合终端deficit与来源尾，以真实持久势跨层结算盒内费用与侵蚀外尾费，避免每层重复领整份Wβ。未证明这一平均总价前，不将二次对数矩改名当作已付。
