# Astra upper 007：四次行列式下界的全域证书

## 结论

在固定 α=3/2、16个未归一化体对角方向的001模型中，严格认证

**G(x)≥g=5563/100000000，所有 x∈[-1,1]^5。**

这将005辅助证书4583/10^8提高至5563/10^8，约提高21.38%。本轮没有认证新的 C5 上界；仍需与足够紧的 I 预算配合。

唯一预登记候选是006的四次 determinant minorant ℓ3_L，η=1/8，a=679/1000，Z=a(1+ηW)。没有额外搜索或增添参数。所有2401有理 Bernstein 系数均已实际生成并检查，最小值出现在(6,6,6,6)，不是根据对称性或角点经验推断。

## 解析链

001给 G(x)=c E_R H(M,x)，其中 c=3/(2^(15/2)π^(5/2))，d=det M，b=xᵀadj(M)x，H=d^(5/4)b^(-7/4)。006已证明全区间0≤d≤1上

0≤ℓ3_L(d)≤d^(5/11)，

并保存共同向下因子 T_L 和四个准确有理系数。此处逐字读取006的有理数，没有从显示小数重建。004的多仿射顶点论证给整个 cube 上 -2≤W≤10，因此 Z≥0。

Young不等式逐点给

H≥(11/4)ℓ3_L(d)Z^7-(7/4)bZ^11。

分布 R 在坐标符号变换和置换下不变。期望右侧是 Xi=xi² 的对称有理多项式 P(X)。我们准确重建该多项式；在正面 X1=1、其余四个 Xi∈[0,1] 上，Bernstein基非负且和为1，所以 P≥最小Bernstein系数=L。c_L=947374313793/10^15<c 且 L>0，于是该面上 G≥c_L L。

G的符号及置换对称性覆盖全部面。由齐次性 G(rx)=r^(-7/2)G(x)，r=||x||∞≤1，得到整个 cube 的同一下界；x=0时G发散为正无穷。此步骤使用G的齐次性，没有误把非齐次证书多项式P当作齐次函数。

准确最小系数为

`645800317046901128982409516324162568551416330275340785910867733668414283063954522673601 / 10997385607996319663598268467014529932224564955168374784000000000000000000000000000000000`。

准确 g=c_L L 是

`611814632209609797266827199895316545126774367932855227801533641645299064016153977450861607091278593 / 10997385607996319663598268467014529932224564955168374784000000000000000000000000000000000000000000000000`，

约0.00005563273436231541。向下取到10^-8得到5563/10^8。准确数与006角点值相等，是本次全部系数验证后的结论。

## 计算规模与正确性

计算 E[d^k W(x)^n]，k=1,2,3,4、n≤7，以及 E[b(x)W(x)^n]、n≤11。determinant 本身不含空间变量，前者总 X 次数 n；后者总次数 n+1。W的每一项对每个 x_i至多一次，b对每个 x_i至多二次，所以符号平均后每个 Xi次数≤6。实际输出127个对称partition系数，代入X1=1及恢复轨道后为1316个非零面单项，tensor多次数(6,6,6,6)，共2401系数。

partition表的每条值是**一个给定单项的系数**，不是整轨道的总系数。恢复空间多项式时枚举不同坐标置换一次。006角点生产曾压缩determinant轨道；本轮没有使用那项优化，完整展开d⁴的171705项，保留d与W(x)的相关性。只对结果的空间单项按对称partition保存。

GMP后端固定为 round012，源SHA256 `7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8`。复用带源指纹的006图矩缓存，空间累加从空开始。最高图次数27；矩整数、所有determinant系数、阶乘和累加均为任意精度整数。50位编码每边占5位，行列式四次幂每边指数≤8，与W相加的保守界≤19<32，因此整数键相加不会串位；最高图次数通过参数与断言核验。

每个完整W次数完成后同步图缓存并原子替换raw结果。若中断，重新从W^0累计空间数据，使用已完成准确图矩缓存；没有把部分空间和重复加回。验证器拒绝未完成W^11的raw输入。

Apple clang21、C++17、`-Wall -Wextra -Werror`实际编译通过。空间生产7.11516秒，图缓存保持147258 canonical项，无新增图矩。这符合006已覆盖同一图集合的预期。Python3.15/Fraction全证书检查约0.075秒。

实际交叉核验：全部004的d、d²、b空间系数准确一致；新d³、d⁴空间多项式恢复角点后与006全部16个新矩准确一致。新d³/d⁴仍共享同一图矩后端，后一项不能伪称独立数学实现。其空间累加路径与006角点轨道累加不同。

Bernstein全表规范序列化SHA256：`04c8cf9a2bc5c5ebeec3e26a1edc6ec3a7ecf06cae00f8443b3fad58d54ee5dd`。保存完整127项P及所有必要空间矩即可准确重新生成2401项；没有重复保存可重建的大表。

## 复现和独审接口

必要文件：registration.md、graph_moments.cpp、verify.py、space_moments.json、verification.json、report.md、manifest.json。raw整数累加和图矩缓存位于/tmp，不归档。后端与004/006是显式锁定依赖。

直接审核有理证书（无需GMP生产）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-5-bounds/astra_upper007/verify.py cube-dim-5-bounds/astra_upper007/space_moments.json
```

从图矩重建：

```sh
clang++ -O2 -std=c++17 -Wall -Wextra -Werror -isystem /opt/homebrew/Cellar/gmp/6.3.0/include -I /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper cube-dim-5-bounds/astra_upper007/graph_moments.cpp -L /opt/homebrew/Cellar/gmp/6.3.0/lib -lgmpxx -lgmp -o /tmp/astra_upper007_graph_moments
/tmp/astra_upper007_graph_moments /tmp/astra_upper007_graph_moments.json /tmp/astra_upper007_graph_cache.tsv
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-5-bounds/astra_upper007/verify.py /tmp/astra_upper007_graph_moments.json
```

空缓存可从头复现；现有缓存仅加速。验证会刷新runtime字段，锁定输出应备份后再运行。独审优先核查006 Taylor共同因子、空间partition语义、127项到2401项的独立Bernstein转换，以及最高次数27的新矩。所有结论以准确有理数为准，不用约5.563e-5的小数反推。

本轮没有宣称 K3 族已到天花板。更高奇阶determinant minorant、非线性Z或条件机制均未被此证书排除；是否继续需另行登记有限批次。
