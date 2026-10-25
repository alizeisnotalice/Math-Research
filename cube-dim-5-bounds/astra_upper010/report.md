# Astra upper 010：一个二阶形状的必要角点排除

本批唯一预登记形状没有通过必要角点筛选，未启动全空间生产，也没有新g或C5界。旧ℓ1下全部尺度a的角点天花板约0.00004641780940800445，严格低于现有全域g=2821/50000000。即使改用真实d^(5/11)，加权Jensen给出的天花板仍严格低于128×10^-6；因此仅提高determinant下界无法让此固定形状跨过该目标。

## 形状与全域正性

令R=M-I，S=Σ_{i<j}r_ij²，T(x)=x'R²x，W=Σr_ij x_i x_j，选择

f=1+W/10+(3/44)S-T(x)/40+W²/80，Z=af。

PSD、diag M=1给0≼M≼5I，因此M²≼5M。于是R²=M²-2M+I≼3M+I，并有T≤4|x|²+6W。全cube上

f≥1-|x|²/10-W/20+W²/80≥1/2-W/20+W²/80
=9/20+(W-2)²/80≥9/20。

所以Z≥0且Young不等式合法。这个形状保留x'R²x和Σr²，确实超出仅依赖W的旧族；本次排除只针对所选系数，不能外推到所有二阶形状。

在角点x=(1,…,1)，T=2Σr_e²+2Σ相邻无序边对r_er_f，W²=Σr_e²+2Σ全部无序边对r_er_f，因此

880f=880+88Σr_e+27Σr_e²-22Σ相邻边对r_er_f+22Σ不交边对r_er_f。

该式有66项，S5不变。选用003锁定ℓ1_L=A_Ld-B_Ud²，计算U=E[ℓ1_L f^7]、V=E[b f^11]，按唯一规则取a=671/1000。实际有理角点下端约0.00004641754145882875，低于009的0.00005642。

更强的全部a必要上限为c U^(11/4)/V^(7/4)。其四次幂以c_U有理上端计算，准确比较后严格小于(2821/50000000)^4，故不能靠重新调a挽救这个旧ℓ1候选。

此外准确计算B=E[f^7]，以f^7加权的凹Jensen给真实分数幂预算

c(E[df^7])^(5/4)B^(3/2)/V^(7/4)
≈0.00010868288291841895<0.000128。

此排除同样由有理四次幂比较给出，不以小数作证明。它没有低于10^-4，所以没有排除该固定形状在更强determinant下界下达到10^-4。以上均为Young证书族的必要角点天花板，不是G上界。

## 有限准确计算

最高图次数max(10+14,4+22)=26。先复制005既有143289项图缓存，再对每个所需canonical图键检查存在；缺项则停止，不允许先算新矩。实际所有键都存在，`new_moments_computed=0`。原缓存与本批副本运行后SHA256完全相同：

`6215f0d5db14cd5150927bbfea8272e50af1d7b1531cc26bfcb7b7fc5e0c48f6`。

覆盖原理也成立：新f单项的任意r图在相同次数W^n展开中都会出现，符号平均只保留适当顶点奇偶性；将空间指数置换为降序后，005会保留对应代表。实际逐键核查加强了这一覆盖论证。

为避免完整展开数千万级单项，保存S5轨道总系数。若P不变且某轨道单项系数为a，总系数a|O|；乘完整不变f时，各旧轨道代表贡献的canonical分布相同，所以直接用总系数乘f的66项并聚合是准确的。前3次幂分别用不压缩的完整多项式展开并canonical聚合，全部准确相等；raw单项数66、1001、8008，对应轨道数5、29、140。

f^7最终19008轨道，f^11最终565985轨道。乘d、d²、b时保留这些因子的完整展开，不能把两个因子同时独立压轨道。期望只需Eulerian图；实际分别评估108186、2995122、1482、9995631个乘积。所有整数系数、图矩、累加使用GMP。每边编码指数≤24<32、总次数≤26，显式逐边核查，没有固定宽度矩溢出。

每个完整幂都在/tmp原子写checkpoint，保存源码SHA及逐行FNV校验；f7单独保留。恢复从完整幂继续，U可重新读取f7重建，不重复叠加部分和。程序复核缓存图矩源指纹和校验，并显式禁止任何缺失图矩求值。

Apple clang21 C++17编译以`-Wall -Wextra -Werror`通过。角点生产实际15.6244秒，Python Fraction/256-bit Arb核验约0.001秒。附加准确核验：E[f]=643/616；0<E[d²f7]≤E[df7]≤E[f7]；E[df7]²≤E[d²f7]E[f7]。所有有理矩、尺度、角点分数和天花板四次幂保存在两个JSON中。

## 复现与归档

kernel源SHA256：`e2e2eeff385f6022db6457824f8d8f40deacfb66818c10d0e73d27431d338953`。后端为锁round012，SHA256 `7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8`。

```sh
clang++ -O2 -std=c++17 -Wall -Wextra -Werror -isystem /opt/homebrew/Cellar/gmp/6.3.0/include -I /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper cube-dim-5-bounds/astra_upper010/corner.cpp -L /opt/homebrew/Cellar/gmp/6.3.0/lib -lgmpxx -lgmp -o /tmp/astra_upper010_corner
cp /tmp/astra_upper005_graph_cache.tsv /tmp/astra_upper010_graph_cache.tsv
/tmp/astra_upper010_corner cube-dim-5-bounds/astra_upper010/corner_moments.json /tmp/astra_upper010_graph_cache.tsv /tmp/astra_upper010_f_power.tsv e2e2eeff385f6022db6457824f8d8f40deacfb66818c10d0e73d27431d338953
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-5-bounds/astra_upper010/verify.py
```

从头重建005缓存需先运行005锁定kernel；本批不扩空缓存。仅核本批结果可直接运行Python验证器，消费保存的四个准确混合矩。新的轨道重权重没有第二种完整高阶独立矩实现，主要依赖已审后端、前三幂独立展开及实际零扩缓存检查。

必要白名单为登记、C++、Python、两个JSON、报告、manifest。大幂表、图缓存、二进制不归档。下一机制不应以此固定f继续盲目提高determinant次数；其它形状与条件机制仍未排除。
