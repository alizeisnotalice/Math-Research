# 第023轮：精确C格的第一轴period3质量一侧变分

对 `a1=1+u^3/2`、`a2=...=a5=(1+u)/2`、`u^5=2` 的原C格，登记并计算第一轴质量

```
w_j(theta)=1+theta*b(j1 mod3), b=(2,-1,-1).
```

两个方向分别取 `theta=+epsilon`、`theta=-epsilon`，`epsilon→0+`。全部 `ell∈{1,2,3,4}^5`、3相位和必要平移像从raw进入一阶并集算法。两方向精确导数相等且严格负，其有理隔离端点的小数显示为

```
-0.534343679305780438274408514328629046238946817331476097129286
-0.534343679305780438274408514328629046238946817331476097126958
```

故本有限functional满足 `R(theta)=C+D*|theta|+O(theta²)`，D<0。这里只断言存在未量化邻域内的局部下降；不宣称R在更大范围为偶函数，不给出具体epsilon区间，也不排除一般weighted inputs、未登记长块或其他尺度。没有正方向，因此未登记或运行任何有限epsilon候选，更没有C5改进宣告。

## 查重、登记和原C身份

本clone和主source `/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/` 均先查period3、b模式、weighted/parity和first-variation登记。未发现本period3族已被实际排除的记录。

已有失败与本族准确区分：006为第一轴period2交替质量、theta=1/8、1/4、1/2三个有限点，且格距来自有理dyadic近似；主source alt014为精确C格上的二维parity `(j1+j2)`、`(j2+j3)`，右导数分别约−0.335058844179和−0.280525163635，并有theta=1/1024有限点排除。这些均不涵盖第一轴period3的 `(2,-1,-1)`。alt014两个导数证书SHA分别为

```
1056f416440f5ffb1fd632a1b5353f4b9316add6a013c2bbc0fc9a7e9a3f351d
41ebf7aabaa261f90be88e648766e3cc3d11de650b37eb9c1ebc5e5eddf41de6
```

新登记保存于 `registration.json`，先于实际数学计算；source与依赖及既有失败证书逐项哈希锁定。只写新023目录，021/022及其他旧文件全部只读；未提交、推送或调用MATLAB。计算使用本机Python整数、Fraction和现有Q(u,v)代数域例程，无程序wall上限；每个方向进程仅一个有限对象、最多200000 union节点，每first-axis slab保存准确prefix。

theta=0时原子全为质量1，格距正是完整003 C源。3a1周期包含3个第一轴格胞，而首轴由半格胞恢复全格胞给2倍；四tail轴反射仍给16倍。若原C半格胞并体积为U，则 `C=32U`，本域体积 `V0=6U`，因此 `16V0/3=C`。本轮的theta=0值由此身份与原完整证书继承，没有把基线重算登记为新实验。jet常数项与原C的全部25个基系数相等仅用于核对上述身份。

完整C source SHA256：`f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2`。

## 完整有限family、质量和starts资格

先枚举1024个长度向量与3个第一轴起点相位，得到每方向3072项 `mass_phase_catalog`，没有从C点42个活盒起步。首轴block起点j、长度ell1时

```
N=Π_k ell_k,
delta_m=sign*sum_(t=0..ell1-1) b((j+t) mod3)*Π_(k=2..5)ell_k,
m(epsilon)=N+epsilon*delta_m,
r(epsilon)=m(epsilon)^(1/5)/2,
r0=N^(1/5)/2, r1=r0*delta_m/(5N).
```

每个块的五轴、十个面共同使用同一r(epsilon)，而非独立半径。N由2和3的幂相乘，所以r0、r1都在原Q(u,v)中，`v^5=3`。程序逐项精确检查 `32*r0^5=N`。

对两方向在 `0≤epsilon≤1/4`，每原子质量至少1/2、至多3/2。period3内质量总和始终3，质量导数严格为0，输入原子仍位于互异的直积格点。block原子是各坐标连续整数索引的Cartesian积，基数严格为N；质量斜率直接求和，与上式一致。

取接收域 `[0,3a1]×[0,s/2]^4`。对全部长度最多4的块，N≤1024、m≤1536，故r<5/2；又a1>3/2，因而2a1>3>5/2。首轴起点j≤−2时其区间上端 `j*a1+r<0`；j≥5时下端 `(j+ell1-1)*a1-r>3a1`。所以所有可能与接收域相交的登记块恰被安全范围 **j=−1,...,4** 覆盖。每相位取两个平移像，共 **6144 raw候选/方向**，不是只取每相位一个不完整盒。

质量仅依赖首轴；固定首轴block及tail长度后，tail的所有平移block质量相同。对tail半胞元 `[0,s/2]`，最近的奇/偶连续索引块逐点最小化最大捕获距离，因此canonical tail区间支配其他starts。这个约化在加权后依然合法。四个tail反射保持输入，故归一化准确为 **16V/3**。不除格胞体积，也不使用旧32反射。

## 一侧union、tie与零宽机制

每个端点用jet `(x0,x1)` 表示 `x(epsilon)=x0+epsilon*x1+O(epsilon²)`。一侧epsilon>0次序按x0、再x1作Q(u,v)精确lex比较；乘法为 `(x0*y0,x0*y1+x1*y0)`。接收端点先作与固定域边界的max/min，再检查非空、去重、支配和完整五维union。

所有raw候选都经历该检查。若某轴有严格负常数宽度，该盒在小邻域仍为空；若常数宽0但一阶宽正，保留其一阶窄条；若常数与一阶宽均0，该轴宽至多O(epsilon²)，其总体积对一阶无贡献。两个端点同位置异速度必须分开排序；同位置同速度可以合并，因为只会差O(epsilon²)。

严格理由是：有限个正质量块的根函数在0附近解析，全部端点的Taylor余项统一为O(epsilon²)。有界轴平行盒对端点位移的体积误差由各面有界截面积乘位移控制；union对称差不超过各盒对称差的和。因此以一阶端点代替真端点、合并同值同速端点，以及按jet支配去盒，都只改变O(epsilon²)体积。这也处理了零宽、并脸、冗余和高阶才开口的盒，而不是假设原42盒一直足够。

实际每方向：6144 raw中5965有严格负常数宽度，179为一阶可行raw，179个不同jet盒，去支配后147盒。本次没有新增零常数宽度的可行盒；这由完整raw判定得出，不是预先排除零宽机制。kept端点中同位置不同速度的组数为 `(6,13,13,13,13)`，因此确有必须一侧分离的并脸tie。

递归slicing与独立坐标覆盖mask算法分别计算union的常数、导数及五个坐标的同步面贡献，三者逐项完全相等。面贡献是同一个共同半径变分的乘积法则分解，不是五个独立优化方向。

## 精确导数证书

令 `15D=Σ_(i,j=0..4) A_ij*u^i*v^j`，两方向均得到以下整数矩阵（行i、列j）：

| i\j | 0 | 1 | 2 | 3 | 4 |
|---:|---:|---:|---:|---:|---:|
| 0 | −136 | −288 | −432 | −80 | 192 |
| 1 | −416 | 468 | 432 | 176 | 0 |
| 2 | −98 | 324 | 0 | 24 | −288 |
| 3 | 116 | 36 | −288 | −232 | −96 |
| 4 | 143 | −180 | 144 | 256 | 72 |

上述D的Q(u,v)精确系数在 `plus_derivative.json` 和 `minus_derivative.json` 中逐项相等，并非由显示小数判断。两方向的五个面贡献也相同，展示约为 `(-0.0238563691599291,-0.127621827536463,-0.127621827536463,-0.127621827536463,-0.127621827536463)`，各项准确系数与隔离区间都已保存，且其和在域内严格等于D。

用原分母10^60的u、v有理根隔离端点，逐项对D向外包围，上端严格负，符号已分离；不存在仅依赖浮点负号的步骤。这里的 `D_minus=d/d epsilon R(-epsilon)|0+<0`，故若改用theta坐标，左导数为 `R'_theta(0−)=−D_minus>0`。

即使±b不由平移互换，也可由实际两个一侧证书推出这个一维有限functional在0为严格局部峰。没有外推一般质量扰动或全局theta最优性，没有把某个模型代价下降当Cn结论。

## 实际运行、检查点和复现

两进程顺序运行，均nice10、单线程、六个数值线程环境变量为1，未设wall hard kill。

| 方向 | PID | 秒 | 峰值字节 | sweep/mask节点 | first-axis slabs |
|---|---:|---:|---:|---:|---:|
| plus | 69609 | 35.4945 | 135217152 | 1672/13939 | 85/85 |
| minus | 69727 | 26.8261 | 135102464 | 1691/14261 | 85/85 |

没有实际fault或到限，没有盲重启。每完成一个first-axis slab原子保存prefix及准确部分体积jet；两个检查点都已完整。若未来需恢复到限任务，应先核对方向、source hash和first-axis cut hash，仅从已保存prefix后继续；本轮没有发生恢复计算。

保存结果核验PID **69955**：约0.1064秒、峰值29425664字节，零新union任务。实际核验两份3072相位质量catalog、两份147个kept盒的2940个端点jet、共22812个显式原子索引/质量求和检查、完整prefix、C常数身份、两方向同一精确导数、同步面之和及严格负隔离区间。021、022全部白名单文件与锁定旧source/失败证书SHA均未变。

首次产生证书使用以下两条命令；结果已完成时只需第三条审核，不应重复union：

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round023/period3_variation.py plus
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round023/period3_variation.py minus
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round023/verify_saved.py
```

source SHA256：`9b91991ad377f1e6f1577b8863fe1917b9bce453a6933f07c07444764d8b8e8e`。
登记SHA256：`b816c73fa53b466412fe2123e51b97b88620c38e0b8fad3722f7de38a80f281f`。
plus导数SHA256：`c3446a276498a9421f87a5c5b77188f770b44a41380b5a4930fdd2e9ff482235`。
minus导数SHA256：`9433a26dd3a63a6bfb2975a3982708a58da7981adf5e22ac779ae7cc69209b6c`。

未解决：未量化局部下降邻域；未计算有限非零epsilon；未登记长度大于4的块或其他质量方向；L1/collar由主代理保留终审。本轮必要文件及SHA白名单见 `archive_whitelist.json`，未提交或推送。
