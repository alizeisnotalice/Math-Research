# 第021轮：a1=151/100 单个五维有限并集的到限恢复

本轮固定旋转 `cos=4/5,sin=3/5`、`u^5=2`、`s=(1+u)/2`，只评估一个新参数 `a1=151/100`。完整登记的有限 witness functional 已计算完毕。实际向内见证比约为 **1.75056778573227106071475**；加入向外误差后，nominal 有限 functional 的区间端点小数显示为 **1.75056778573227106071475 至 1.75056778582070712136233**。精确有理端点在 `results.json`，小数只供阅读。

nominal 上端严格小于完整 `lower_round003/a2_equals_a3_results.json` 的 `C_cert` 下端，差距至少约 **0.24374867574777868324305**。本候选不能超过 C_cert。该结论只排除这个登记的有限 functional；没有证明全部尺度参数最优，也没有给出 C5 算子上界。L1 与 collar 资格由主代理另行判断。

## 输入、恢复边界和失败原因

先读仓库 AGENTS；本次只写 `lower_round021/`，017–020 均只读，不提交、不推送、不调用 MATLAB。计算使用本机 Python 标准库 `Fraction` 和整数，所有大小比较均为精确有理比较。

原程序 source SHA256：`62d307aa30983bbc38bb5747086123f50a74a744c16f75ddf0686f9e88839660`。
原登记 SHA256：`fe53db70a93d3bab80214433daac0a6c0275df3fa6b155a0e172911bd835cd86`。
原 PID 55866 已退出。它已保存完整几何、59/59 planar x-slabs（43011 cells）、完整三维直方图（26100 cells）。旧 `checkpoint.json` 的 44035 nodes 是 planar 最后一个 slab 的快照，并非退出时计数。

原全程预算共享 200000 nodes。库存1024 + planar43011 + unrotated26100 = 70135；最后 AND 需262×649=170038对，合计240173。原 check 先加一后判断，在200001时抛异常，已处理的129865对 AND 没有保存累加器，因此不可按旧 pair 索引续算。实际退出耗时约0.191秒，原因是 node 上限而非120秒。

`repair_registration.json` 在新计算前登记：只复用三个成功产物，从 AND prefix=0 恢复丢失的汇总；零几何或直方图重建。排序为 plane mask 数值升序，再 unrotated mask 数值升序；每个 plane-mask row 后原子保存精确 prefix、下一个二维索引、整数总量和 joint histogram。新预算最多20万 mask 对或120秒；不启动第二候选。原脚本、原登记、原checkpoint/limit以及三份成功输入的 SHA 全锁定且未修改。

## 五维机制与归一

raw 库存完整1024项，逐项保留判定：762项017统一跨度排除、198项旧已认证 planar 零面积、10项020精确零面积，剩54个 unrotated shape 归并成36个 `(K,n)` family。18个新 planar case 全部在 `new18_accounting` 中：10个有贡献候选family，8个由精确零面积排除；没有只套用旧26个family。

每个 family 的 `N=Kn`；`r_lower` 为80位二进制有理下界，严格满足 `(2r_lower)^5<N`。三未旋转轴的向内半径为 `floor(G*r_lower/a1)/G`、`floor(G*r_lower/s_upper)/G`（后者两轴）；旋转平面使用同一后者半径，`G=5*2^40`。所以五轴共同使用同一严格 `r_lower`，不存在不同轴各取最优半径。每个三轴盒保留 shape 身份，在盒交集以后才并 family mask。

225个 planar atom 的整数晶格索引与坐标均逐项去重。K 个未旋转 Cartesian 索引乘 n 个不同 planar 索引给出 Kn 个不同单位原子。向内盒的324条端点保护不等式、向外盒包含向内盒、共同半径条件均由精确有理核验。

旋转格点平移 `R(1,1)=(1/5,7/5)`、`R(3,4)=(0,5)` 保持完整 capture signature。stripe `[0,1/5]×[0,5]` 在平移及 y模5后，25份铺满 `[0,5]^2`。该 stripe 面积为1，`sum plane histogram=G²`；未旋转半周期立方体体积1/8，`sum unrotated histogram=G³/8`。

对每对 signature，只有 `plane_family_mask & unrotated_family_mask != 0` 才贡献，且该格子只计一次。plane 所有 case 在同一个 family mask 中联合，未将各 case 面积相加。整数汇总给出

```
V = 68279344437595937324167945193578864489431850775933089400625503
    /627710173538668076383578942320766641610235544446403451289600000
R_inner = 8*(151/100)*s^4*V
```

这里25倍 full-plane 质量与 periodmass=25抵消；仅对三个未旋转轴反射给8倍；物理 Jacobian 为 `a1*s^4`。没有在旋转平面重复乘反射因子。

## 向内见证与 nominal 误差

`inner_union_ratio_interval` 仅包住同一个实际向内有限并集的 `s` 代数误差，其两端约为

```
1.75056778573227106071475235229398219273528921583053047570598
1.75056778573227106071475504794594273908275236972527931869007
```

nominal 共同半径 `(Kn)^(1/5)/2` 的有限 functional 另有向外误差。设 family 的未旋转盒内外差体积和为 `du`，外盒并体积上界为 `vu`；225窗口 annulus 在全 `[0,5]^2` 的面积不超过 `900*(rho_hi²-rho_inner²)`，除25得到 stripe 误差 `dp<=min(1,36*(rho_hi²-rho_inner²))`。该 family 的乘积差不超过 `du+vu*dp`。对family求和仅用于误差上界，不作为并集体积。各项在 `shell_rows`。

保存的 `normalized_nominal_minus_inner_union_upper` 与根区间共同给出 `nominal_finite_family_ratio_interval`。对完整 C_cert 的比较为 `nominal_high < C_cert_low`，不是比较显示小数。验证器还用 C_cert 的 `u^5=2,v^5=3` 根隔离有理端点和25项系数重建其整个已存区间，逐项相等。

## 准确新旧比较

7个此前 partial planar case `(27,10),(36,10),(36,11),(36,12),(48,9),(48,10),(48,11)` 的全period marginal，均在**本轮实际 inward rho**处与020精确面积多项式完全相等；`(48,8)` 为25，`(8,1),(8,2)` 与017/018整平面证书一样为25。精确rho、面积和逐项判定保存于 `exact_marginal_checks`，未把020 nominal radius 当本轮 inward radius。

旧016完整旋转输入使用 `a_ref=1+u^3/2`、26 families，其区间约为1.73995838416010979至1.73995838418300350。本轮整个 nominal 区间严格高于旧016，差的精确区间端点小数显示约为 **0.01060940154926756 至 0.01060940166059733**。参数与family库存同时改变，故不能将该增量归因于单独加入7个case；本轮没有运行额外旧参数或消融候选。

## 实际运行、核验、复现

AND计算 PID **57687**：nice=10，单进程，Python active threads=1，六个数值线程环境变量均为1；170038新mask对，约0.1463秒，峰值25165824字节，完成prefix=170038、next indices=[262,0]。没有重算任何几何或直方图 cell。

独立保存证书核验 PID **57820**：nice=10、单线程，约0.0053秒，峰值23314432字节；核对输入SHA、1024库存逐项对应、324端点保护、225atom去重、完整prefix与joint整数总量、完整C_cert根及系数区间、两个新旧比较的精确正残差。该核验没有再次跑AND；它不是第二个几何并集引擎。原直方图的构造正确性由锁定原源码、产物及全部新增marginal交叉检查支撑。

已有完整结果时，日常审核只运行验证器：

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round021/verify_certificate.py
```

如在空的恢复输出目录中复现AND，应保留锁定输入及两份登记，使用同样环境运行 `continue_and.py`；如已有 `and_checkpoint.json`，它从精确prefix继续。完整prefix再次运行不会重复贡献，但本轮没有再次启动它。不要为复现最后一步而运行原 `finite_union.py` 重建成功直方图。

接续源码SHA256：`e822d72707522961fcaf2217eb4d5cbde0c62ad6e359f9d9a6f263b883f7333e`。
修复登记SHA256：`f3ecdd7deec0d69f065ffe647e5bacd86e0d7b1d93374c9b89d7abc194b95dd8`。
验证器SHA256：`192e0c763eef6d1409505b4c5ff04902115687b742704fcc7a18f8f04e59061f`。

必要归档白名单及逐文件SHA见 `archive_whitelist.json`。`start.flag`、Python缓存不归档。本轮未提交或推送，交主代理审核后整合。
