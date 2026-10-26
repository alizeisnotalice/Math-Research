# 第022轮：静态四维 section 与物理第一轴仿射接口

对固定旋转 `cos=4/5,sin=3/5`、固定 `s=(1+u)/2,u^5=2` 和登记的54个shape/36个family，构造了一个适用于整个 `a∈[3/2,2]` 的静态四维联合section。第一轴保持物理坐标，完全去掉 `floor(G*r/a)`。逐section独立取最大值的安全放松，加上uniform nominal误差后，给出精确有理上界，其小数显示为

**1.91366364740209885946918330720016630568827303528609641274158**。

该上界严格低于完整003 `C_cert` 的下端，精确正差的小数显示约为 **0.0806528141663869451361920916892**。因此这个固定有限functional在整个登记窗口没有超过C_cert的机会。这里没有选取、搜索或宣称最优a，也没有把结论扩大成所有尺度、所有family或一般C5算子上界。L1与collar资格由主代理终审。

## 对象和来源

本轮只写 `lower_round022/`，旧021和所有旧轮次未改。复用021锁定的262个plane family signatures及其36个family的共同严格 `r_lower`。未读取或使用021的649个三轴unrotated masks生成新结果。两个tail轴从54个shape各自的021 tail区间拆出，shape身份一直保留到与plane family的联合之后。

021的1024项raw库存已全部会计：762项017在整个 `[3/2,2]` 窗口的跨度排除、198项旧认证planar零面积、10项020精确planar零面积、54项剩余shape。新18个planar case均继承会计：10个有贡献family、8个精确零面积排除。这里依赖这些锁定旧证书的适用范围；本轮没有重新进行旧证书的完整几何证明审计。

新输入和来源SHA在计算前保存于 `registration.json`。主程序源码SHA256为 `e782fe9e4e3ce1bd8ab27f2001f1fb9278ef782017e9093ef3d370d3cd95d73c`；登记SHA256为 `767001f6a7e4ab5f40279baeda2aaf3d5d405c4317126aa15a61f8e7be252d90`。所有运算使用本机Python整数与Fraction。每进程180秒hard wall、nice10、单线程、最多20万新4D mask对；检查点记录精确prefix，可恢复。

## 精确section与仿射引理

设shape的第一轴block长度为ell，标准连续晶格索引集合为I_ell。在物理半周期 `[0,a/2]` 内，该shape的共同半径r所允许的第一轴区间为

```
J_(ell,r)(a) = [max(0,max(I_ell)*a-r), min(a/2,min(I_ell)*a+r)]
```

若右端不大于左端则为空。这四类端点都是**有理仿射函数**，因为r取021已经认证的严格有理 `r_lower`。每个family有 `(2*r_lower)^5<Kn`；tail与plane的向内阈值仍不超过同一 `r_lower/s_upper`。没有为不同坐标选择不同半径。

直接化简可得 `h_(ell,r)(a)=max(0,min(a/2,r-(ell-1)*a/2))`。奇数ell的J为左锚区间 `[0,h]`，偶数ell的J为右锚区间 `[a/2-h,a/2]`。因此同一ell中r较大的区间包含r较小者。对某个固定四维section S，只需保留每个ell出现的最大r，形成descriptor `D_S`，且

```
L_S(a) = min(a/2, max_odd_ell h_(ell,r)(a) + max_even_ell h_(ell,r)(a))
```

空类的max取0。压缩前后严格相等，保留shape身份的目的正是防止tail盒的不同shape互相拼接。

对descriptor的原始第一轴端点，再加0和a/2，列出所有两两仿射交点。任意两个相邻交点之间，原始端点的顺序、裁切/空集状态、区间覆盖和并集端点身份都固定，故L_S为有理仿射函数；它在交点连续。因此每段的两个端值准确确定斜率和截距。每section的最大值必在该有限分段端点达到。本轮生成158个descriptor、合计364段，完整系数、定义域和section权重在 `results.json:descriptors`。这不是017的1946个事件区间扫描，也不是a点网格。

## 仅一次四维联合权重构造

两个tail坐标仍使用021的向内网格 `G=5*2^40`，坐标域 `[0,1/2]^2`。54个tail盒给出900个网格cell和421个**shape masks**；总整数质量 `G²/4`。plane的262个family masks逐个展开为其全部shape bits，然后与tail shape mask求AND。共 **262×421=110302** 个mask对，得到4088个四维shape signatures（含0），总质量 `G^4/4`。

每个section权重为其整数质量除以 `G^4`，因此权重总和1/4。一个section包含某shape，当且仅当两个tail坐标属于该shape的同一个tail盒，且plane捕获其family所需的n个不同atom。所有shape都在同一个section共同mask中，第一轴取这些区间的并集，重叠只计一次。没有把各family面积求和当并集。

给定权重W_S，冻结四轴、第一轴精确物理区间的实际向内functional为

```
R_inner(a) = 8*s^4*Σ_S W_S*L_S(a).
```

第一轴物理长度已含a的尺度，故Jacobain只余 `s^4`，不得再乘a。旋转格点的平移 `(1/5,7/5)` 和 `(0,5)` 保持plane capture signature，25个stripe铺满 `[0,5]^2`；25倍plane面积与periodmass25抵消。8只来自三个未旋转轴的反射，不对旋转平面再乘反射。

225个plane atom索引和坐标去重、Kn个不同Cartesian原子及共同半径，沿用锁定021证据；021所有白名单文件的SHA本轮全部再次核对未变。新joint histogram还逐shape核验：其54个一阶marginal必须等于**该shape tail盒面积 × 该family plane marginal面积**。这同时检查了shape与family对应关系。

## uniform上包与nominal误差

对每个section独立取 `M_S=max_(a∈[3/2,2]) L_S(a)`，则对窗口内每一个共同a，

```
R_inner(a) <= 8*s_upper^4*Σ_S W_S*M_S.
```

不同section可以在不同a达到最大值，所以右边是安全放松，并不是某个可实现的见证比。冻结functional的这个上界约为 **1.91366364717699216220156942035**。

nominal半径r_N=(Kn)^(1/5)/2位于021存储的 `[r_lower,r_upper]` 中。第一轴r增大dr使裁切区间增量长度至多2dr；窗口上端A=2给出任意第一轴区间长度至多A/2。对某family各shape，令tail向内/向外盒面积为b_i/b_o，则其三轴接收集差的uniform上界是

```
du_f = Σ_shapes [2*(r_upper-r_lower)*b_o + (A/2)*(b_o-b_i)].
vu_f = min(A/8, Σ_shapes (A/2)*b_o).
```

225个plane窗口的annulus bound，除以25个stripe后，为

```
dp_f = min(1,36*(rho_upper²-rho_inner²)).
```

故名义family乘积与冻结乘积的差至多 `du_f+vu_f*dp_f`。对family求和**只用于误差上界**。其总误差在乘8s_upper^4前约为 `2.11211035809294399734504449066e-11`，各项均在 `shell_rows`。得到整个窗口统一的nominal有限functional上界1.913663647402…。

与C_cert使用完整有理区间作 `C_low-window_upper>0` 精确比较；没有用显示小数决定符号。完整C_cert根隔离与25项系数在锁定021验证记录中已重建相等，本轮重新验证其SHA。

## 运行、验收与恢复

唯一计算PID **63173**：nice10、单线程、六个线程环境变量全为1；约0.311秒、峰值24313856字节。110302新mask对全部完成；checkpoint prefix=110302，next=[262,0]。没有中途失败，没有重复生成成功直方图，没有复跑021的151/100并集，没有第二个4D histogram。

核验PID **63280**：约0.1056秒、峰值24739840字节；零新joint mask对。核验输入锁、021全部白名单SHA、权重总量、完整prefix、54个tail与joint marginals、158个descriptor精确权重、364段与独立通用区间并算法在端点及三个内部有理点的一致性（1820次），并精确核对放松上界、误差和C正差。分段内部全体成立还依赖上面的仿射端点顺序引理；这些有限点核验本身不替代该引理。

完整结果和checkpoint已保存，日常审核只需运行验证器：

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round022/verify_sections.py
```

`sections.py` 先复用已有tail histogram，再从精确4D prefix继续；完成的4D histogram可直接用于后续读取仿射接口。本轮未再次启动它。hard wall时保留最近的原子检查点与成功产物，每批资源上限不是研究总时限。

本轮只交付静态section/PWA接口及足够低的安全上包。未合并全局事件求最优a，未生成新模型或新参数网格，未提交、推送或调用MATLAB。必要文件和逐项SHA见 `archive_whitelist.json`。
