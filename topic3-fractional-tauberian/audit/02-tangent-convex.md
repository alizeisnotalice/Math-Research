# 切向多面体与凸体部分核对记录

核对基线为 PDFKit 恢复的原稿第 13--16、18--22 页：
`original-p13.txt` 至 `original-p16.txt`，以及 `original-p18.txt` 至
`original-p22.txt`。本记录不把旧修订稿当作数学依据。

## 第 6 节：中心方体切向多面体

|原编号与页码|原命题及假设|本稿处理与依据|状态与新标签|
|---|---|---|---|
|原 (26)，p.13|$Q_0=[-1/2,1/2]^n$，$c=\min\{\beta,n-\beta\}$，$\sigma=2c$，定义 $N_c$；$k=\lfloor\sigma\rfloor$，并包含 $k=0,n$ 端点。|按同一参数与排序约定重排为正常分段公式；恢复 PDFKit 在求和处的排版断裂。|已核对；`tan:normal-form`、`tan:polytope`。|
|引理 6.1、原 (27)，p.13|对所有 $v$，$N_c(v)=\inf_{u\ge0}\{cu+\sum_j(|v_j|-u/2)_+\}$。|按 $z=u/2$ 分段斜率证明；显式处理 $\sigma\in\mathbb N$ 与 $\sigma=n$。|严格证明；`tan:parameter-elimination`、`tan:infimum`。|
|定理 6.2、原 (28)--(34)，pp.13--15|固定 $Q_0$、$t\downarrow0$；外扩与内缩中心见证方体给出全部一阶可行中心，并排除遗漏中心。|保留两条对数展开、紧参数集上一致 $O(t^2)$、见证方体的 $\ell=1+O(t)$ 与中心 $x=O(t)$、缩放光晕统一有界支集，以及边界零测后的支配收敛。原稿定理的结论是固定 $Q_0$ 的极限，不是全体 Tauberian 常数等式。|严格证明；`tan:unit-cube-limit`、`tan:outer-expansion`、`tan:inner-expansion`、`tan:uniform-support`。|
|定理 6.2 体积计算，p.15|$|S_{n,c}|=2^n/[k!(2c)^{n-k}]$。|重算有序正锥的 Beta 积分，再乘以 $2^nn!$；$k=n$ 单列为 $\ell^1$ 球。|严格证明；`tan:cube-limit`。|
|推论 6.3、原 (35)，p.16|由固定 $Q_0$ 测试，给全体常数的 $\liminf$ 下界。|保持量词：常数的上确界中选 $E=Q_0$，故只得 $\liminf\ge|S_{n,c}|$。没有改写为极限、上界或 sharp 等式。|严格推论；`tan:global-liminf`、`tan:global-liminf-formula`。|
|注 6.4，p.16|当 $c>1/2$，旧候选 $c^{-n}$ 被新下界严格超过；真正问题是全体集合是否以 $Q_0$ 为一阶极值。|保留反证方向；该原开放首项问题由第9轮gl:main闭合，近极值分类仍待研究。未采纳旧审计中“原稿把固定模型跳为全体 sharp 等式”的错误描述。|准确状态说明；`tan:status`。|

## 第 8 节：固定中心对称凸体

|原编号与页码|原命题及假设|本稿处理与依据|状态与新标签|
|---|---|---|---|
|引理 8.1、pp.18--19|原稿假设 $K=-K$、$0\in\operatorname{int}K$；两个位似体的对称差至多为 $\varepsilon|B_0|$。|完整保留尺度比较、同尺度平移归约、内接欧氏球的切片下界与紧环带排除。证明实际只用欧氏半径和凸性，故将引理推广到任意有界满维凸体并以欧氏位移表述。|严格证明；`cvx:near-homothetic`、`cvx:near-homothetic-bound`。|
|定理 8.2，pp.19--20|固定 $K$ 的非中心阶为 $1-\alpha$，中心阶为 $(1-\alpha)^n$；给出 $B_\alpha,L_\alpha$ 的显式下界。|重写上界的见证比较链，并将中心/非中心的量词分开；下界仍取 $E=K$，用内含和外包两个机制。长度尺度严格为 $B_\alpha=\alpha^{1/\beta}$、$L_\alpha=\alpha^{-1/(n-\beta)}$，不是归一化指数 $\theta=\beta/n$ 的幂。|严格到既有双密度稳定性为止；`cvx:endpoint-order`、`cvx:endpoint-order-formula`、`cvx:endpoint-lower`。|

## 第 9 节：卷积体、投影体与同似一阶变分

|原编号与页码|原命题及假设|本稿处理与依据|状态与新标签|
|---|---|---|---|
|命题 9.1、原 (40)--(42)，p.20|一般凸体 $K$，锚定同似基 $x+sK$；$E=K$。|逐行代入定义，保留完整同似表示与固定尺度 $s=1$ 的卷积体嵌入。|严格证明；`cvx:convolution-embedding`、`cvx:embedding`、`cvx:homothetic-representation`。|
|引理 9.2、原 (43)--(44)，pp.20--22|任意凸体 $K$；协方差函数沿方向的切片一阶导数。|保留 Fubini 截线公式、支配收敛、统一线性夹逼、径向函数一致收敛和 Hausdorff 收敛。缩放明确是 $(1-\alpha)^{-1}C^K_\alpha\to|K|\Pi^*K$，不是错误的 $(1-\alpha)^{-1/n}$。|严格证明；`cvx:covariogram-limit`、`cvx:hausdorff-limit`、`cvx:covariogram-volume`。|
|推论 9.3、原 (45)，p.22|由固定尺度卷积体嵌入，给 $E=K$ 的中心光晕下界。|保持为固定 $K$ 的 $\liminf$；未提升为同似光晕等式或全体集合常数。|严格推论；`cvx:projection-liminf`、`cvx:projection-liminf-formula`。|
|问题 9.4，p.22|一般凸体同时允许尺度和平移的一阶切向体；原稿把它列为开放问题。|用径向函数和边界缺陷泛函完成任意有界满维锚定凸体的固定模型切向极限；切向体可不对称。|固定模型已闭合；全体集合极值问题仍开放，见下行。|

## 新的固定凸体模型结论

对任意有界满维凸体 $K$ 且 $0\in\operatorname{int}K$，定义
\[
D_K(v,a)=\int_{\partial K}[-\langle v+ay,\nu_K(y)\rangle]_+\,d\mathcal H^{n-1}(y),
\quad
N_{\beta,K}(v)=\inf_{a\in\mathbb R}
\left((n-\beta)a+\frac{D_K(v,a)}{|K|}\right).
\]
通过径向函数的精确体积积分、Minkowski 泛函在几乎处处可微方向上的隐函数求导、
径向面积公式及支配收敛，结合第 8 节的近重合引理，可得
\[
t^{-n}|H^K_{\beta,\mathrm c}(1-t,K)|\longrightarrow
|\{v:N_{\beta,K}(v)<1\}|.
\]
证明写入 `cvx:first-variation` 与 `cvx:full-tangent-limit`。其范围严格限于固定测试集合
$E=K$；切向体一般不必中心对称。不涉及全体可测集的 Tauberian 上确界，因此不声称任意凸体、
椭球或方体是全局极值体。

## 保存的本机复现入口

- `verification/baseline.py` 实际检查双密度凸性、矩形首项及171个切向体积参数案例。
- `verification/polytope_tangent.wl` 以Wolfram15精确算术枚举正方形和三角形的对偶顶点，检查支持函数和反射恒等式，并计算圆盘方向积分。
- `verification/polytope_tangent.m` 以MATLAB R2026a的linprog、convhull及数值积分独立复核同一组几何例子；完整结果摘要保存在results/。

## 保留的边界

1. `cvx:full-tangent-limit` 不推出
   $C_{\beta,\mathcal B^K_{\mathrm c}}(\alpha)$ 的 sharp 首项；仍缺任意 $E$ 的统一匹配上界。
2. `cvx:full-tangent-limit` 已推广至非中心对称固定凸体，但这不改变其“固定 $E=K$”范围；
   卷积体推论仍只给固定尺度的严格下界。
3. 本轮未使用 Petty 投影不等式宣布极值体，也未把原稿关于二维普通卷积体的文献性讨论
   扩展到含尺度参数的分数次问题。

## 固定模型的几何对偶与可计算例子

新增 `cvx:reflection-duality` 由散度恒等式给出 N_{n-β,K}(-v)=N_{β,K}(v)。新增 `cvx:polytope-duality` 给出S为加权切片立方体线性像的极体，原始/对偶LP均可行且最优值有限。正方形、非对称三角形的精确顶点与面积以及圆盘积分实际由Wolfram15和MATLAB R2026a交叉核验。证明与数据分别见正文和proofs/polytope-duality.md。

## 第5轮：中心对称半阶的仿射极值

新增half:projection-identity、half:ellipsoid-extremum，证明中心对称K、β=n/2下N=hΠK/|K|。既有Petty投影不等式以定理A单独编号，所用规范与等号条件已回到原作者PDF核验。原9.4的这一固定模型子类得到椭球极值，任意β及固定基全体输入的极值问题保留。脚本half_order.wl与half_order_matlab.m实际核对绝对值配对恒等式、7个维数的Γ常数、4000个方体方向及非对称反例的代价。

## 第9轮状态更新

原稿固定模型定理与liminf推论的量词保持原样；新gl:main以共同框架、质量修改、面层容量上界和凸纤维对称平均补齐任意可测E的统一上界，确定中心方体全体首项。此推广不包括任意凸体基，也不确定全体次首项系数。
