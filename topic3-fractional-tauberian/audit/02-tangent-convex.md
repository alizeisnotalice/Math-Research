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
|注 6.4，p.16|当 $c>1/2$，旧候选 $c^{-n}$ 被新下界严格超过；真正问题是全体集合是否以 $Q_0$ 为一阶极值。|保留反证方向和开放问题；未采纳旧审计中“原稿把固定模型跳为全体 sharp 等式”的错误描述。|准确状态说明；`tan:status`。|

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
|问题 9.4，p.22|一般凸体允许尺度与平移的切向体及仿射极值问题。|本轮保留原稿的固定尺度下界和待证问题。|开放；后续逐轮研究。|
