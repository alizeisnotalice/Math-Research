# A086：以独立四子集支付 rank≥4 pattern 余项

本轮提交一个完整的 aggregate 引理和精确有理证书接口，替换 076 中需要逐个
rank≥4 pattern 付常数的步骤。主审保留最终证明验收。四个新小符号配置一次运行
全部通过；没有实际 body32 常数、full6D 密度、body kernel、全域 g 或新 Cn 结果。

## 数学结果

给每个 rank≥4 的固定参数 pattern L 分配其字典序首个独立四列 J⊂L，再为每 J
固定两个实际列的 completion N_J。所有物理微分都在这个固定六列 basis 下进行。
逐 Leibniz 项取绝对值后，才能把 assigned patterns 扩大到全部 L⊃J 的非负和。
未解变量的 χ/η 和为 1；两个 completion 列保留绝对导数和
S_r=|(cχ)^(r)|+|(cη)^(r)|，没有把它错误地替换成 |c^(r)|。

若 R 为固定 cutoff 半径，0<t≤1，则

\[
S_r\le C_rt^{-1-r},\quad
\pi(C_0,C_1,C_2)=(1,1+4/R,2+8/R+12/R^2).
\]

四个选中尾的 0、1、2 阶导数仍各为 O(t)。因此每项在 k 次微分、有 q 次命中
选中尾时保留 t^(2−k+q)，从而整个 sum 有明示的 k=0,1,2 上界

\[
|D_{u_1}\cdots D_{u_k}p_t^{\ge4}|
\le t^{2-k}\sum_{J\text{ indep.4}}
\frac{h_k\prod_r\|N_J^{-1}u_r\|_1}{\pi^6R^8|\det N_J|}.
\]

h_0=1、h_1=max(4/R,1+4/R)，h_2=max(20/R²,2+8/R+12/R²,h_1²)。
逐 pattern 指数规模的常数和已被独立四子集和替换；证明允许非负过计数，但原
pattern 分配仍唯一。没有对 rank 事件、indicator 或 basis selector 微分。

进一步保留时间幂并在 (0,T] 积分，可用同一个固定 PSD 物理矩阵 M_4 同时控制
所有方向的 Hessian。它是 Σ_J q_J N_J^-T diag(d_tail I4,d_comp I2)N_J^-1，
其中

\[
\begin{aligned}
d_{tail}&=(a_2+3a_1^2)T^3/3+a_1b_1T^2,\\
d_{comp}&=(b_2+b_1^2)T+2a_1b_1T^2,\\
a_1&=4/R,\quad a_2=20/R^2,\quad
b_1=1+4/R,\quad b_2=2+8/R+12/R^2.
\end{aligned}
\]

有 −M_4≤D²R_{≥4,T}≤M_4。Young 余量有显式正平方和，因此没有将逐方向
min 拼成矩阵。有理 A、R、T 给出精确有理 π⁶M_4；tr M_4 是该余项的一个
各向同性支付额。它不是完整 Green 的 Γ。face 参数化 E 必须使用 EᵀM_4E。

理论对全空间 x 成立，不借助局部 S2 距离条件；应用 076 时仍需与它的同模型
cutoff 参数一致。重复、平行及零列均不会破坏固定基或正和。R_{≥4,T} 实际为
C²，零端的 value/gradient/Hessian 分别由 t²、t、1 可积支配。

## 接口及有限计算证据

aggregate.py 使用 Fraction 精确线性代数，提供 independent_four_bases、
scalar_constants、majorant_from_bases 和 certificate。最高层函数自行枚举全部
独立四子集及固定 completion；不枚举实际模型的 2^m patterns，不调用密度/kernel。
低层 majorant_from_bases 明确要求输入 records 已完整覆盖所有 J。脚本的直接 CLI
拒绝自动执行；本轮唯一调用者是预登记小 checker。

四配置全部 PASS：

| 配置 | 实际工作与范围 |
|---|---|
| C1 | 一个新的非轴有理 6×8 toy frame，含一对平行列和一条零列。检查其全部 256 toy patterns，108 个 rank≥4 patterns 唯一分配到 25 个独立 J；固定 completion 的 determinant 为 ±12 或 ±24。没有密度计算。 |
| C2 | 六个已解因子总微分阶数≤2 的 28 个命中多重指标，各对两个 completion 与两个未解变量的 16 种选择精确因式分解；四条尾未丢失。 |
| C3 | 一个固定非轴 N（det=12）和两条已登记物理方向，六个抽象二次 jet 的一阶及混合二阶微分与 N^-1 Leibniz 式严格相等。 |
| C4 | R=3/5、T=2/7 下的 28 个时间指数，15 项正平方和，固定物理矩阵拉回恒等式，以及同一 toy frame 的 aggregate 精确有理矩阵。 |

toy 的 π⁶tr M_4 为 51983388671875/4374822312，完整有理矩阵及 K_0,K_1,K_2
见 results.json。该数字仅属于登记的 6×8 toy，不是 body32 Γ，更不是 C6 界。
这些检查支持有限代数环节；全称证明由 theorem.md 的固定积分、正和、支配微分和
时间积分完成，不能用有限 fixture 的 PASS 替代。

## 执行和封存总账

登记上限四个新配置（用户上限五个），一次 check_fixed.py 调用，tool chunk
cbada7，exit0。实际四配置、四通过、零失败/重试。checker 实测 wall
0.748749792 秒、CPU 0.615081 秒。程序使用默认 170 秒总 alarm，清理/写结果时
不取消；独立硬 CPU 上限 25 秒，满足任务 CPU 30 秒预算。没有子进程或等待回收。

body/native kernel=0、full6D density=0、实际 body32 常数实例化=0、旧 076 三配置
重跑=0、MATLAB=0、递归 agent=0。没有 stage、commit、push。输入 076 theorem 和
manifest 仅只读绑定；本轮没有改其他代理的文件。

本轮写代码后曾按主审优先级暂停给 A082 runtime 独审，WIP.json 记录当时实际
配置为零。A082 已单独封存；恢复本轮后才写定 theorem/pre_run seal 并执行四配置。
历史 WIP 文件保留本地但不作为最终状态或归档工件。两轮结果没有混账。

## 仍未支付的部分

本轮只关闭 076 rank≥4 余项的 aggregate 支付接口。实际 in-plane 完整 fan 与
振幅导数、quotient coercivity/weighted/signed 层、其他 remainder、上尾、真正
scalar vertex 数据和全域覆盖仍需另行给出。没有计算实际 frame 上 Σ_J，也没有
宣称所得常数在实际六维问题中足够小。最终 mathematical acceptance 由 main 决定。
