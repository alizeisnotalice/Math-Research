# R045：全部二轴角的共同源价格与正常槽剩余预算

2026-11-05。固定 \(n\ge2\)。新局部证明候选，供独立审核；不在本稿宣布全输入八阶上界。所有实际源、见证、CDF 盒及槽分配均取自同一个预先固定有限末层。仓库及旧冻结稿只读。

本稿补 R042 正常 packet 之外的全部二临界边。核心是同一源角先相加所有有界目标的实际质量，再取一个正部禁层；从同一份完整 \(J\) 中只扣一次这层和一次余量。超长二轴角只需避开固定的**正常小管**，不要求避开最后覆盖全部有界层的大管，故阈值次序无循环。

## 1. 实际输入和共同正常集合

记 \(a=v^n\)、\(k=n^2/27\)、\(S=\sum_i r_i^n\)。沿 R040 的原有限严格见证族，采用 R042 §1 的实际
\[
 F_i\subset G_i,\quad R_i=(|G_i|/a)^{1/n},\quad
 X_i=x_i+W_i+[-R_i/2,R_i/2]^n,\quad
 \mathscr D_i=(x_i+W_i^0)\setminus F_i,
\]
\[
 J_i=(R_i+|W_i|^{1/n})^n-(R_i+|F_i|^{1/n})^n.
\]
所有源及角对求和均取外部实际源指标 \(i\notin I\)；\(S\) 仍可包含原全选族。对 \(x\in F_i\)，原指定严格见证 \(Q(x,r(x))\) 包含完整 \(G_i\)，且 \(r(x)>R_i\)。不同 \(F_j\) 两两不交，\(F_j\cap G_i=\varnothing\) 对 \(j\ne i\) 成立。相交度有固定上界，相交尺度可比；原／有效／中点源盒各宽介于 \(cvR_i,CvR_i\)。所有常数允许依赖原固定选择常数。

保留同层定量粗先验
\[
 \sum_i\varepsilon_iR_i^n=O(v^5S),\quad
 |\mathscr D_i|\le C\varepsilon_i aR_i^n,\quad
 \sum_{\rm deep}(\alpha_i+\Delta_i/R_i)R_i^n=O(v^6S),
\tag{1}
\]
\[
 \sum_i\xi_iR_i^n=O(v^{5/2}S),\quad0\le\xi_i\le C,\quad
 S(\partial I)=O(v^5S).
\tag{2}
\]
\(\Delta_i\) 始终是端移，不是删除质量。另用 R042 (1.3)–(1.4) 的强 CDF：对任意 \(Y_i\subset(x_i+W_i')\setminus F_i\)，完整 \(J_i\) 同时供应 \(Y_i\) 的鲜导数成本、固定正系数的
\[
 \frac va|\mathscr D_i\setminus Y_i|
 \quad\hbox{和}\quad R_i^{n-1}\Delta_i.
\tag{3}
\]
这里是一次实际体积余量分配，不能对不同边集分别重新使用全部 \(J_i\)。

对于深处二临界边 \(e=\{i,j\}\)，写 \(L_1,L_2>0\)、非临界投影交体积 \(b_e\)、\(R_e=\max(R_i,R_j)\)，并采用 R042 (1.5)：
\[
 p=\frac{L_1}{v^3R_e},\quad q=\frac{L_2}{v^3R_e},\quad
 z=\frac{\operatorname{clip}((R_i-R_j)/2;\,-\max L_\ell,\max L_\ell)}
           {v^3R_e},
\]
\[
 x=\frac{1+v}{n}p,\quad y=\frac{1+v}{n}q,\quad
 \zeta=\frac{1+v}{n}z,\quad
 \widehat b=\frac{b_e}{(1+v)^{n-2}R_e^{n-2}},\quad b=\min(1,\widehat b).
\tag{4}
\]
固定足够小的 \(\nu>0\)，定义本稿与 R043 **共用的正常集合** \(\mathcal N\)：
\[
 |x-1/3|<\nu,\quad |y-1/3|<\nu,\quad|\zeta|<\nu,\quad 1-b<\nu.
\tag{5}
\]
先取 \(\nu\) 小于 R042 正常紧窗及中心分离常数，最后取 \(v\) 小。于是根差未截断，\(R_i/R_j=1+O(v^3)\)，且
\[
 c_n v^3R_e\le L_1,L_2\le T_0v^3R_e,\qquad
 b_e\ge(1-\nu)R_e^{n-2}.
\tag{6}
\]
故本正常集合严格满足 R043 (2) 的所有几何要求。其中心是 \(n/[3(1+v)]\)，没有另换成 \(n/3\) 定义第二个集合。同源同轴对同符号至多有一个正常目标，仍由实际中心分离得出。

记
\[
 \mathfrak E=\sum_iR_i^{n-1}\Delta_i,\qquad
 A_6(v)=\frac{k}{(1+v)^2}v^6,\qquad \Gamma=\frac{10}{9}k.
\]
免费面量记为 \(\mathcal F\)，使用原已证的两个接口
\[
 \sum_{\rm deep}\alpha_iR_i^n\le C\mathcal F,\qquad
 \text{一次真实免费面条恢复的正体积}\ge c_0\mathcal F.
\tag{7}
\]
这些条与全部 \(X_i\) 和被标记部分分离，且重数有界。后者不是从 \(J_i\) 再造的源成本。

## 2. 与所有固定窗无关的非正常强制性

令
\[
 f(s,t)=s_+t_+(s_++t_+),\quad
 D(x,y,\zeta)=\tfrac12[f(x+\zeta,y+\zeta)+f(x-\zeta,y-\zeta)],
\]
\[
 G(b,x,y,\zeta)=b/27-bxy+D.
\tag{8}
\]
对 \(x,y\ge0\)、\(0\le b\le1\)、\(|\zeta|\le\max(x,y)\)，有 \(D\ge xy(x+y)\) 和 \(G\ge0\)。在 (5) 外，存在只依赖 \(\nu\) 的 \(c_\nu>0\)，使
\[
 G\ge c_\nu(b+D).
\tag{9}
\]
重要的是 \(c_\nu\) 不依赖后选的大窗 \(T\) 或小角下界 \(m\)。

证明：凸性给 \(D\ge xy(x+y)\)，故 \(xy\le(D/2)^{2/3}\)，并有
\[
 G\ge b/27-\sup_{x,y\ge0}xy(b-x-y)
       =b(1-b^2)/27\ge0.
\]
若 (9) 不成立，取 \(b+D>0\) 且 \(G/(b+D)\to0\) 的反例序列。\(b=0\) 时该比值为 1，可先排除；若 \(D/b\to\infty\)，用
\(bxy/D\le2^{-2/3}(b/D)^{1/3}b^{2/3}\to0\)，亏损比趋于 1，不可能。若 \(D/b\) 有界而 \(b\to0\)，则 \(D\to0\)、\(xy\to0\)，亏损比仍有统一正下界。余下 \(b\) 有正下界、\(D\) 有上界；近取等迫使 \(b\to1\)，再迫使 \(xy\) 有正下界。因此 \(x+y\le D/(xy)\) 有界，\(\zeta\) 也有界。取极限后，唯一取等为 \(b=1,x=y=1/3,\zeta=0\)：前两变量由三数 AM--GM，最后由正宽处的严格凸性。与 (5) 外矛盾。\(b=D=0\) 无需控制，直接满足 (9)。

正常窗内保留 R042 已证的显式二次下界
\[
 G\ge c\{(1-b)+(x-1/3)^2+(y-1/3)^2+\zeta^2\}.
\tag{10}
\]
不将非正常边的完整 \(D\) 压成横权 \(bD\)；即使 \(b\) 很小，真实整横截面的源成本仍在 (9) 中。

## 3. 阈值次序与所有二轴边的分组

先固定 \(\nu\) 和原 CDF 余量常数。正常目标的标称禁层全部位于一个固定小管 \(\mathcal C_i(K_N)\)，其中 \(K_N\) 只依赖 \(n,\nu\) 和原选择常数，不依赖下述 \(T\)。令 \(c_*>0\) 为 (3) 的一个固定可用余量下界。

取固定大 \(T>T_0\)，使第 7 节的长角源成本转换常数 \(C_L/T\) 小于 \(c_*/100\) 及第 5 节由 (9) 固定保留的完整成本系数 \(c_\nu'/100\)。\(c_\nu'\) 与 \(T,m\) 无关，故可在此使用；\(C_L\) 也只依赖正常小窗和原几何常数，与 \(T\) 无关，其来源在第 7 节逐项说明。随后取
\[
 0<m\le\min\{n/12,\ k/(16T)\},\qquad
 C_{\rm soft}\ge\max\{1/m,16T/k,12/n,1\},
\quad \omega_e=\min(1,C_{\rm soft}\min(p_e,q_e)).
\tag{11}
\]
所有参数先于 \(v\) 固定，\(\omega_e=1\) 在正常边上。

全部深处二轴边分成：

* 正常 \(\mathcal N\)；
* 有界非正常 \(\mathcal B\)：\(m\le p,q\le T\)，且不在 \(\mathcal N\)；
* 有界小角 \(\mathcal D\)：\(\max(p,q)\le T\)，\(\min(p,q)<m\)；
* 长角 \(\mathcal L\)：\(\max(p,q)>T\)。

边界二轴边另作 \(\mathcal E_\partial\)。对 \(\mathcal D,\mathcal L,\mathcal E_\partial\) 不在有界共同禁层内追加目标。它们的全部交损仍在后面结算；不删去任何原盒基线。尤其不按 \(G_v\) 筛掉二轴边：非 \(G_v\) 端点仍参加全部共同层、退款及余量结算；\(G_v\) 只在第 9 节的实际槽分配中使用，其漏片按角权收费，绝不按 gross 交量剔除。

## 4. 同一见证、同一正部、同一余量

对每源按轴对及面对符号分组，在同一两个排序和同一源横向核心内，先相加 \(\mathcal N\cup\mathcal B\) 中全部目标的实际帽质量。所有目标真实 \(F_j\) 不交，且与完整自身 \(G_i\) 排他。R037 (11) 的同一个实际轴帽控制给
\[
 s\ge \sum_{j\in\mathcal N\cup\mathcal B}A_{ij}(y)-CR_iE_i(T),\qquad
 A_{ij}(y)=\frac{T_j^\perp P_{ij}(Q_{ij}-y)_+
                  -|\mathscr D_j\cap\mathcal P_j(y)|}
                    {naR_i^{n-1}}\ge0,
\tag{12}
\]
\[
 \sum_iE_i(T)R_i^n=O_T(v^6S).
\]
\(A_{ij}\ge0\) 是帽中实际 \(F_j\) 质量除以容量导数，不能把该非负性忽略。所有源角共用足够大的核心剪裁 \(H_Tv^3R_i\)，使组间层不交。

在 (12) 正部以下取真实禁层 \(Y_{i,\mathrm{all}}\)。**以完全相同的 \(E_i(T)\)、核心及排序**，只留下正常目标，得 \(Y_{i,N}\)。于是
\[
 Y_{i,N}\subset Y_{i,\mathrm{all}}\subset
      (x_i+W_i')\setminus F_i,\qquad
 Y_{i,N}\subset\mathcal C_i(K_N).
\tag{13}
\]
同角正常与非正常层无需互不相交；这里只对一个相加高度取一次正部。

定义完整原始横截面名义非正常成本
\[
 B_{\mathcal B}=
 (1+v)^{n-1}\frac va
 \sum_{\rm sources,\ groups,\ sorts}
 \int_{\rm source\ core}\int
       \sum_{j\in\mathcal B}
       \frac{T_j^\perp P_{ij}(Q_{ij}-y)_+}{naR_i^{n-1}}
                 \,dy\,d x_\perp.
\tag{14}
\]
它含实际目标整横宽 \(T_j^\perp\) 和实际源核心横宽；不乘接触小权 \(b_e\)。这就是共同成本中非正常目标的线性名义份额。

设
\[
 R_{\rm all}=\frac va\sum_i|\mathscr D_i\setminus Y_{i,\mathrm{all}}|,
 \qquad
 R_N=\frac va\sum_i|\mathscr D_i\setminus Y_{i,N}|.
\]
由 \(u\mapsto u_+\) 的单调性及 1-Lipschitz 性，
\[
 \boxed{R_N\le R_{\rm all}+(1+v)^{-(n-1)}B_{\mathcal B}
              \le R_{\rm all}+B_{\mathcal B}.}
\tag{15}
\]
此式没有 \(T\) 相关乘法常数或未记的小量；其原因是 (14) 保留了完整原始横截面成本。若事先将它压成 \(bD\)，(15) 就不能用于这里的吸收。

当 \(T\) 固定且 \(v\) 足够小，有界全部层位于某个固定大管
\[
 Y_{i,\mathrm{all}}\subset\mathcal C_i(K_{\rm all}),\qquad
 Z_{\rm all}=\frac va\sum_i|[(x_i+W_i')\setminus F_i]
                                \setminus\mathcal C_i(K_{\rm all})|
 \le R_{\rm all}.
\tag{16}
\]
\(K_{\rm all}\) 在 \(T,m\) 后选，可任意增大。后面不要求长角源条避开这个大管。

## 5. 有界联合退款与固定的非正常折扣

将 (12) 乘鲜导数并积分，实际 \(J\) 同时给
\[
 \sum_iJ_i\ge B_N+B_{\mathcal B}
       -\mathfrak F_N-\mathfrak F_{\mathcal B}
       +c_*(R_{\rm all}+\mathfrak E)-C_Tv^8S.
\tag{17}
\]
\(B_N\) 是与 (14) 相同的正常名义成本；\(\mathfrak F_N,\mathfrak F_{\mathcal B}\) 是各目标实际缺量的 Tonelli 核上包。高度侵蚀为 \(O_T(v^2\sum E_iR_i^n)=O_T(v^8S)\)。这里没有抽走固定比例的正常成本。

对正常退款，可重放 R042 (3.4)–(4.4)，因为它允许同一较大的 \(E_i(T)\) 和核心，仅使绝对 \(C_Tv^8S\) 常数改变。用 (15) 得，对任意固定小 \(\epsilon>0\)，
\[
 \mathfrak F_N\le \Gamma v^7W_N+\epsilon\mathfrak G_N
       +Cv(R_{\rm all}+B_{\mathcal B}+\mathfrak E)+C_{\epsilon,T}v^8S,
\tag{18}
\]
\[
 W_N=\sum_{\mathcal N}b_eR_e^2,\qquad
 \mathfrak G_N=v^6\sum_{\mathcal N}G(b,x,y,\zeta)R_e^n.
\]
式 (18) 中乘余量及 \(B_{\mathcal B}\) 的 \(C\) 可取与 \(T\) 无关：正常帽核只用固定正常窗的 \(Cv^2/a\)，(15) 常数是一；\(T\) 仅进入已经显式允许的绝对八阶误差。

对有界非正常退款，核有界于 \(C_Tv^2/a\)。目标缺集先分为 \(Y_{j,N}\) 与其余。正常层每胞体积至多 \(Ca v^5R_j^n\)，而 \(\mathcal B\) 中 \(p,q\ge m\) 给每条边的 (14) 成本至少 \(c_m v^6R_e^n\)。所以
\[
 \mathfrak F_{\mathcal B}\le C_{T,m}v
                  (R_{\rm all}+B_{\mathcal B})+C_Tv^8S.
\tag{19}
\]
这个系数可依赖固定 \(T,m\)，因为它只须在最终 \(v\) 小时吸收，不用于选择 \(T\)。

正常领先实际成本与交量之差仍精确给
\[
 A_6(v)W_N-c\,\mathfrak G_N+O_T(v^8S).
\tag{20}
\]
完整高维因子是 R042 修复版 (5.2) 的 \(n^2(1+v)^{n-4}v^6R_e^n\)，未用不等尺度的 \(R_j/R_i\) 作无量纲省略。

对 \(\mathcal B\)，实际完整 (14) 与理想完整 \(D\) 比较时，根比误差是 \(Cv\) 倍完整成本，核心和目标宽度误差的总和是 \(O_T(v^{17/2}S+v^8S)\)。单条原宽不需一致接近一；用 (2) 的加权平均支付。又原宽上下界给 (14) 与 \(v^6\sum_{\mathcal B}D R_e^n\) 之间固定双侧比较常数，均与 \(T,m\) 无关（只在取小 \(v\) 的阈值上依赖 \(T\)）。结合 (9)，可以固定 \(\gamma_\nu,c_\nu'>0\)，均不依赖 \(T,m\)，使
\[
 \sum_{\mathcal B}I_e-B_{\mathcal B}
 \le (1-2\gamma_\nu)A_6(v)W_{\mathcal B}
      -c_\nu'B_{\mathcal B}+C_Tv^8S.
\tag{21}
\]
这里 \(W_{\mathcal B}=\sum_{\mathcal B}b_eR_e^2\)，软权等于一。由 \(b\) 截至 \(\widehat b\) 的误差，使用
\(\sum(\widehat b-1)_+R_e^n=O(v^{7/2}S)\)，也是更高阶绝对误差。保留 (21) 中固定的完整 \(B_{\mathcal B}\) 份额，正是为了吸收 (18)、(19) 及后面长角对坏层的支出。不能将整份非正常亏损只换成一个角权折扣后再次花费。

## 6. 有界小角与外部边界

对 \(\mathcal D\) 不计任何源成本或退款，直接有
\[
 pq\le T\min(p,q)\le(k/16)\omega_e
 \le\tfrac14\,\frac{k}{(1+v)^2}\omega_e
\tag{22}
\]
对足够小 \(v\le1\) 成立。第二步在软权未饱和时用 \(C_{\rm soft}\ge16T/k\)，饱和时用 \(\min(p,q)<m\le k/(16T)\)。因此小角的 gross 本身有固定折扣。

边界二轴角的原单源界是
\[
 \sum_{\mathcal E_\partial}I_e
 \le C v^2\sum_{\rm neighbors\ of\ boundary}\varepsilon_iR_i^n
          +Cv^4S(\partial I).
\]
由 \(\varepsilon_i aR_i^n\le C|\mathscr D_i|\)，将缺集分成 \(Y_{i,N}\) 和其余，前者每胞不超过 \(Ca v^5R_i^n\)，后者用 (15)。邻边界总权 \(O(v^5S)\)，有限相交度给
\[
 \sum_{\mathcal E_\partial}I_e
 \le Cv(R_{\rm all}+B_{\mathcal B})+Cv^9S.
\tag{23}
\]
这包括边界端自身；若其没有正常层，则其缺集全在余量中。没有调用边界 normal profile 的定性小量。

## 7. 所有长而退化的二轴角

取 \(R_i\ge R_j\) 的所有者，令物理无量纲长短交长
\[
 P=\max L_\ell/r_i,\quad Q=\min L_\ell/r_i,\quad
 H=\min(Q,PQ/v^2),\quad 0<Q\le P\le C_Av.
\]
\(\mathcal L\) 上 \(P>(1-O(v))Tv^3\)。

### 7.1 正常层在任意长帽中的占比

设目标原始帽 \(\mathcal P\) 在两临界轴的边长是固定倍数 \(Pr_i,Qr_i\)（已截至真实目标宽度，保留固定比例），其余轴是目标原始整宽。则
\[
 |Y_{j,N}\cap\mathcal P|
 \le C\left(\frac{v^3}{P}+v^2\right)|\mathcal P|.
\tag{24}
\]
常数只依赖正常小窗及原宽界，不依赖 \(T,Q\) 或最后的大管。

证明：一个正常层至少有两个轴的坐标支撑长度不超过 \(C_Nv^3r_j\)。若它的轴对含长帽轴，该轴给比例 \(C_Nv^3/P\)，不论短轴多小。若不含长轴，至多一个层轴能是短帽轴，因此至少一个其他轴在帽中取的是完整宽度 \(\ge cvr_j\)，层的该轴只占 \(C_Nv^3r_j\)，给 \(Cv^2\)。对固定数正常角求和。\(n=2\) 时每个层轴对都含长轴，第二种情况为空。有效面位置或 \(E_i(T)\) 的侵蚀只移动或缩小支撑，不改变这个长度估计。

故取 \(T\) 充分大并令 \(v\) 小，正常层最多占每个保证目标帽的四分之一。若
\[
 |\,\mathscr D_j\setminus Y_{j,N}\,|
        \ge c\,a r_i^n PQ/v^2,
\]
该边全部交损可记为 \(Cv(v/a)|\mathscr D_j\setminus Y_{j,N}|\)。有限度和 (15) 给这些目标余量例外的总支出
\[
 Cv(R_{\rm all}+B_{\mathcal B}).
\tag{25}
\]
这一步不用整缺集的 \(\varepsilon_j\) 直接判坏，因其正常两臂不能被当成可任意再付一次的钱。

### 7.2 面缝／端移例外只花 \(Cv\) 的真实正预算

其余边若 \(\alpha_i+\Delta_i/r_i\ge cH\)，用 \(PQ\le CvH\) 得
\[
 I_e\le Cv(\alpha_i+\Delta_i/r_i)r_i^n.
\]
合计不超过
\[
 Cv(\mathcal F+\mathfrak E).
\tag{26}
\]
端移由 (3) 的余量吸收；免费面只能由 (7) 那次真实条恢复支付。不得从 \(\sum J_i\) 另外扣出一份“免费面成本”。

### 7.3 正常例外之外的实际源条与两管次序

对剩余边，重复 R020/R044 的同一个实际见证容量矛盾：长轴有效面深取 \([c_2P,2c_2P]r_i\)，短轴取 \([c_3H,2c_3H]r_i\)，其他轴取距两面均 \(>K_Nv^3r_i\) 的核心。常数 \(c_2,c_3\) 只依赖原几何，先固定。若条中存在 \(F_i\) 点，它的同一个见证包含完整 \(G_i\)，并保证目标原始帽各临界交长至少为固定倍数 \(Pr_i,Qr_i\)，截于真实宽度后仍如此。由 (24)、目标余量非例外和非临界整宽捕获，至少留下固定正倍数 \(a r_i^nPQ/v^2\) 的实际 \(F_j\) 质量。

真正浅面的原点态容量另一方面给剩余容量至多
\(C a r_i^n(c_3+c_{\rm exception})H\)。取固定小系数使它严格小于上述质量，因为 \(H\le PQ/v^2\)。故整条都位于 \(D'_i=(x_i+W_i')\setminus F_i\)。它的源成本满足
\[
 \frac va|\mathcal R_e|
 \ge c r_i^n PQ\min(1/v,P/v^3)\ge cT I_e
\tag{27}
\]
对充分小 \(v\) 成立，常数 \(c>0\) 与 \(T,m,K_{\rm all}\) 无关。

取 \(T>4K_N/c_2\)。长轴和其余核心轴都远离正常小管，仅短轴可能靠近一面，故
\[
 \mathcal R_e\cap Y_{i,N}=\varnothing.
\]
所有源条即使重叠，最多有固定相交度重数。因此 (27)、(15) 给
\[
 \sum_{\mathcal L\ {\rm remaining}}I_e
 \le\frac{C_L}{T}R_N
 \le\frac{C_L}{T}(R_{\rm all}+B_{\mathcal B}).
\tag{28}
\]
\(C_L\) 与 \(T\) 无关：(24) 的正常支撑常数是固定的，穿透及容量比例仅依赖原宽界，源核心只用 \(K_N\)，条重数是原固定相交度。这说明第 3 节选 \(T\) 的条件没有暗含 \(T\) 自己。

这一步**不要求** \(\mathcal R_e\) 在 \(\mathcal C_i(K_{\rm all})\) 外。它与有界非正常层重叠的那部分已由 (15) 的完整 \(B_{\mathcal B}\) 记账，其 \(C_L/T\) 份额由 (21) 吸收。最后才选择覆盖全部 \(Y_{\rm all}\) 的 \(K_{\rm all}\)，所以不存在 \(T>K_{\rm all}(T)\) 的循环。

## 8. 全二轴共同价格

令 \(W_O\) 为所有不在 \(\mathcal N\) 的二轴边（含边界）的 \(\sum b_eR_e^2\omega_e\)。有界非正常的软权为一，其余软预算非负。合并 (17)–(28)，先取固定 \(\epsilon\) 小以吸收正常 Young 项；再依第 3 节取 \(T\)，用 (21) 的固定完整成本份额吸收 \(C_LB_{\mathcal B}/T\)；最终令 \(v\) 小，以吸收 \(C_{T,m}v\) 的所有余量支出。可以保留固定的 \(c_R R_{\rm all}\)、正常 \(\mathfrak G_N\)、端移及一份完整非正常成本。

特别存在固定 \(\gamma_\nu,c_R>0\)，使
\[
 \boxed{\begin{aligned}
 H_2+c_R Z_{\rm all}
 &\le [A_6(v)+\Gamma v^7]W_N
        +(1-\gamma_\nu)A_6(v)W_O\\
 &\qquad +Cv\mathcal F+C_{n,\nu,T,m}v^8S ,
 \qquad H_2=\sum_{q(e)=2}I_e-\sum_iJ_i .
 \end{aligned}}
\tag{29}
\]
\(\gamma_\nu\) 可预先取小于 (21) 剩余折扣和小角折扣；长角及边界交量已由余量支付，故给它们任何非负的 \((1-\gamma_\nu)A_6\omega b\) 上包都合法。式中每个 \(J_i\) 只在 (17) 使用一次。

对任意预定 \(0<\eta<c_R\)，可只保留 \(\eta Z_{\rm all}\)。这供 R044 从同一份预算中分配；不再把 \(Z_{\rm all}\) 或 \(J\) 全额各给另一个模块。

## 9. 实际槽空缺的领先正预算

采用 R036 同一 \(C_{\rm soft}\) 的全二轴 adaptive-clean 四份分配。正常边软权为一。记
\[
 \mathcal C=\frac14\sum_sR_i^2|F_s|,\qquad
 D_N=\frac14\sum_sR_i^2|F_s\setminus O_s|\ge0,
\tag{30}
\]
其中 \(O_s\) 是**同一正常集合** \(\mathcal N\) 的正常 clean 占用并。R036/R043 的已证漏片收费、根权传输给
\[
 W_N+W_O\le\mathcal C+O(v^3S),\qquad
 W_N=\mathcal C-D_N+O(v^3S).
\tag{31}
\]
对实际重数先取同槽占用并；不假设一条正常边独占整个 helper 槽。由 (29)、(31) 得
\[
 \boxed{\displaystyle
 H_2+c_RZ_{\rm all}+c_Dv^6D_N
 \le[A_6(v)+\Gamma v^7]\mathcal C
        +Cv\mathcal F+Cv^8S}
\tag{32}
\]
对某个固定 \(c_D>0\)。确切地，两项价格之和至多
\[
 [A_6+\Gamma v^7]\mathcal C
   -[\gamma_\nu A_6+\Gamma v^7]D_N+O(v^9S).
\]
因此 \(c_D\) 来自真实的“总容量减正常占用”与非正常严格折扣，不是重复使用 (21) 已用于坏层的同一整份成本；(21) 在第 5 节已特意把两份固定余量分别留好。

二维的横积为空积、\(b_e=1\)；各源角只有这一轴对，核心剪裁和 (24) 第二分支为空。实际角位占用使用 R027 已审的二维唯一角槽及根权 incidence，四份容量仍为 (30)，误差不大于 \(O(v^3S)\)。第 2 节标量、同见证加法、(15)、长条容量和所有吸收完全适用。因此 (29)–(32) 是每个固定 \(n\ge2\) 的同一局部命题，不从 R042 的 \(n\ge3\) 声明直接越界引用。

## 10. 免费面共享、与 R043/R044 的接口及边界

等价地，对任意预先指定的 \(\chi>0\)，先保留上述固定 \(c_R,c_D\)，最后取 \(v\) 足够小，使 (32) 的 \(Cv\mathcal F\le\chi\mathcal F\)。故其同一联合接口也可写成
\[
 H_2+c_RZ_{\rm all}+c_Dv^6D_N
 \le[A_6(v)+\Gamma v^7]\mathcal C+\chi\mathcal F+C_\chi v^8S.
\tag{32'}
\]
\(\chi\) 仅影响起效阈值，不要求总装提供一个无法匹配的大固定面系数；例如可先按同一次恢复收入取 \(\chi<c_0/8\)。

式 (32) 右侧的 \(Cv\mathcal F\) 由 (7) 那**同一次**真实正恢复预算支付。随后全局小源盒外删除若再花 \(C'v\mathcal F\)，应共同要求
\[
 (C+C')v\mathcal F\le\tfrac12c_0\mathcal F.
\tag{33}
\]
保留的另一半 \(c_0\mathcal F\) 不被本稿虚构成第二次收入。内部图与边界价格仍按原固定顺序选定；本稿没有改变 \(I\) 或源递推。

若 R043 对同一正常集合 (5) 和同一 \(D_N\) 交付其好所有者高支撑 gross \(O(v^8S)\) 及坏根权
\(S({\rm Bad})\le C(D_N+v^{9/2}S)\)，而 R044 对同一最后大管交付其坏所有者净价，则可先固定
\(\eta<\min(c_R,c_D)/4\)，由其两份供货共同支付全部高支撑交量。这里只指出预算接头；那些新几何仍须各自独审，本稿不把口头状态当作它们已经成立。

本稿的新已证候选是 (29)、(32) 的全二轴共同价格；它只用本层显式列出的粗先验和实际几何，故对任意预先固定的有限末深度 \(N\ge6\) 可同层重放，特别包括新八层族。须先在新族上交付 (1)–(3) 及槽分配所需的原固定窗口输入，不能直接借用旧末层对象。不单独宣布 \(C_n\) 的全域八阶定理。其定量性不调用 R037 的无速率正常轮廓，也不用任何“六阶饱和遂自动给八阶”的跳步。所有固定参数顺序为：原有限末层及 CDF → 正常窗 \(\nu,K_N\) 和全域 \(c_\nu\) → 大窗 \(T\) → \(m,C_{\rm soft}\) 和共同核心、\(E_i(T)\) → 最后大管 \(K_{\rm all}\) → 预算份额 → 小 \(v\)。只允许最终 \(Cv^8S\) 常数依赖固定 \(T,m\)，用于后选 \(T\) 的吸收常数已在 (15)、(24)、(27)、(28) 明确保持独立。

## 11. 计算证据口径

正常紧窗二次强制性直接复用已实际成功的 R042 coercivity.wl/log（\(c=1/1000\)），不重跑相同计算。全域非正常强制性 (9) 已给解析反例序列证明，其常数不由有限采样推断。另已向 Terra 提交一个固定有理窗的独立严格折扣交叉核验；该有限输出未在本文冒充全域证明。本文的关键新几何是 (15) 的一次正部比较、(24) 的 uniform 短边帽比例及 (28) 两管记账，均由实际集合包含证明，不能由标量程序代替。
