# 第431轮　细带共同源的门槛位移不等式

- A 覆盖全部等权深度轮廓。第421轮选出的细带在自身门槛上的近最优损失至多为逆接收份额，且向下移动一个原稳定步长后仍有同阶控制。
- 对共同源 $E\cap D$ 得到显式的体积—接收份额二择估计；在现有步长下仍有平方根损失。
- 一般系数差仍为 $(\mu-1)/8$。

2026-11-13。设第349轮稳定门槛为 $p$，步长 $\delta=1/(\sqrt n\log^2 n)$，原目录 $|U|\ge(1-1/n)m H_n(p)$。第421轮带的下端 $q\ge p$，负责体积 $A=\theta |U|$。由单调性同时得到
\[
 mH_n(q)/A\le[(1-1/n)\theta]^{-1},\qquad
 mH_n(q-\delta)/A\le(1+O_I(1/\log n))/\theta .
\]
这两式不要求 $H_n(q-\delta)/H_n(q)$ 本身接近一。结合第426轮，已可在所选细带取得 $v_q\gtrsim_I \theta m\delta$；此前在未分带稳定门槛使用的估计由此获得带内版本。

为追踪共同源，沿第430轮记 $Y_i=|(E\setminus D)\cap R_i|/|R_i|$，$\alpha=v_qH_n(q)/A$，$\beta=mH_n(q)/A$。当 $Y_i\le y$，原行严格成功保证 $E\cap D$ 在同一盒对 $q-y$ 严格成功。因此对 $0<y<q$、$0<t\le1$，
\[
 1\le |E\cap D| H_n(q-y)/A+\alpha+t\beta+(1-q)\alpha/(ty).
\]
这给出对所有允许步长同时成立的下界，而非仅有依概率收敛。令 $y=\delta$。若前项至少 $1/2$，共同源已至少为 $c_I\theta m$；否则取 $t=1/(8\beta)$（$\beta\ge1$），得 $\alpha[1+8(1-q)\beta/\delta]\ge3/8$。于是 $v_q/m=\alpha/\beta\gtrsim_I\delta/\beta^2\gtrsim_I\theta^2\delta$。第二支甚至弱于上段直接传递，故仅凭共享修补及单调性不能声称改善了平方根指数。要取得次多项式损失，需要排除小共同源且赤字落在稳定步长以上的第二支，或在更宽门槛范围控制 halo 衰减。
