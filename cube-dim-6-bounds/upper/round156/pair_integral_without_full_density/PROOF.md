# C156准备：不展开完整density的真实联合律上包

状态：Sol仅推导/静审，未求值、未运行prepare或科学预检。Terra
执行一次有限诊断。目标是完整真实beta2的7e-8/8e-8门，准确门
来自146，不用显示小数构造来源。没有全9961 density再生成。

## 1. 单个真实联合律与flat seed

固定130原锥C和无序body pair。令v_j=(b_j,c_j)>=0，E_j iidExp(1)，
(A,B)=sum_(j=1..6)v_j E_j。零列以概率1积分掉；重复列/方向
保留其独立E，不合并成错误的单Exp。下面也允许整个law秩0/1。

选两个独立列i,j，D=|det(v_i,v_j)|>0，其二维density为
D^-1 exp(-E_i-E_j)在正确的column cone内、外为0，故L-infinity
不超过1/D。其余独立概率measure卷积不增加L-infinity。由于
k>=0，I0=int_0^infinity k(a)da，

 E[k(A)k(B)] <= min(p^2,I0^2/D).

最佳两列flat族取D最大值；若无独立列只保留p^2。任何全支撑
线性坐标变换也同时改变density sup与kernel的L1尺度，不能凭
投影/缩放白赚一个det改进。

若四个不同列分为两个独立pair，D1,D2>0，pair density的L2范数
平方准确为1/(4D)。Cauchy--Schwarz直接给

 ||rho_pair1*rho_pair2||infinity <=1/(4sqrt(D1*D2)).

剩余概率卷积仍不增sup。这是正卷积，不生成完整六列density。
四列subset中只需在三个不相交pairing中最大化D1D2。

## 2. k_h全半轴logconcavity的独立闭合

令s=sqrt(h)>0。原130 closed式
k(a)=[erfc(a/(2s)-s)-exp(2a)erfc(a/(2s)+s)]/2，a>=0。
逐项微分时两个Gaussian导数项准确抵消，得到

 k'(a)=-q(a), q(a)=exp(2a)erfc(a/(2s)+s)>0.

该式在a=0取右导数，k(0)=p=erf(s)，k'_+(0)=-(1-p)。不在
Gamma截断表达下对a=0盲交换导数。原正表达给
k(a)<=p exp(a-a^2/(4h))，故k趋0且可积；q的Gaussian尾亦可积，
因此k(a)=int_a^infinity q(t)dt。

下面独立证明erfc的log凹性。设T(x)=int_x^infinity exp(-t^2)dt，
R(x)=exp(-x^2)/T(x)，x>=0。对x>0，因t>x，

 T(x)<(1/x)int_x^infinity t exp(-t^2)dt=exp(-x^2)/(2x).

于是R(x)>2x；直接微分R'=R(R-2x)>0，x=0时仍有R'(0)=R(0)^2>0。
所以(log erfc x)'=-R(x)递减，erfc在[0,infinity)严格log凹。
q的log为线性函数2a加一个仿射正参数下的log erfc，故log q凹。

survival步骤也直接证明，不靠未核条件的Prekopa。记L=log q。
对u>=0，L(a+u)-L(a)随a非增；故

 k(a)/q(a)=int_0^infinity q(a+u)/q(a)du

随a非增。q/k随a非减，而(log k)'=-q/k，故log k在整个[0,infinity)
（0取右导数）凹。所用比值均正、有限，收敛已由Gaussian尾付清。

令f=k/p，g=log f，g(0)=0。当a+b>0时凹性给
g(a)>=a*g(a+b)/(a+b)、g(b)>=b*g(a+b)/(a+b)，所以

 f(a+b)<=f(a)f(b), a,b>=0.                                 (1)

a=b=0时两端都为1，因此包括全部边界情形。

这里不是把logconvex的Gamma y^-1/2权误用为logconcave混合。

## 3. 两个固定切线与完整联合MGF

对log k的切线可写k(a)<=C_r exp(-tau_r a)。只固定两个锚点：

 t0=0: tau0=(1-p)/p, C0=p；
 t1=s: Q1=exp(2s)erfc(1/2+s),
       K1=[erfc(1/2-s)-Q1]/2,
       tau1=Q1/K1, C1=K1 exp(tau1*s).

K1=k(s)>0，tau1=-k'(s)/k(s)>0。全半轴凹性使两条切线在所有
a>=0合法，不只在锚点附近。允许A、B分别选r,t in{0,1}，用原
独立E的完整MGF准确积分，得四个全联合上包

 G_rt=C_r*C_t*product_j(1+tau_r*b_j+tau_t*c_j)^(-1).          (2)

没有假定A、B独立，没有有限截断或遗漏概率。单列因子也有

 T_j=E[f(b_jE_j)f(c_jE_j)]
 <=S_j:=min(1,min_(r,t) C_r*C_t/
              [p^2*(1+tau_r*b_j+tau_t*c_j)]).               (3)

零列直接S_j=1，避免除零或虚构density。

## 4. seed乘剩余正因子

对任一选定source列集合J，由(1)及各组E的独立性，

 E[k(A)k(B)] <= E[k(A_J)k(B_J)]*product_(j notin J) S_j.

J为两独立列时用第一节min(p^2,I0^2/D)；J为四列时用
min(p^2,I0^2/(4sqrt(Pmax)))。J为空时为p^2*product S_j。
每raw pair最多1+15+15个seed上包，再与(2)四个上包取min，
最多35个明确有限合法候选。不搜索连续锚点，不扩大失败网格。
全支持、零列、重复方向、秩1输入均已覆盖；不成立的独立seed
跳过，始终有p^2与MGF因子界可用。

## 5. 全域alpha与预算方向

完整真实beta2 <= B := (46080/pi^3)*sum_(21 C) alpha_C
                   sum_(496 i<j) min(上述有限caps).

130原group-reduced21fan与全部pair-sum一起付46080群因子。每个
单项只是该约化表达式内的项，不是原始固定pair在全R6的事件。
全部10416项恰一次，核对46080 sum alpha=原R0；公共pi只乘一次。

Terra对B产生向外区间[BL,BU]。BU是合法真实beta2上包；BL仅是
这个指定上包表达式的下界，不是实际beta2/M/C6下界。若BL>=8e-8
可严格排除本固定seed/两锚点族达到8e-8预算；不排其他切线、
条件卷积、全density积分或真实矩更小。若BU改善146保存合法cap，
只先登记新的beta2候选，C6仍需完整M消费及globalg资格。

## 6. 唯一低成本诊断

只新增I0闭式一次、k(s)一点和其log切线：复用128准确p区间、
130准确kappa区间；新sqrt(h)、pi常数/sqrt(pi)各1，exp共3，
erfc共2，不新算erf或I30。I0=h(1+p)-p/2+s exp(-h)/sqrt(pi)。
四列P的sqrt按相同准确P缓存；若I0_L^4>=16P*p_U^4，则该seed
的min已准确等于p^2，直接省去sqrt。最多15*10416个不同P请求。

固定Arb512、一个nice10 worker、线程1，无时间/CPU cap。0新
density/SymPy/MC/native/6D盒/旧343456functional/旧基线调用。
所有pair搜索与MGF/种子预算在精确hash双gate后一次完成。
