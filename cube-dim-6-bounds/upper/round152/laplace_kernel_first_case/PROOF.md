# C152准备：从正Gamma截断核到一个真实m2贡献

状态：PREPARED_THEORY_ONLY。Sol未运行任何求值、数值或科学预检；
单个新Arb证书由Terra在其独立clone经主与独审双gate后执行。

## 1. 原核的精确归一化

Y~Gamma(1/2,rate1)，h>0，p_h=erf(sqrt(h))。130原核为

 k_h(a)=exp(a) E[exp(-a^2/(4Y))*1(Y<=h)], a>=0.

原未截断期望为exp(-a)，故此处exp(a)不能漏掉；若h趋无穷，
k_h(a)趋1而不是exp(-a)。130/147已证明0<=k_h(a)<=p_h，且k_h
递减。另由Y<=h得到高斯包络
k_h(a)<=p_h exp(a-a^2/(4h))，保证下面所有实lambda的积分有限。

## 2. 独立正积分推导Laplace核

定义I_lambda(h)=int_0^infinity exp(-lambda*a)k_h(a)da。先对正的
Gamma/a被积函数使用Tonelli，再完成a方向Gaussian平方：

 I_lambda = (1/sqrt(pi)) int_0^h y^(-1/2)e^(-y)
              [int_0^infinity exp(-(lambda-1)a-a^2/(4y))da]dy
          = int_0^h exp(lambda*(lambda-2)*y)
                     erfc((lambda-1)*sqrt(y))dy.                 (1)

这一路径直接来自正截断Gamma表达式，不假定A、B独立，也不依赖
对147原程序的执行。设c=lambda-1，t=lambda*(lambda-2)=c^2-1。
在y>0上

 d/dy[exp(t*y)erfc(c*sqrt(y))]
   =t exp(t*y)erfc(c*sqrt(y))-c e^(-y)/(sqrt(pi*y)).

最后一项在0附近可积；积分并使用边界值1和
int_0^h e^(-y)/sqrt(pi*y)dy=p_h，得lambda!=0,2时

 I_lambda=[exp(h*lambda*(lambda-2))*erfc((lambda-1)*sqrt(h))
                   -1+(lambda-1)*p_h]/[lambda*(lambda-2)].       (2)

lambda=0,2不是发散，而是可去退化。直接由(1)积分erf(sqrt(y))：

 I_0 = h*(1+p_h)-p_h/2+sqrt(h)*exp(-h)/sqrt(pi),
 I_2 = h*(1-p_h)+p_h/2-sqrt(h)*exp(-h)/sqrt(pi).                  (3)

因此I_0+I_2=2h；lambda=1另有I_1=1-exp(-h)。这些是解析身份，
本包只数值认证lambda30，退化分支不额外运行。靠近退化点时若
未来使用(2)，应采用已认证的局部展开或足够区间精度，不能把
零除或严重消去直接当作结果。

## 3. 真实cone0、pair(16,17)的完整尺度

147真实prototype_result.json而非一个未绑定toy：原130 cone0的
alpha0=1/171992678400，零基body pair为(16,17)。四个零列以准确
概率质量1积分掉；余下

 A=E5/60+E6/30, B=E5/60, E5,E6 iidExp(1),
 E5=60B, E6=30(A-B).

Jacobian绝对值为1800，支撑0<B<A，密度
rho(A,B)=1800 exp(-30A-30B)。其质量准确为
1800/(2*30^2)=1。A、B相关，不能按独立输出乘边缘期望。

但exp(-30(A+B))*k_h(A)*k_h(B)是对称被积函数，故交换A、B
将完整正象限分成两个等积分三角形（对角零测），得到

 E_pair=E[k_h(A)k_h(B)]=900*I_30(h)^2,
 I_30=[exp(840h)*erfc(29sqrt(h))-1+29erf(sqrt(h))]/840.          (4)

h固定1/128。由130完整正表达式，本单个无序pair/单个cone的
mu2贡献为

 contribution=(46080/pi^3)*(1/171992678400)*900*I_30^2.          (5)

不额外乘496，不将900当作独立输出概率因子，不自动乘任何其他
pattern multiplicity。本量是完整beta2的一个正部分，不是完整
beta2的上端，亦不单独产生M或C6界。

## 4. 首个唯一数值证书的门控与验收

Terra固定Arb512，只求一次lambda30的完整signed分子，再除840；
全表达式向外包围，不能逐项裁剪。计算E_pair、公共pi因子及
(5)的准确rational端点。检查I_lower>0、I_upper<p_lower/30、
E_upper<p_lower^2，和I/E/部分mu2区间宽均<=10^-80。后两个概率
上界来自k<=p与rho质量1，不是经验截断。

只有一次sqrt(h)、erf、exp、erfc求值（各1），pi常数1；0逐点
k_h值、0数值quadrature、0Gamma/Bessel/MC/native、0旧343456
functional、0旧6D网格、0新density生成。lambda0/2代码分支仅
供读审，不在本有限首例中求值。固定一worker、nice>=10、线程1、
无wall/CPU限，排他结果/started标记，不自动重试或扩精度。

所有数值仅用于这个真实原型；完整9961模式的后续积分接口及
全域误差付账见GENERAL_2D_INTERFACE.md。下界150保持USER_PAUSED。
