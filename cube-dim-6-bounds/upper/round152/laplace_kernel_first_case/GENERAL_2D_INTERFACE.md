# 接149 full compiler的完整beta2二维积分接口（本轮不执行）

本件不再生成density，不读取149正在运行的输出作理论来源。
消费147已证的正pushforward与有限signed exponential-polynomial
形式，待149按其完整封存/审计接口交付9961 canonical patterns。

## 1. 覆盖、权与目标

每个原130 cone/pair由iidExp(1)^6的两非负线性形式产生联合law。
完整beta2为(46080/pi^3)*sum_C alpha_C sum_(i<j) E[k_h(A)k_h(B)]，
h=1/128。149须提供10416个原cone/pair到9961个key的准确覆盖，
每个原项恰一次。输出交换或iid列置换可以合并，但只保留ray
方向而丢掉列尺度不保持原law；任何缩放都必须显式付变换尺度。
令w_k为该key的全部正alpha/multiplicity之和。应核对
sum_k w_k=496*sum_C alpha_C；公共46080/pi^3只乘一次。

完整真实beta2上端严格小于7/10^8足够通过146固定候选的g0模型
门，仍需同模型globalg0；小于8/10^8足够改善128。这里只规定
目标，当前单个m2原型不等于完整beta2。完整链不得使用144乐观
bL代替真实上端。

## 2. 支撑ray变换与单项资格

对一个已认证扇区，取其准确非负独立边界rays u,v，令
(A,B)=u*s+v*t，s,t>=0，Jacobian=|det(u,v)|。对每一完整密度项
P(A,B)exp(-lambda.(A,B))，准确展开

 P(u*s+v*t)=sum_(p,q)c_pq*s^p*t^q,
 d1=lambda.u, d2=lambda.v.

必须逐项、逐扇区认证d1,d2的符号。lambda坐标有负数不自动
有问题；沿其真实支撑rays的衰减才是资格。147已审m3/m6项均有
正边界衰减；不能据这几例替全9961模式预付。

若d1,d2>0，每个monomial积分可准确写成

 c_pq*|det|*p!/d1^(p+1)*q!/d2^(q+1)
   *E[k_h(u1*S+v1*T)*k_h(u2*S+v2*T)],

S,T独立Gamma(p+1,rate d1)、Gamma(q+1,rate d2)。Gamma测度本身
是正概率；c_pq仍可能为负。所有这些signed系数必须原样合成，
误差按绝对系数/normalizer付账。这里没有宣称每个signed项是
正density，没有取正部，没有丢负项或把部分展开重归一化。

## 3. 正Gamma权的二维单调盒与完整尾

因rays非负、k递减，K(s,t)=k(A)k(B)沿s,t递减。在Gamma矩形
[s0,s1]x[t0,t1]上，用K的upper-corner可靠下端、lower-corner
可靠上端，乘该矩形真实Gamma概率区间即可围住完整积分。
shape为整数，CDF由有限exp-polynomial式给出，不能用抽样替代。

截断[0,S]x[0,T]后余概率准确为
1-F_(p+1,d1)(S)*F_(q+1,d2)(T)，尾积分在
[0,p_h,U^2*该余概率上端]。Gamma整质量始终1；有限矩形概率
从不重归一化。可按全域加权width选择下一盒并沿明确轴中分，
每批只增加已登记的新点；原130的343456个6D endpoint不重跑。
二维化只解决维度，尚不保证某一有限盒预算有竞争力。

## 4. 完整density矩形质量：处理signed展开/非正ray decay的另一接口

若某项d<=0，不能把它强行解释成正rate Gamma，也不能在无穷域
逐项套不收敛的Laplace质量。可在有限ray矩形上先完整计算真实
density质量，再乘K的单调角值。

对整数p，J_p(d;l,r)=int_l^r x^p exp(-d*x)dx：d=0时是准确多项式；
d!=0时用
J_0=(exp(-d*l)-exp(-d*r))/d，
J_p=(l^p exp(-d*l)-r^p exp(-d*r))/d+(p/d)*J_(p-1)。
完整矩形质量为所有密度项、所有c_pq的signed和
|det|*sum c_pq J_p(d1)*J_q(d2)。应先得到该完整和的区间，再
使用已证明的正density/概率性质与[0,1]相交；不对单项裁零。
区间过宽表示消去/精度问题，不能据浮点符号删项。

此法的尾必须由完整pushforward支付。例如Z=A+B=sum_j d_j E_j
且d_j=b_j+c_j>=0。任选准确theta>0满足theta*max d_j<1，有
P(Z>L)<=exp(-theta*L)*product_j(1-theta*d_j)^(-1)。ray域s>S或
t>T分别蕴含Z>||u||_1*S或Z>||v||_1*T，故用两个Chernoff尾
上界的和乘p_h,U^2可付完整遗漏概率，且最多截到1。不能借用旧
6D uniform-t尾而不证明其对新域的覆盖。kernel本身的Gaussian
包络也可进一步加强尾，但不能免去完整概率/符号账本。

## 5. 9961模式与公共pi的统一误差

全域建议预算epsilon_beta<=1e-9，全部以真实beta2单位合成：
尾2.5e-10，有限盒/积分余项5e-10，density/CDF/kernel算术及公共
pi因子2.5e-10。这是待实际反馈的规划，不是已取得的误差。

若各pattern完整期望区间宽为delta_k，kappa=46080/pi^3的区间
为[kL,kU]，则总width可用

 kU*sum_k w_k*delta_k +(kU-kL)*p_h,U^2*sum_k w_k

统一控制；各pattern非负且<=p_h^2的资格来自完整density，不是
signed展开项。Gamma signed-term法的delta_k必须先按全部
|c_pq|及normalizer累计，再与pattern权合成，负系数项区间方向
不能反用。任何resonant多项式、极限端点、扇区边界和重数都由
完整compiler语义负责，不能漏负密度项或遗漏任何原pair。

一个可携带的分配为delta_k<=epsilon_int/(kU*9961*w_k)，其中
w_k>0，全部取准确有理保守端；公共pi与尾留独立额度。若有限
主体严格小于6.9e-8、全部余项不超过1e-9，则完整上端严格小于
7e-8。真正验收直接比较完整beta2准确上端与7/10^8，不从局部
pattern成功、原型小值或显示小数推断globalg/C6。

## 6. 调度接口

先只执行C152单m2闭式Arb首例。成功后，主可依据149完整封存
metadata选择一个已审非可分离pattern测试上述一般方法，再按
实际weighted width/消去成本登记后续有限批。不得现在重编译
9961模式、重复149那100项、启动盲2D/6D网格或自动扩worker。
