# C146：固定144准确raw系数、完整H认证与合法成本

状态：PREPARED_UNRUN，主与独审双gate后仅认证一个新候选。144
乐观目标低值不是M上包。本包不重跑任何旧目标、kernel或优化。

## 1. 准确系数与一般cost

固定144保存的49维2^-44候选x。定义tilde_beta为其优化时使用的
正参数表：除h=1/128,k=2外均为128原合法矩上包，该一项等于
130真实mu2下端bL。定义a_i=x_i/tilde_beta_i为准确rational商。
a并不一般是dyadic，绝不对a量化或重新归一化。

令W_a(N)=a0+sum_(h,k)a_hk Phi_k(N_h)，H(a)=E W_a^(-1/5)。所有
a>=0、a0>0。对真实det权d，一般Holder给

 M6 <= [E d W_a]^(1/6) [E W_a^(-1/5)]^(5/6),
 M6^6 <= C_U(a) H(a)^5,
 C_U(a)=a0*beta0_U+sum a_hk*beta_hk_U.

若随后消费123既有完整时间/空间接口，
I^5 <= (3/16)^6 * C_U(a)*H_U^5 /384^5；cost只出现一次幂。

这里真实矩满足E[d Phi_k(N_h)]<=beta_hk_U，cost并未假设为1。
当前合法beta2取128原完整上包和130完整新上包中较紧者，即原
128/124 cap约1.83312547745e-7；130约4.66948892556e-7没有改善它。
其余合法cap等于tilde_beta。故changed index32的x记为x_b时，

 C_U(b)=1-x_b+x_b*b/bL,
 x_b=10680174170037/17592186044416.

preflight另用全部49个准确a逐项计算sum a_i*合法cap_i并验证该
恒等式。bL只定义系数；它不被当作真实cost cap。若未来量化a，
必须按实际量化系数重算cost，本包不执行该操作。

## 2. 完整分布与新H的认证

消费128已认证的nested-count分布接口与固定Joint程序：8个h的
两套完整CDF端点law各有9个正类别概率，准确和为1；不将CDF
区间独立混成类别。Phi非负非减，a非负，故W随count非减。
prefix为全部735471个0<=n0<=...<=n7<=16状态，余K=n7=17..32
的16个binomial尾行均完整计算。

下端：pU完整law下正prefix的和；超出prefix取0仍保持可靠下界。
上端：pL law下完整prefix，加上全部16尾行；尾行丢弃前七组
非负项，只保留a0和最后h组，负五分之一次幂因此只会上升。
没有概率删尾或重归一化。

原128 Joint读取准确a后仅把W的值向下/向上包围到192位整数
网格，再对W^-1/5产生整数端点vL,vU，逐项准确检查
vL^5*WU<=2^(192*6)<=vU^5*WL。这是数值函数值的向外包围，
没有把实际a改成dyadic；cost仍用准确商a计算。

15个独立worker各处理不相交连续prefix索引区间；前6片49032项，
其余9片49031项，共735471。每项保存新root整数对、状态范围/
顺序hash、完整概率质量及加权整数和。worker逻辑root请求总数
735471；实际Arb计算可因同一新W的局部cache命中减少。pool退出
后collector处理16尾请求，因此并发科学worker至多15。所有结果
用旧准确概率分母2^6144和加权分母2^6336合成新H_L,H_U，核对
两套prefix+tail概率质量均准确为1与完整状态顺序hash。

## 3. 真实beta2上端门槛

记Q_U=H_U^5。对g0门Theta=(8/9)^6*24576^5/65000^6，充分门为

 b < bL/x_b * [Theta/Q_U-(1-x_b)].

改善128完整Mu^6的充分门将Theta替换为保存128准确Mu^6。
collector输出两个门的准确rational值、当前合法cap的准确比较
和M候选六次方C_U(current_b)*Q_U。也输出以H_L^5替换Q_U的
相应必要门；两端只涉及该固定raw候选，不是全族最优门。

当前合法cap先实际测试，不假定真实beta2已改善。144原浮点H^5
约2.50321e-8，乘当前合法cost后只约1.33958e-7，表明这个候选
预计需要明显更紧的真实beta2上包；浮点值不作任何准确判定。
若当前合法成本下不竞争128，完整证书仍用于确定下一积分所需
真实上包，而不追加盲优化或kernel批。

## 4. 固定预算与停止

一次主launch，15worker、每worker单线程、nice>=10，无wall/CPU
上限。新候选逻辑root请求恰735471+16=735487；Arb192，root
局部cache上限沿用65536，15份prefix压缩证据各少于100MB。无
优化、旧baseline、kernel/erfc/heat-kernel/native/MC/MATLAB调用。
任何错误停止且不自动重启，不擅自增候选/精度/状态域。主与独审
必须批准准确代码、输入与源hash；尚无release或science marker。
本候选认证不证明globalg0，未配合法全域g前不宣布新C6界。
