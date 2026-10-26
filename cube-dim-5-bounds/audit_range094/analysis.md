# 094 独立解析审阅

## 固定模型与 091 比值

审阅对象为未单位化 body16 方向 v=(1,±1,±1,±1,±1)、共同 Cauchy 时间 t、α=1、逆 Fourier 因子(2π)^−5。固定 B、K=ker B、A=Bᵀ/16，BBᵀ=16I，fiber Jacobian为1024。075 的逐点正有限 fiber 表示是共享已接受解析输入；此处独立核后续 lift 操作。

若 Ba=u，则 a−Au∈K。对固定 u 平移积分变量 w，Hausdorff 测度不变，故 p_t(x+u)=1024^−1∫_K f_t(Ax+w+a)dw。无需 a=Au，也无需随 u 连续选择 a，更无需 Cauchy 矩。ρ=(sqrt(1+a²/(4t²))+|a|/(2t))² 满足 (ρ−1)²t²=ρa²。二次式 ρ[t²+(z+a)²]−(t²+z²) 的判别式准确是 4[ρa²−(ρ−1)²t²]=0，首项ρ−1≥0；a=0单独给相等。因此 s_t(z+a)/s_t(z)≤ρ。将位移反向得下界，再对产品正积分即可转给 p。

u1=0 时，尾绝对值排序 w1=r≥w2≥w3≥w4≥w5=0，符号σ、同值按索引。取v0=(1,−σ)，vk依次翻转排序前k个尾坐标，a0=−r/2、ak=(wk−wk+1)/2。首坐标Σak=0；排第j个尾坐标准确为 σπ(j)[r−Σk<j(wk−wk+1)]=uπ(j)。零和并列只产生零 gap，不破坏提升；r=0全系数零。

四个 dk=wk−wk+1≥0、Σdk=r。asinh 在非负轴凹，故产品 log 比上界为 2asinh(r/4t)+2Σasinh(dk/4t)≤2asinh(r/4t)+8asinh(r/16t)。这与只检查盒角点或密度凸性无关。双侧常数R(r/t)对r/t递增。真实盒任意位移u=x−c和s∈[0,1]有||su||∞≤h∞，故同一R(h∞/a)支付所有anchor-to-box线段及全部t≥a。正时间积分保留下上界，无穷尾由既有t^-5峰值可积保证。此步骤没有改变eps、anchor、face体积、坐标或时间权。

R(z)<e^z 对z>0来自asinh(s)<s；若比较匹配的旧Q-asinh，准确平方完成uᵀQu≥(5/24)||u||∞²、Jensen给logR≤10asinh(z/10)，而4sqrt6<10，故091的严格解析比较正确。但092安全性无需数值假定新RU小于每个旧RU；直接保留旧aggregate range/mean/self证书并取交即可。

## 085 均值接口

独立审阅085证明：最小l1提升||a||1=r，且u1=0时Σa=0，所以|a_i|≤r/2。令b_i=|a_i|/r，Σb_i=1、b_i≤1/2；加权Cauchy–Schwarz将产品二阶导上界归结为Σb_i(S_i²+b_iT_i)。单坐标有

  gβ(q)=((4+2β)q−2β)/(1+q)² ≤ K(β)=(β+2)²/[4(β+1)],

差为[(β+2)q−(3β+2)]²/[4(β+1)(1+q)²]；K递增且K(1/2)=25/24。正fiber积分及有界score的支配求导给 D_u²p≤(25/24)||u||∞²p/t²。Taylor准确积分余项的1/2给25/48；先平均线性项，再用|D_μ B(c)|≤||μ||∞m/a，不得改成平均|u|或把偏心anchor当中心。因此 C=min(RU,1+||μ||∞/a+(25/48)M2 RU/a²) 合法。M2是同一真实盒相对同一anchor的E||u||∞²。

094消费085/088已接受的全部M2资格和其冻结逐节点记录，不重新计算旧成功moment实验；新增消费仍逐节点核h∞、μ和记录身份。此复用边界在结果中明确。

## 量化、凸弦与层级

每band准确输入是 mL=H_a,L−H_b,U>0、mU=H_a,U−H_b,L。尾段b=∞的H_b=0。mL/RU为下界，mU RU为上界，mU C为均值上界；不是有限设计点近似积分。新range和mean分别和同anchor/真实盒的旧证书交，保持0<L≤M≤U。函数f(z)=z^(5/4)递增凸，盒平均≤[(U−M)f(L)+(M−L)f(U)]/(U−L)，U=L时取常数；方向不是Jensen上界。self上界可取旧self与新weighted chord之min，即使新range未包含旧range，两个积分证书也可取min。

所有canonical节点表示完整coordinate permutation orbit；children合并完整parent orbit。有效U=min(selfU,Σchild effU)、L=max(selfL,Σchild effL)合法。数值最优cover可混用历史parent和children，但physical frontier仍独立定义几何覆盖；不得把两份都加进总积分。每个有限更新prefix中未更新self继续付旧088区间，故不必等所有节点改善才拥有完整证书。094沿用087已证的相同082空间闭覆盖，检查初末全部edges/frontier完全不变并独立核volume=1，不重跑旧边界例。

## 全面 I 与 C5 门

076的E是整个正面完整短时间L^(5/4)范数。G=H+R为扩展非负函数；奇异线零面测不阻止Minkowski给 I≤(UH^(4/5)+E)^(5/4)。本任务覆盖全部F，不额外加064tube。若B_U≥UH^(4/5)+E，则(B_U−E)^5≥UH^4精确检验其方向；I_U^4≥B_U^5及cost_U^4≥(5B_U/g)^5同理。g=1/10000属于主既有完整接受输入；094不重生产g。32门为5^5B_U^5<32^4g^5。显示cost下降不构成该门通过，也不能据未过门推断真实C5≥32。

## 代码范围与风险结论

实际逐行审阅091 ratio.py/report及092 prepare.py/engine.py。正平方根的整数取整方向、正整数幂与最终向外量化正确；192-bit表示证书网格，不冒充原H更高精度。092逐node只读取旧H，采用新RU并作aggregate交，代码与mechanism一致。partial block失败时state不推进，冻结完整block才更新self，保留完整付款。任何进程崩溃后的自动恢复并未作为本轮实验声称。

目前解析未发现足以阻断092完整I的gap。全数值结论必须由另行冻结的whole_scope与实际whole_results支撑，不能以作者自报pass替代。原始H数值内核及085/088 M2生产属于共享接受输入；094独立的是091解析审阅、89项新ratio代数资格、全节点新包络/凸弦及14prefix合成。
