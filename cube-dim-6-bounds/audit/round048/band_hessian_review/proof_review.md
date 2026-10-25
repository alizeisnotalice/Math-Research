# F048 正 band Hessian：解析与源码独立审查

## 1. 固定对象及 cone 归一化

只处理同一已认证 body32 fan、同一正 band 0<a<T。本批准确参数a=1/250，T=6/125，点为旧oldc。新frame、更换fan或只改q不在合法输入承诺内。

在一锥中令y_j=q_j u_j，alpha=|detR|/∏q_j，b_j=x·r_j/q_j。该锥密度贡献为(2π)^−6 alpha Re∏(t−ib_j)^−1。因为DD5[(t−ib)^−1]=i^5∏(t−ib_j)^−1=i∏，其时间积分的实部是

    DD5[f;b0,…,b5]，f(b)=atan(b/a)−atan(b/T)。

这也可由Re[−i/(t−ib)]=b/(t²+b²)直接积分得到。零点处按解析延拓，不能用sign(0)或删项。整个正时间band保证f在|Im b|<a的条带解析。coarea/cone权重和最终(2π)^−6未因求二阶导改变。

## 2. 有标签节点与全零/重复 DD7

取位于解析条带内、围住全部实节点的固定矩形contour，Cauchy差商表示为∮f(z)/∏(z−b_j)/(2πi)。对标签i求导插入一个相同极点。再次对同一标签求导时，平方因子求导产生2；对不同标签求导各插一个因子，系数为1。因此

    ∂ij DD5=DD7(base+bi+bj) (i≠j)，
    ∂ii DD5=2 DD7(base+bi+bi)。

不同标签数值相同也不能误乘2。此推导本身涵盖重复与全零；无需先要求节点不同再进行不受控极限。

归一化jet由f'(b)=a/(a²+b²)−T/(T²+b²)的Taylor系数产生。对单个u，分母乘幂级数给g0=u/(u²+b²)，g_k=−(2b g_(k−1)+g_(k−2))/(u²+b²)，g_-1=0；f的k阶normalized jet为(g_(k−1),a−g_(k−1),T)/k。代码正确除k而非再除k!。b=0时偶jet为0，奇jet为(−1)^((k−1)/2)(a^−k−T^−k)/k，全部保留。

重复base分支排序准确Rat并做八条目Hermite表：首尾准确相同则取相应normalized jet；否则以真实节点差相除。不同准确节点的Arb差若包含0则整体失败，没有丢项。Base全零、两重复零、四重零、六重同一非零和多个簇均在解析域内。

源代码使用atan((T−a)b/(aT+b²))计算f0。由于a,T>0，两个atan(b/a)、atan(b/T)同号且差严格位于(−π/2,π/2)，分母aT+b²>0，故主值恒等式无遗漏π分支。

## 3. 互异base快式

令w_i=∏_(j≠i)(b_i−b_j)^−1，c_i=f(b_i)w_i，S_i=Σ_(j≠i)(b_i−b_j)^−1。直接对重心和求导给

    d_i=f'(b_i)w_i−Σ_(j≠i)(c_i+c_j)/(b_i−b_j)。

再求导，mixed为(d_i−d_j)/(b_i−b_j)，diagonal为

    f''(b_i)w_i−2f'(b_i)w_i S_i+c_i S_i²
      +Σ_(j≠i)(c_i+2c_j)/(b_i−b_j)²。

源码每个无序node pair分别更新两边S的正负号及各diagonal的c_i+2c_j项；f''来自2·jet2。与上述公式逐项一致。此分支允许单个零节点，只禁止数值未分离的准确互异节点。出现重复base时不走快式，而改走DD7。

## 4. 真实物理链与全负号

b_i=Σ_k r_(i,k)s_k x_(sigma(k))/q_i是线性映射，物理Hessian为J^T Hnode J，没有第二阶节点映射项。存i≤j时，node diagonal贡献r_(i,k)r_(i,l)，node mixed贡献r_(i,k)r_(j,l)+r_(j,k)r_(i,l)，再除q_iq_j。该公式对物理k=l也自动产生正确双项；对k≠l不重复加额外2。随后乘s_ks_l并映射到sorted(sigma(k),sigma(l))。21字段索引0,…,20覆盖全部上三角，diagonal和mixed无缺项。

f为奇函数，DD5为偶函数；其node Hessian在全局负号下仍为偶。两个Jacobian负号相消。作者仅取最后sign为+的32代表、权2，恰覆盖64signs；未另作负号canonicalization。本独立全群参考显式使用64signs、权1，因此也交叉验证这项归一化。

## 5. 与作者不同的完整两分量参考

参考不构造node Hessian，不使用J^T Hnode J。直接对物理方向alpha,beta扰动b_i，记d_i^alpha=∂alpha b_i。定义

    A_i^alpha=Σ_(j≠i)(d_i^alpha−d_j^alpha)/(b_i−b_j)，
    C_i^(alpha,beta)=Σ_(j≠i)(d_i^alpha−d_j^alpha)(d_i^beta−d_j^beta)/(b_i−b_j)²。

w_i的物理一阶导为−w_i A_i^alpha，二阶为w_i(A_i^alpha A_i^beta+C_i^(alpha,beta))。因此一项f(b_i)w_i的二阶导准确为

    w_i [ f''(b_i)d_i^alpha d_i^beta
          −f'(b_i)(d_i^alpha A_i^beta+d_i^beta A_i^alpha)
          +f(b_i)(A_i^alpha A_i^beta+C_i^(alpha,beta)) ]。

alpha=beta直接使用同式，无任何事后倍数。参考的f使用两个单独atan，f'、f''使用直接有理公式，区别于作者atan商和递归jet。15个node pair分别构造物理差分比，再累积六个pole-product项，只输出x2x2、x5x6。

每分量实际遍历720 permutations×64signs×21cones=967680项，全部计算，没有global-negative配对；只共享两分量的节点与权重计算。最终同样乘alpha和(2π)^−6。准确输入固定oldc，没有新空间点。该独立参考若遇zero或repeated base或精度分离失败，按预登记整体停止；这不是作者引擎的限制。本点实测二者均0，故参考完成且无跳项。

Arb768所得两个完整区间均被作者Arb512区间完整包含。独立完整21分量参考没有做，始终false。其它19项的支持来自解析推导、源码核验和scalar锚，不能称为逐项全群独立实算。

## 6. 八项独立零/重标量参考

另选8个已保存锚，不调用作者anchors。对八个copied节点，按每个不同pole x的重数m展开∏_(y≠x)(x−y+z)^−m_y；系数准确有理。留数为normalized jet k乘该展开的m−1−k次系数之和。这是直接confluent pole-product方法，不使用作者Hermite三角表，也不使用其8×8 Hermite插值权重求解。

jet k≥1由符号微分u/(u²+x²)后除k!准确求值；仅atan值用Arb1024。8个参考涵盖全零diagonal和不同标签mixed、混合/四重零、互异快式diagonal/mixed、六重非零及多重复簇。每个完整区间被保存native512区间包含。准确留数权重还逐项核对对0至7次单项式的作用。

## 7. 输入、整数和ABI范围

C048 parser要求规范int64数字、正分母、0<a<T、21cone、q≤10^6、ray坐标绝对值≤10^6、positive det≤10^9、无额外token。物理点int64最小值拒绝，避免取负溢出。六项int128点积的绝对值至多6·(2^63−1)·10^6，分母积至多(2^63−1)·10^6，均远小于2^127；最终必须再通过int64门控。Rat排序的两个int64乘积也各处于int128范围。ray二次系数至多2·10^12，在int128中生成再过int64检查。

这些格式/范围检查不证明任意输入fan的覆盖或q真实性。数学验收只针对已锁7a57caad…旧fan；wrapper的哈希门是可信输入契约的一部分。不能把本引擎直接套到尚未完成全fan的新frame。

native以真实FLINT3.6.0 headers编译，static_assert检查ulong/fmpz8bytes、arb_struct48bytes，并在运行时检查版本。source/binary/header archive和三个runtime dylib哈希均吻合。独立参考也用相同已锁headers/runtime并检查该ABI；不声称跨FLINT版本或任意ABI兼容。所有数学量用准确有理整数或Arb，double只用于时钟和报告导航显示。

端点通过arb_get_interval_fmpz_2exp完整导出，比较使用准确Fraction，未从短radius字符串重造。旧band、6grad、p_a、p_T共9字段与C016原件完整dyadic逐字相同；原Euler区间含0。这是旧接口一致性验证，不能推导全部Hessian正确或二阶Euler。

## 8. 验收范围

接受本matching positive band Hessian的解析/源码接口与有限证据：8独立scalar关系、2个独立全群component。工具没有盒上余项或全域oracle，没有全Green Hessian替换，没有新grad(p)/二阶Euler，也没有Cn改善。
