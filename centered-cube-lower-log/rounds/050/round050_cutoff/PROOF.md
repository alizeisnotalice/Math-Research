# 有限尺度周期比较到Lebesgue可积截断的转移

先实际运行probe.py：一维正周期阶梯密度、三个共同尺度，对N4/8/16直接算截断窗口与边界积分，观察到归一化weak profile差按1/N缩小。数值网格不作向外证书，以下证明独立于它。

固定n，任意非负Z^n周期函数f，其一个周期积分1。令Λ为J个正的全边长尺度，全部L≤2。N≥2为整数，F_N=f·1_{[-N,N]^n}∈L¹(R^n)，质量V=(2N)^n。记Mper=max_Λ A_L f，Mcut=max_Λ A_L F_N；周期weak范数在T^n上，截断weak范数用Lebesgue测度并除V。

在core C=[−N+1,N−1]^n，每个中心窗口半边长≤1，所以A_LF_N=A_Lf，除边界零测集外处处成立。C各边长为整数个周期，因此任意λ>0的core超水平集体积正好是(2N−2)^n倍周期超水平集质量。记r_N=(1−1/N)^n，Wper=||Mper||_{1,∞}，Wcut=||Mcut||_{1,∞}/V。

对core之外，有限max≤sum及Fubini给

∫_{C^c}Mcut≤Σ_L∫F_N(y) A_L1_{C^c}(y)dy。

若y∈[−N+2,N−2]^n，半边长≤1的窗口完全位于C，因此指示平均为0；其余位置该平均≤1。整数周期积分故有

∫_{C^c}Mcut≤J[(2N)^n−(2N−4)^n]≤2JnV/N。

不要求f有统一高度界；仅用非负质量、Tonelli及Bernoulli不等式。n1或N2时内部为空也成立。

同一严格阈值的core、外部拆分给

r_N Wper≤Wcut≤r_N Wper+2Jn/N。

没有使用错误的weak-L¹函数和三角不等式；外部只用L¹付款。又周期平均质量1，有限并集Markov给Wper≤J，1−r_N≤n/N，所以保守地

|Wcut−Wper|≤3Jn/N。

在第50轮自然增长网格J=O(Kv)、n=Θ_p(K²v)、N=n²时，3Jn/N=O(1/K)=o(1)。该界对所有质量1周期来源一致，故同轮的mixture≤best periodic product+o√n，可以转为这些真正Lebesgue可积截断的归一化weak范数比较，额外误差o(1)。所有bank及混合均用同一个N，质量V相同。

这里只转移有限共同尺度集合，不恢复连续尺度，也不提供各product来源的一般维数上界。没有循环调用待证中心极大弱型界。
