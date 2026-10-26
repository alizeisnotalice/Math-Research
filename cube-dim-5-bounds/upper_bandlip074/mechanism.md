# 074 冻结机制：已审全frame Lipschitz比 + 正时间段均值

固定五维alpha=1/body16/Fourier逆因子(2pi)^-5；唯一中心为067 low。
新eps=1/512。只消费067已存89段中的第16至88段，共73段，严格覆盖
[1/512,infinity)；中心H直接取067 checkpoint第16个端点（零起算）。
无新H、stencil、jet、锥或精度改变。076未锁结果完全不消费。

固定三个真实裁片h=1/64、1/128、1/256，定义为
{1}×prod[max(0,c_j-h),min(1,c_j+h)]。h1/64两坐标裁切，h1/128
一坐标裁切，h1/256对称。所有矩、体积、h_eff以真实端点准确计算。

075已锁并独审的定理给所有x,u,t>0：
 |D_u p_t(x)| <= ||u||infinity p_t(x)/t，
 p_t(x+u)/p_t(x) <= exp(||u||infinity/t)，反向同型。
此定理由独立Cauchy单坐标方向导数的纤维转移、±body覆盖cube全部顶点而来，
不取Cauchy矩。只读报告SHA为
3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225。

另保留067基N的旧合法比
 R_old(h,t)=(sqrt(1+h²/(20t²))+h/(2sqrt5 t))^10。
实际R(h,t)=min(exp(h_eff/t),R_old(h,t))，不得假设exp对所有参数更小。
R两支均递减，整片以及连接中心的线段适用。每段[a,b]使用左端a。

令mu=E u、eta=E(uᵀQu)，Q=I4/4-J4/36。
072证明的五基密度Hessian上界仍使用(21/4)uᵀQu/t²。
线性项同时有两条已审score界，故lambda=min(||N^-1(0,mu)||1,||mu||infinity)。
真实片均值满足
 mean B_j <= m_j C_j，
 C_j=min{R_j,1+lambda/a_j+(21/8)R_j eta/a_j²}。
这里m_j为保存正时间段中心质量，未重新求值。对称片lambda=0。
各片以正m_j合成L=sum m_j/R_j、U=sum m_j R_j、M=sum m_j C_j，再用
 volume [L^q+(U^q-L^q)(M-L)/(U-L)]，q=5/4，认证真实连续power积分上端。

预登记两门：chord_upper/(volume*H_eps(center)_lower^q)≤11/10、≤21/20。
这仅比较中心常值基线，不把它叫真实积分相对误差，不作为全域I/g/C5证书。
三个配对均不同于072旧h/eps组合；本批不重新运行067或072生产器。

单程序wall预算180秒，170秒软检查点、175秒SIGALRM硬停止，留5秒收尾。
每完成一个时间段原子保存前缀，fingerprint锁源码、登记和输入；同源resume。
固定3片、73段、384-bit Arb；成功或失败都保存，不自适应加参数。
