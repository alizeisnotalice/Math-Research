# 175 有理细化斜率的独立审计准备

原112同一准确range，U>0、0≤L≤U；s是φ(x)=x^(5/4)准确弦斜率，
L=U时取导数。令a=U^(1/4)、t=(L/U)^(1/4)，则
`s=a R(t)`，`R=(1+t+t²+t³+t⁴)/(1+t+t²+t³)`。
几何级数约分只在t<1用作推导；延拓到t=1为R(1)=5/4，
t=0为R(0)=1，均无除U−L。D=1+t+t²+t³≥1，
`R'=t³(4+3t+2t²+t³)/D²≥0`。
所以0≤y≤a、0≤z≤t给yR(z)≤s。

对B=2^64、N=2^160，独审直接二分有理不等式，而不形成producer的floor后isqrt输入。
对整数k，`k⁴≤floor(n/d)`当且仅当`k⁴d≤n`；因此最大可行k恰是需要的整数四次根。
`root_ratio`先用bit长度给严格不可行上端，再二分，并保存left≥0、next>0。
Y使用n=Ui·2^256、d=N，Z使用n=Li·2^256、d=Ui，两个新根各有两条交叉四次幂夹逼。
左残差为0合法，next必须严格正；Li=0返回Z=0，Li=Ui可返回Z=B。
实际固定输入仍严格Li<Ui，不强行扩大旧strict-active集合。零差mean亦合法。
Y的整数部分不要求限于64位，B约束的是小数格点。

正的齐次Horner多项式给D_int=B³+B²Z+BZ²+Z³及P_int=B⁴+B³Z+B²Z²+BZ³+Z⁴。
新下端为Y·P_int/(B²D_int)，取gcd约分并与已认证s164交叉乘比较；
其max仍是同一准确弦下端且≥1/16。旧s164只复制，不重跑其root/threshold witness。
总gain仍相对原112 mean/self一次计算，独立gcd约分及shift/subtract长除保存商和余数。
不能从166self重复扣相对112的总gain；新self须≥selfL且≤self166。
旧factor、mean、H/R/J/log、旧端点power、旧gain和旧DFS均不重跑。

误差预算独立核对：D²的系数为1,2,3,4,3,2,1；减去
`t³(4+3t+2t²+t³)`恰余`1+2t+3t²>0`，故0≤R'≤1；单调性给1≤R≤5/4。
两个向下根格点满足0≤a−y<1/B、0≤t−z<1/B，故
`0≤s−yR(z)=(a−y)R(t)+y(R(t)−R(z))≤(5/4+y)/B`。
又`R−(1+t⁴/4)=t⁴(4−D)/(4D)≥0`，因为0≤t≤1时D≤4。
这些结论保留端点和perfect根等号；没有重跑主已注册的SymPy两残差，
也不能据解析误差预算预称有限floor支付/完整树给严格Cn改善。

完整原97972节点/52360根沿已核geometry/orbit/完整16-child分割，
更新97971个active self，唯一nonactive self保留但孩子仍递推。
独立递归DFS检查全部父子、根集合、无重/无环及全部节点覆盖，
取min(newself,childsum)，仅根求和；逐项≤166的保存self/eff。
完整同模型n5、16未单位化body方向、Lebesgue、共同Cauchy时间、alpha1、(2pi)^-5、H512、
S1/S2/闭面边界/全尾与029一般L1/全部共同正尺度沿已接纳链继承。
不把固定H齐次扩到Q，不自行合成I/Cn。

当前只准备代码、源绑定和理论；真实195942新根、97971新斜率/gain与DFS均未执行。
生产全部join并冻结后另绑定actual bundle，主审准确plan及该binding后才执行15worker单批。
正常数学运行无wallcap；全部ready后开gate。
故障先写STOP让worker在行间关闭gzip，grace 0.25秒后仅对仍存活任务terminate，
再共等2秒后kill，再共等2秒及每任务0.01秒有界join；记录真实exitcode或未退出PID。
不把故障清理时限当数学wallcap，不自动重启；被强制结束的未闭合字节只保留，不伪称成功前缀。
