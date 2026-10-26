# A072：Green半凹性的轴族局部定理与rank3交叠诊断

本轮证明了**任意六个独立Cauchy方向所成frame的Green函数，在其真S2外局部半凹**，并给出适用于有限rank3交叠的明确充分条件。没有证明原body32的半凹性，也没有得到body32或其他Cauchy frame的正上曲率发散反例。

固定8个新符号/几何配置及1个保存检查全部通过。它们只涉及一个新axis6解析路径、一个8列交叠几何点和一个抽象函数反例，范围详见下表。F076负责一般frame/body32的理论推广，本包未接管或重复该方向。主流程保留最终证明审定。

## 可以直接引用的新局部结论

[theorem.md](theorem.md) 包含自含证明，分清以下三层。

**三维normal deficit准确身份。** 对任何coercive三维对称Cauchy商密度q_t及符号psi，

    integral_0^infinity t³(q_t(0)-q_t(z))dt
      =3/(8pi²) integral_S2 |theta.z|/psi(theta)^4 dtheta
      =:N_psi(z).

右侧是正范数，故负号后的cusp具有非正distribution Hessian。证明使用非负Fourier差和Tonelli，准确保留3/(8pi²)系数；没有把完整body32局部Green等同于这一profile。

**有限交叠充分条件。** 如果实际函数在一个球上具有准确分解

    f=R-sum_j A_j(x)N_j(P_j(x-x_j)),

其中有限个A_j>=0且C²、N_j为L_j-Lipschitz范数、R有统一上Hessian常数C_R，那么f半凹，常数可取

    C_R+sum_j[sup(N_j)*||Hess A_j||infinity
                       +2L_j||grad A_j||infinity].

它允许strata、normal neighborhoods与相同span任意交叠；不要求disjointness。证明以凸范数的正mollification和准确乘积Hessian式进行。前提中的准确分解、振幅符号、统一余项不能由点态leading profile代替。

**axis6及可逆六列frame。** 任一点off-axis-S2至少有三个非零坐标，记作y，余三坐标为z。准确写成

    p_t=A(y)H_y(t)t³q_t(z),
    A(y)=pi^-3 product y_i^-2,
    H_y(t)=product(1+t²/y_i²)^-1.

由此

    G(y,z)=G(y,0)-A(y)N_axis3(z)+R(y,z),
    R=A(y) integral_0^infinity
               t³(q_t(0)-q_t(z))(1-H_y(t))dt.

在y远离零的每个紧邻域，1-H_y及其至二阶y导数均为O(t²)。小时间对z求至二阶后仍由O(t^(2-|alpha|))支配；大时间利用q_t偶性及O(t^-5)Hessian界，差、梯度、Hessian均给t³之后的O(t^-2)支配。因此R为联合C²并有有界Hessian。上述充分条件即证明axis6半凹性；G_N(x)=|det N|^-1G_axis6(N^-1x)传递到任意可逆六列frame。

若另外有0<m<=G<=M与局部Lipschitz常数L，则G^q、1<q<2的上Hessian常数可取qM^(q-1)C_G+q(q-1)m^(q-2)L²。它说明半凹性路线与连续面积分目标相容，但这里没有body32这些定量常数或任何全I认证。

## 实际固定fixtures与准确结果

|配置|本次实际范围及结果|
|---|---|
|1|Laplace有理身份和三维Fourier归一化，准确得到3/(8pi²)。|
|2|新axis6点(1,2,4,0,0,0)，全部15个pair残差严格正，最小为1；点在选定rank3空间。未使用069旧点。|
|3|axis3法向slice(a,b,0)，a,b>0：2pi²N=(a²+ab+b²)/(a+b)。Hessian为2/(a+b)³乘[[b²,-ab],[-ab,a²]]，半正定，det=0；无参数采样。|
|4|唯一真实axis6 Green路径(1,2,4,s,0,0)，0<s<1/2，以真实product密度时间积分做准确partial fractions，非DD导数有限部。|
|5|1-H_y的准确t²因子及其y、一重y、二重y、混合y导数的t=0商极限，核验统一余项证明所需的阶。|
|6|I6加v7=e1+e4、v8=e2-e4，x*=(1,1,1,0,0,0)：全部28个pair残差最小为1；两选定rank3空间并秩4、交维2。只做几何，0该8列Green/Hessian求值。|
|7|抽象f(s)=-|s|+s²log(1/|s|)：gradient jump=-2，但s>0时f''=2log(1/s)-3趋+infinity。它不是Cauchy或body32反例。|
|8|带正C²权重的平滑范数乘积Hessian身份和范数Hessian的正多项式分子，支持有限交叠充分条件；无网格。|
|9|纯保存数据检查43个已存残差的准确正号、path值/导数分数、交叠秩及各项适用边界；0新几何残差、0新Green求值。|

配置4的准确结果为

    pi^5 G(1,2,4,s,0,0)
      =(|s|+7)/[180(|s|+1)(|s|+2)(|s|+4)], |s|<1/2.

所以G(x0)=7/(1440pi^5)，左右导数分别为+1/(128pi^5)、-1/(128pi^5)，梯度jump=-1/(64pi^5)。右侧二阶极限为217/(11520pi^5)，正但有限。这验证的是负cusp与有限正余项曲率的兼容性；半凹不等于处处凹。

配置6中W1=span(e1,e2,e3)、W2=span(v7,v8,e3)，交为span(e1+e2,e3)。交内的x*却off真S2，严格说明一般有限frame的rank3交叠不能自动归入S2。该几何例不证明body32存在同样交叠，也不判定此8列frame的Green半凹性。

配置7准确反驳的是“负范数leading cusp+C² off singular set便足够”的逻辑推断。负Dirac质量无法抵消附近非零点的无界正经典二阶导数。它强调需要一致C1,1/上Hessian余项，不能被报为目标Cauchy命题的反例。

## 运行、计数与180秒总墙钟

先于任何科学执行登记固定9配置并绑定代码SHA。每个程序由run_guard以新的专属process group启动，科学等待最多170秒，预留10秒终止/回收；总墙钟限制180秒包括清理。超时仅SIGKILL自身组并wait/reap，保存检查点，无自动重跑。本次均正常完成。

|程序|PID|tool chunk|退出|包含启动/回收墙钟|程序内计时|
|---|---:|---|---:|---:|---:|
|8配置producer|77804|1cd201|0|0.502653292秒|0.216608958秒|
|保存checker|77849|57c01a|0|0.024403375秒|0.000817250秒|

实际新工作为43个tiny pair残差、1条真实Green符号路径、1个normal范数符号slice、1个抽象反例与1个交叠几何点。不能把这说成“所有kernel工作为0”；准确口径是body32 fullkernel=0、全fan/64orthants=0、旧baseline=0、旧069 fixture=0、body32 pairs/triples=0、grid=0、MATLAB=0。没有递归agent、安装环境、stage/commit/push、失败或重试。所有固定来源及程序绑定前后不变，子进程均已回收。

069来源manifest为527cb3097e55d346e2867c7f58b0ed66f8a8a590ffe26b07979e9f6e111f8ce8，已同时绑定B069独审PASS。它只提供off真S3的C²，不授予S3上Hessian结论；040正J2发散亦不被误读为实际正上曲率发散。

## 仍待F076/主流程闭合的范围

body32额外耦合列破坏了轴族的精确切向/法向独立分离。必须另证实际多rank3分解、每个同span贡献的非负振幅及其C²控制，并给出穿过所有交叠和相对normal尺度的一致余项上曲率。只有局部leading coefficient或有限个负jump样本不够。

本包没有发现或排除一般Cauchy frame的半凹性反例；唯一实际Green路径属于已解析证明的独立六列族。没有新增body32完整G数据、global g、全I或Cn结论。可引用theorem的第1节normal deficit、第2节有限重叠充分条件、第3节轴/线性像定理，并须同时保留第6节body32缺口。

冻结后的原执行入口为python -B run_guard.py producer和python -B run_guard.py checker，但原目录拒绝重跑。复现需要新登记/输出目录、相同源SHA与fixture，并维持含清理180秒限制。最终manifest绑定必要文件；source路径记录实际原件位置，异机重放应显式重绑相同哈希原件，不能修改本次收据。
