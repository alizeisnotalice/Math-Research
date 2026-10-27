# R064正式源的局部接口接入锁

**Local interface ACCEPT。** 本附记只对照正式源对已冻结相对坐标引理的调用，不审核target径向列、target计数、actual-history能量或完整C_n结论。

正式源：`/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round064/radial_calibration/PROOF.md`，实际SHA-256为`17ce38ab5af46c0358046bdcc2aa57a05844c0799e9a97280d959887f59615d6`。

本次实际读取的§3正文段SHA为`c91acca59ecae204291f873a1854e82d19b07e1a1f9faa9f972f2bf6f64f97ea`；§6正文至“Before cap restriction and before soft weighting”之前的phi/face-mixing段SHA为`43e5ab8850ca4c2a5ff2b0e26c94037895d37de6cf13cae7f9af168c959c4584`。这两段身份供后续局部增量对照，不扩大审核责任。

实际全文读取局部相关§2–3、§5乘积桥及§6至face混合，确认以下接入：

1. s=3/2，source q∈[2/3,909/1000]，target r∈[1+δ,sq_max−δ]且δ=2^(−30)。h<δ/2使完整匹配区间[(r−h)/s,(r+h)/s]位于source band，active R=sq−h≥r−2h>1、r<3/2。inactive fullband点可能R<1，但源明确其columns为零、不调用本局部矩引理。
2. ε=ζt及εq_max≤h/s、h≤1/1000保持a_0>0、g=R−a_0=(s−1+ε)q−h>0，且g−R/3=εq−2h/3≤0。+ζρ项保留，未用旧absolute d域。ε≤h<1/2保证g'=1/2+ε<1。
3. baseline q_b=(r+h)/s是右匹配端点，确有R=r。没有把r/s误用为mean-one baseline。q区间宽度为2h/s=4h/3；x=909L/1000给n/ρ≤2727L/2000<(3/2)L。因而256导数上界与AC积分给512Lh=2^(−55)，无需K_w条件。源以更宽的2^(−50)支付该shift，方向安全。
4. 源重写了完整新链R'=3/2、g'=1/2+ε，F_u≤1/5、|F_R|≤1/15及F≤R/30。其B'、B/q与|∂logJ|<3/ρ推导符合冻结引理；onset kink按a.e./AC处理。个体|J−1|<1/(3ρ)一致推出最终固定eps_J正常数的归一化界。
5. 源baseline明确W是target Uniform[-r,r]，满足R=r、R≥1、g≤R/3≤R/2、ρ≥4R，故EJ=1及3/6/12合法。pinned target坐标为确定项[−g/ρ,0]，没有另赋Uniform law。
6. 源phi矩明确用FULLSOURCE Z~Unif[-ρ,ρ]，两个精确式、capture inverse密度及J_i≥1−eps_J后的去indicator上界正确。相对条件给Q_pin≤log(6/5)，R/q≤s；Q_5=34187/187500与本次新Fraction核验一致，mean rate<999/20000、variance rate<1/150。
7. 每个captured signed source face内才使用独立free坐标。源分别以共同均值上界1/5+(999/20000)x与variance≤x/150应用Chebyshev，再混合概率；没有把mixture variance写成独立variance和。χ−999/20000=1/20000，x≥8000使1/5≤x/40000，故其cutoff gap≥x/40000及40000²/(150x)概率上界正确。
8. 逆概率估计先在fullsource、固定face、capture-and-inverse law中完成；源随后另付原cap删除，未假设cap-conditioned坐标仍独立。ell_pin的q导数在源中正确为负。后续利用该负号证明的count calibration不是本附记的审核责任。

冻结依赖实际保持：

* `LEMMA.md`：`c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`。
* `verify.py`：`440980a84518ae4154619abb818fb19432f2362ff471978afbecf81933aaaa95`。
* `verification.json`：`a01b9cf2774cd287abedc1e2f5622771c5ad1283c80e02ba3475df30c1fd4742`。

它们与full作者首次真实read并绑定的身份一致，本附记没有修改其中任何字节。此次只做源读取、身份核查与解析调用范围对照，没有重新运行SymPy、Fraction检查、作者certificate或Terra数值批。新局部数学与实际新计算证据仍完整保留在原三个冻结文件，不能将本附记误作另一场计算run。
