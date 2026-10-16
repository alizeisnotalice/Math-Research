# 第233轮：增长重门与全亏损亚线性来源预算

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮在原far正交通上得到无需终端亏损截断、无终端逆质量、无额外late项的估计，并将原线性来源拥挤费ρ改进为ρ^(2/3)乘对数因子。保任意原输入的全部合法尺度与soft参数，不以固定尺度模型替代原支。

## 核心证明及原费用

原heavy时刻的完整分母满足e_t>ell q c^(-n/2)t^n，ell=1+logn；这不是任意存活子集H的下界。保真实Ω先积分receiver/hard来源，静态列价F≤C0/√wstar，C0=16/3。原权w_t≤(g*/t)^n wstar，与同源Cauchy及增长heavy分母联合，得到节点价

G_jM=(C0/√ell)√ρ c^(3n/4) exp(x/2)(1−exp(−s))，其中ρ=||ν||/(qa^n)，x=nlog(M/M0)，M0=a/c²，s=x−U>0。

同一节点另有QG/Gini价Q_jM=(ρ/4)c^(2n)exp(−x)s。两价取min，保正规节点及唯一封顶点求和，且共同ν,q,a,β跨j时，得到

Xi(ρ,U)=min{ρ exp(−U),18[C0²/(4ell)]^(1/3)ρ^(2/3)[1+log_+(√(ellρ)exp(−3U/2)/(4C0))]^(1/3)}。

原总费用≤2Λ√(K/3)·24e⁴N_n·EΣV||βV||Xi(ρV,U)，N_n≤8+4√nlog(n−1)。维数、菜单、三腿、朝向、比例带全部保留；ρ=0和空时窗独立处理。全U≥0、全终端亏损仅指所选原far交通，未取消其它母量分支。若来源随j改变，保逐j价不能套合并式。

新增空间精化保α_M(y,v)=∫Ω capture·√(P_node(x−y)P_test(x−v))dx≤1。上式growingfloor中的||ν||可替为m_M(y)=∫Zτ α_M(y,v)²ν(dv)；各腿/soft/朝向的对象分别保留，不能免费抽出共同α。另有全carrier固定floor层势，以删除后的正势下降上界连续交通；这不是重领删除负信用。

**仍未付：**EΣ||β||Xi的平均来源预算、精确平方affinity质量的跨尺度/历史总价、√K、near/浅层/终端互补域及输入质量桥。来源亚线性不等于一般根阶闭合，不能通过Young吸收把新维数费隐藏。

## 独审与本机验证

审查者先独立筛查命题，再读作者证明和父稿，核原heavy资格、增长floor、BV跳方向、完整j/M求和与全部前缀。另独立执行作者附件635项及47项检查，计数与主脚本分开。

主脚本141个配置：36个floor费用分布、54个product权heavy时钟（同步全Z/BV与growingfloor）、6个精确有理平方尾阶梯、45个真实物理节点菜单，另6个解析价比较。40/70dps复核；n至10^6，小n用独立Fraction物理尺度循环对照，大菜单保解析正遗漏尾。加性Δ与−logwstar分开，审计发现的初版混用已修正并完整重跑。有限检验不认证FULL、空间来源资格或一般上界；普通浮点不是有向区间证书。

最大有限节点“sum/两价min”约0.312132，growingfloor模型最大比约0.035080；不解释为一般锐常数。n512,U0,ρ=10^12时，新解析价格/旧ρ约0.00417976，仅比较两条已证上界。

## 归档与复现

[research.zip](research.zip)：11成员，85,695字节，SHA256 `1bec82bfc31c3af184aa78f0e7f3d72b85f5987f6fe2d7166bee24bdf6e097c0`。保存详细LaTeX、两份证明、独审、必要脚本与结果；不编译PDF。本轮仅新增本README和ZIP，另更新项目索引。

解压后`python3 reproduce.py --verify-only`核验；完整包外复现：

```sh
PYTHONPATH=/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor python3 reproduce.py --output-dir <新建包外目录>
```

使用D02/N04核对共同归一化与势/误差费用。没有可归属本轮的实际额度扣减数据。下一最小任务是保真实接收位置与共同历史，对平方affinity来源质量建立平均几何收费；不重复已经核定的重门代数，不把新来源矩当常数。

检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round233_20261017_final_checkpoint`；ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round233_20261017_research.zip`。工作文件、检查点、ZIP及备份已逐成员核验；最终包外重放回执另行记录。

冻结包外完整重放141配置成功，完整JSON逐字节一致，10项manifest运行前后通过，约35.219秒；Python3.14/mpmath1.3.0。回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round233_20261017_final_replay/receipt.json`。逐成员备份核验及重放完成后才移除本轮工作目录。提交前项目根821项、rounds目录232项，Git pack约366.94 MiB；本轮新增ZIP仅85,695字节。
