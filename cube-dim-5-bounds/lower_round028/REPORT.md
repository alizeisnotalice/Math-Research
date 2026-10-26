# 第028轮：period3位置机制的非局部静态输入与有限事件准备

本轮完成整个 `0≤epsilon<a1` 窗口所需的**全213 tail-shape静态数据**，准确接续025本地fine prefix，成功旧cell零重算。全数据为139843个shape masks，精确压缩为91个首轴descriptor；首轴279个仿射端点产生965个准确去重的闭窗口安全事件、964个开cell。

首个安全sectionwise swept-interval上放松为 **约2.2260292432086652**，严格高于完整C，差约 **0.23171278164017930**。该放松太粗，不能排除整个非局部窗口，不能当成一个可实现见证比，也没有据此选出唯一contact或候选。本轮**实际epsilon有限functional/union评估数为0**，没有参数网格、964-cell真实union扫描、有限epsilon新候选或Cn宣告。

完整事件/PWA身份可把旧025准确负导延伸到首个全局cell：首正事件约为 `0.00405257139915791107`，该区间内 `R(epsilon)=C+epsilon*D025`，对 `0<epsilon≤该首contact` 严格低于C。这是继承的仿射排除证书，不是新基线或新epsilon union实验。其余963个开cell未作实际参数union评估。

## 全窗输入、原子顺序、distinct与safe starts

固定原C格 `a=a1=1+u³/2`、四tail间距 `s=(1+u)/2`，`u⁵=2,v⁵=3`。原子质量全为1，首轴位置为

```
p_j=j*a+epsilon*z_(j mod3), z=(1,0,-1).
```

一个周期的三条gap准确为 `a−epsilon,a−epsilon,a+2epsilon`，所以epsilon<a时原子严格有序且不同，`p_(j+3)=p_j+3a`，周期质量恒为3。共同半径 `r=N^(1/5)/2`、N=Πell，与epsilon无关，五轴同r。

全窗j=−1,...,4的资格不能照搬局部“最大位移很小”估计。这里：j≤−2时，逐相位可见 `p_j≤−2a`，而r≤2<2a，因此首轴receiver上端在0左侧。j≥5时，block最后索引k≥ell1+4，且 `p_k>(k−1)a≥(ell1+3)a`；又 **r<ell1*a**，所以receiver下端大于3a。后一个不等式对ell1=1由 `r≤u³<a=1+u³/2`（u³<2）得到，对ell1≥2由r≤2<ell1*a得到。全部213个正tail形状均准确核对这些不等式。

因此全部1024长度×6 starts=6144 raw覆盖整个窗口。tail固定为空/零宽的811个长度给4866 raw；原213正tail形状完整保留在静态数据中，没有沿用旧51形状局部过滤。

ell1=4可另作**全窗**排除：first/last原子相隔3a，因为block跨一个完整period；而2r≤4<3a。213项中的90项（540 raw）因此始终首轴为空。余123项、738 raw进入全窗capture事件；长度3等在较大epsilon确可出生。例如ell=(3,1,1,1,1)、j=0，其首轴宽度为 `2*(v/2)-2a+2epsilon`，在 `epsilon=a-v/2` 从零打开；tail盒为正。这说明旧局部strict-negative排除不可用于本窗口。

epsilon=a时两条gap为0，三个相位原子碰撞。本轮只把这个端点用于有限receiver公式的**连续upper闭包**，不能按“3个互异单位原子”直接登记新witness。所有合法候选必须仍满足epsilon<a，并重新做共同半径严格收缩与collar/L1核验。

## 实际旧prefix与零旧cell重算

先读并核对真实输入，没有凭口述数字猜路径。旧prefix为

```
/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-continue/
cube-dim-5-bounds/lower_round025/tail_checkpoint.json
SHA256 c92eb7cbb1555ddae5f1ceb95f880d0bf63f794084de20b542b8fc81bc08664b
```

真实状态为11/26 first-tail slabs、193336 cells、58755 masks，213正tail形状，每轴26段。几何文件SHA为 `d92de1328b64901e0685a08b67bea9b745401d81f37d26f8d23c4c2409f23bf7`。旧源和prefix均未改。

新计算只写028与外部 `/Users/zhengzhihao/Desktop/T/cube5-lower-scratch028/`。全histogram放在本地SQLite，逐x-slab把该slab的权重与completed-prefix放在同一个事务中提交；中断只能留下前一个准确slab。使用直接Q(u,v)系数卷积、稀疏累加和小型SQL缓存，避免把全部细mask塞入旧无界F乘法缓存。

余下263640 cells分两次有限进程：先新增193336 cells推进到22/26，再新增70304 cells完成26/26。每进程均小于200000新cell，**没有重算任何旧193336 fine cells，也没有重算025的51形状coarse数据**。最终456976-cell全213 histogram有139843个不同mask，四tail物理总量准确为 `((1+u)/4)^4=(s/2)^4`。

必要完整数据以 `full_tail_histogram.json.gz` 流式导出，16839842字节，SHA256：

```
394d43e3bef7e8bb0aeed25e5d7f6ee56d5e41ee36c0d697421f3a33a134636b
```

原大fine prefix仍本地保留；新SQLite是scratch恢复状态，不归档为仓库产物。压缩完整histogram与流式reader足以使用、审查全部数学数据。

## 新213输入的独立pair-factorization核验

没有复用025的51-shape核验count。先实际构造两份全213的2D signature数据，各650项，准确新AND对象为650×650=422500。另行登记后分三批199550、199550、23400对，均小于200000；每left-signature row事务保存精确prefix。

完整pair-factorization输出与流式读取的全139843个fine-slab mask**每个精确系数完全相同**，总物理tail量也相等。它和fine-slab producer使用不同的几何分组、mask交集和求和路径，但共享Fraction与Q(u,v)卷积内核，不声称独立底层算术。没有新epsilon/5D参数union。

## 91个准确descriptor与全窗PWA接口

每个tail section保存全部shape身份。相同ell1的首轴接收区间只依赖first/last原子位置和半径；较大N给较大r，区间对每一个j及epsilon都包含小r区间。因此在一个固定section内，只须保留ell1=1、2、3各自出现的最大N。被选择的family本身仍在该section中，所以没有拼接不同tail盒或混用半径。

完整139843 masks可准确压成91个 `(maxN_ell1=1,maxN_ell1=2,maxN_ell1=3)` descriptor，并把其**互不相交tail section**权重聚合；这是准确等式，不是upper relaxation。权重总和仍为 `(s/2)^4`。每descriptor只剩最多18个首轴raw区间（3长度×6 starts）。

在有序窗口中区间为

```
[max(0,(j+ell1-1)*a-r+epsilon*z_last),
 min(3a,j*a+r+epsilon*z_first)].
```

端点准确仿射，tail权重固定，所以名义有限functional是epsilon的有限连续PWA函数。91-descriptor数据已为后续局部interval upper或单个contact认证提供接口，但本轮没有逐965事件进行实际R值扫描。

## 完整安全事件inventory

收集738个全窗保留raw的首末capture端点替代式，再加入固定0、3a域面。准确去重得到279个仿射式，斜率仅−1、0、1。不同斜率的交点为 `(d-c)/(m-n)`，仅除以整数，仍在Q(u,v)中；包括裁切、空盒出生、不同case接触及端点交换。

先登记实际279式、38781对，然后准确比较并按完整系数去重排序，得到**965个闭窗口安全arrangement事件、964个开cell**。每个事件保存精确系数、来源form-pair和有理隔离区间。该集合是完整安全超集，可能包含在实际union中不活跃的交换；本轮没有声称965点都是最终R的有效拐点。实际R只能在这些点中变换仿射式。

全部6144 raw仍在 `capture_forms.json` 逐项会计，ell1=4与固定tail零体积仅按各自全窗证明排除，没有套用局部51名单。没有旧1946事件重扫，也没有epsilon密集网格。

## 首个安全max-relaxation与真实缺口

对每descriptor内每个raw区间，分别取整个闭 `[0,a]` 上的最小lower端点和最大upper端点：

```
L_min=max(0,min(c,c+m*a)),
U_max=min(3a,max(d,d+n*a)).
```

该swept区间包含每一个epsilon的实际raw区间。先并这些swept区间，再对tail descriptor权重求和，得到

```
R(epsilon) <= (16/3)*Σ_descriptor weight*length(union swept intervals).
```

这里不同区间、不同receiver点可独立“选epsilon”，所以只是安全放松，并非某个可实现参数或见证比；不存在按各family面积相加当并集。每descriptor内仍作区间union，tail权重来自不交section。

91项上包准确为约2.2260292432086652，相对完整C的残差有理隔离**下端严格正**，约0.23171278164017930。故这个全局放松不足以排整个窗口；不能从upper高于C推断存在获胜参数，也不能据此挑出唯一未排contact。更细的964×91 interval-bound对象尚未启动。

## 首cell的继承排除，不是新参数实验

完整全局端点arrangement的首个正事件为

```
epsilon_1≈0.00405257139915791107234190609501.
```

在0与epsilon_1之间，没有任何raw端点排序/裁切/出生事件，因此R准确仿射。原C身份和025的D025≈−0.11835180380160859给出 `R(epsilon)=C+D025*epsilon`。有限interval union连续，把身份延伸到首contact；其 `R(epsilon_1)-C` 准确小于0，展示约−0.00047962913512514751。

`first_cell_certificate.json` 保存该身份、准确首事件与负残差。这里没有重算旧C或在epsilon_1重新运行union；只是从完整PWA资格继承已有C/D。它不排除其余963个开cell。

## 实际运行、故障及复核

所有数学进程nice10、单线程，六个数值线程环境变量均为1，无wallcap；所有批次对象有限、可恢复。各核心进程如下：

| 对象 | PID | 秒 | 峰值字节 |
|---|---:|---:|---:|
| 新fine cells 193336 | 79627 | 20.2370 | 216940544 |
| 新fine cells 70304 | 79989 | 8.5030 | 91209728 |
| 全histogram检查/压缩导出 | 80742 | 5.7562 | 77938688 |
| 全窗输入/descriptor/pair count plan | 80969 | 5.0639 | 58163200 |
| 38781事件比较 | 81846 | 1.5767 | 145096704 |
| 91项安全上放松 | 81892 | 0.1582 | 45170688 |
| independent pair 199550 | 82265 | 21.0691 | 52674560 |
| independent pair 199550 | 82629 | 25.8918 | 54755328 |
| independent pair 23400 | 83338 | 4.0025 | 54870016 |
| 流式139843系数完整比较 | 83909 | 5.3515 | 54362112 |
| 准备材料最小复核 | 84056 | 5.7637 | 132874240 |

有一次真实代码解析故障：原 `prepare_window.py` 的keyword误写成冒号，SyntaxError发生在任何plan数学计算之前，没有损坏static prefix或产物。保留失败source、原登记和 `plan_fault.json`；另写修正v2并在执行前重新登记，没有改写旧登记或盲重启。原 `.py` 本地不改，归档同字节/同SHA的 `prepare_window.failed.txt` 作为锁定历史输入，避免把不可执行的失败版本混入可执行源码；如需历史重放可从该文本恢复原文件名，明确不应执行。

最小复核再次准确检查38781事件对的完整性，从压缩artifact流式读取全部139843 masks，独立重建91个descriptor权重并核对总量、上包-C正残差与首cell继承证书；旧025/026全白名单哈希未变。它没有实际epsilon functional评估。

现有证书审核命令：

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round028/verify_preparation.py
```

若需要重放独立全static核验，`verify_pairs.py resume`的实际对象已登记为422500对，每进程≤200000、按left row精确恢复，最后 `verify_pairs.py compare` 流式逐系数核对。当前prefix已完整，日常审核不必重算。

## 交付与边界

只写028及外scratch028；未改025/026或其他研究者文件，无MATLAB、提交或推送。完整压缩静态数据、准确capture/events、91-descriptor接口、安全upper、首cell身份、source/故障/批次receipts及最小审计列入白名单；SQLite与重复的完整event checkpoint仅作本地恢复数据。

任意后续候选仍须epsilon<a、互异原子与正gap、共同r严格内缩及固定collar/等质量平滑L1传递；epsilon=a只能连续upper。当前不是完整窗口排除，没有选出唯一新候选，也没有一般输入或Cn上界。主代理保留proof、L1/collar终审。
