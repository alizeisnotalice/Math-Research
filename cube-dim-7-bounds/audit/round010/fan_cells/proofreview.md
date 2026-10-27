# R010 ordered cells/pulling独立Sol数学与结构审计

2026-10-27。新独立clone起点 `75eedb49b0e95dff748081d658d5a9f53a593974`，仅写 `audit/round010/fan_cells/`。全文读Sol cells_proof.md/registration和Terra generator、fullchecker、report、receipt及输入/输出数据schema。本审不执行rank、determinant、cell/pulling或R009近两百万cofactor任务，仅读源码和有限证据接口，最后仅核SHA256元数据。

**结论：接受partial-cell到最终ordered-frequency triangulation的证明与实现规格；未发现实质阻断。** 全compatible-ray convex hull、rank/full-dimension、29层split、最终顶点、36零集facets及positive pulling的覆盖/内部不交链均成立。它们与已经完整认证的R009 ray catalog及本轮完整保留/检查记录相接后，足以建立完整ordered frequency triangulation；R009单独ray list原本不够。

实际135 cells、212 top simplices、629遇到的face-DAG nodes归属Terra实际构造/fullchecker及主集成人另行portable运行的证据，不是本审重算结果。精确volume1是独立算术校验，不单独承担覆盖证明。没有H、g、I或Cn数值评价。主集成人保留最终采用和证明责任。

## 1. 前提、截面与中间vertex全包含

输入是R009完整46个primitive nonnegative w rays、七coordinate walls、29cut-prefix normals及64body forms。完整性必须来自R009 active-rank证明和实际全量认证；本轮复制hash或每ray witness只核validity，不能单独重新证明没有遗漏。R009结果已经另审，本轮没有重新扫描1,947,792子集。

坐标ξ=Mw、detM1、lambda=Σk k w_k=Σξ_i>0。对任意partial sign prefixσ，Pσ由w≥0、lambda1和已施加cut弱不等式定义，故是bounded polytope。任何vertex处的active homogeneous rows加lambda须rank7，而它们annihilate非零ray使homogeneous rank≤6，所以可选六个独立active rows。它们都属于R009的36normal全集，该vertex direction必在完整catalog。

另一方面每个compatible catalog point w_i/N_i本身满足Pσ。有限bounded polytope等于其vertices convex hull，因此

`Pσ=conv{w_i/N_i : i compatible withσ}`。

这一论证允许中间compatible点冗余，绝不假设它们全已是partial vertices。对lower-dimensional非空partial cell仍成立；empty set按rank0处理。

## 2. rank7判full6D及29split完整覆盖

所有归一columns落lambda1，正scale不改linear rank。其affine differences在ker lambda，而任一basepoint的lambda1，因此affine dimension=linear rank−1。raw primitive columns rank7恰好对应section full6D，准确且无需floating LP。

下一cut将parent compatible bitset按dot≥0及dot≤0分两支。每支的true vertices也属于global catalog，因此各child conv恰好是parent与相应halfspace交集，不是一个只筛旧vertices而漏掉intersection vertices的通常错误。

对full-dimensional parent，若cut穿interior，两child都full-dimensional，interiors按严格相反sign互不交。若不穿，唯一低维child在cut超平面上，是retained child的boundary subset；舍弃它不丢closed union。至少一支full-dimensional，否则parent不可能full-dimensional。归纳29次给闭simplex全覆盖；不同final cells至少一个cut在其relative interiors有严格相反sign，故内部不交，只交lower-dimensional boundaries。

同一rank7 bitset不能有conflicting signs：冲突cut必须在所有该columns上为0，与rank7和normal非零矛盾。source显式reject冲突，而非掩盖duplicate。所有29parent/child masks、ranks和sign paths保留，足以由portable checker逐项检查这一归纳链的实际输入。

## 3. 最终vertices与整数interior witness

最终cell已包含全部29cuts。任何compatible catalog ray都有六个独立active-row witness；这些rows在该final cell都是有效weak inequalities，在此ray为0。加lambda即rank7，所以它是该cell section的true vertex。故final compatible set恰好是其extreme vertices，中间冗余问题已经消失。

整数witness=Σ_(cell vertices) primitive w_i。对任一wall或signed cut，其值在所有vertices上非负；rank7保证该非零normal不可能在所有vertices上为0，所以至少一项严格正，sum也严格正。所有walls和cuts同时严格，使witness位于真正interior。无需通过浮点样本猜cell sign。

## 4. 36个zero sets/rank的facet完备性

对final cell的任一face F，原36normals各自在整个cell有固定weak sign。F_n=F∩{n=0}因此是exposed subface。r=rank(F)时，proper、非empty且rank r−1的F_n正是cone facet；identically zero的normal给F自身，必须排除；相同facet bitsets应dedup。

无遗漏：取F的一个proper facet G的relative-interior点。F在其linear span内仍由原cell inequalities界定；若该点没有新original inequality变equality，它就会是F的relative interior点。故有原n非identicallyzero于F但在G为0，其zero face为proper并包含G；facet maximality迫使zero face=G。

递归face始终是原cell与若干原zero hyperplanes之交，vertex set恰为落在交集上的cell vertices，没有凭空新增vertices。bitset决定其cone，rank决定section维度；同一bitset在不同cell表示同一geometry，global memo合法。冲突cache必须reject。

source记录629个**实际遇到/递归所需**face nodes，并未保证完整face lattice都成为DAG node。一个node存下全部facet bitsets，但pulling只访问不含anchor的facets；未访问的其他facet不必都有自身node。这已经足够验证本算法的递归覆盖，不应把629解读为整个polyhedral complex的全部faces数量。

## 5. Positive pulling的准确覆盖、内部不交和共享面

rank1的face只能含一个primitive catalog ray，否则存在duplicate collinear ray。一般face F选global-order最小vertex a，找所有不含a的facets G，对Pull(G)各simplex加a。G的span有rank r−1，a在其support hyperplane之外，所以每个r-column输出simplex独立，且无Steiner ray。

在compact lambda-section，从a通过generic interior point延长的ray唯一退出到某个不含a的facet；含a的facet不可能在正延长方向首次阻挡一个interior point。不同generic exit facets只在lower-dimensional边界相交。对facet递归后，其simplices也覆盖且内部不交，coning保持这一性质。boundary由closedness/limit覆盖。因此不靠volume等式，就已证明整个F的positive triangulation。

共享face兼容也可按维度归纳。若face E不含a，它包含在某个不含a的parent facet中；coning simplex与E的交只来自该facet端，递归用相同global order给Pull(E)。若E含a，限制变成a对E自身不含a facets的cones；这些facet均由parent facets截出，递归限制仍用同order，所以也给Pull(E)。这证明ordered cone内相邻cells的conforming restrictions，不只是每cell单独覆盖。

不同final cells不可能共同含一个independent7-ray top simplex：相反sign的cut须在这七columns上全为0，与其spanR7矛盾。source将此类duplicate当error，不靠global volume抵消它。

## 6. volume1校验与primitive kernel尺度

U=diag(1,…,7)W diag(1/N_i)的columns sum1。它将standard affine simplex映到目标simplex，sum-functional保留，所以完整det在transverse quotient的factor为1，|detU|等于六维section volume ratio。故

`vol_ratio=7!|detW|/ΠN_i`。

这里standard section的Euclidean sqrt7/6!因子在ratio中消掉，不应再乘一次。全coverage/内部不交已由前面证明，才推出Σvol1；实际准确Fraction总和1检验det、N和normalization。单靠Σvol1仍可能掩盖重叠和缺口，故source/checker另外验证29层及pulling完整性是必要的。

cell exponent `ell_C=Σ_noncut prefix +Σ σ_j cut_j`正确。35noncut forms在整个ordered cone≥0，cut signs在cell固定，因而每ray `q_i=ell_C·w_i=Σ_(64body)|v·M w_i|>0`。等式对所有positive combinations成立；没有漏weight、单位化或128条有向body的额外factor2。

Rξ=MW、detM1给 `alpha=|detW|/Πq_i`，`b_i=x·M w_i/q_i`。primitive ray和section-normalised ray的det/q_i scale同步相消，b_i不变。若使用u columns，频率变换有1/7!；detU若配primitive q_i会错误混合N_i和7!。Terra存primitive W、N、q、alpha及volume，checker逐项使用准确integer/rational公式，约定一致。

## 7. B7全频域覆盖与DD6接口

fundamental ordered cone的七wall非平凡；任何independent7-ray simplex的interior都是strict w_i>0，从而ξ1>…>ξ7>0。B7 signed-permutation images覆盖R7，full-dimensional interiors互不交，sorting walls为零测。其group action preserves原ψ和absolute determinant，所以全frequency positive cone decomposition来自对212 ordered top simplices使用实际group phasesgξ_i。

ordered cone内common pulling order证明conforming triangulation。跨所有B7 sorting walls的全局face-to-face refinement未另作证书；Fourier积分仅需已经证明的闭覆盖、内部不交和零测边界，因此不需把此额外严格fan术语当公式前提。不能对一般x将单cone项乘645120，也不能按ray stabilizer除top-cone multiplicity。opposite pairing只在真one-per-pair index set上乘2Re，全部images已保留时不再乘2。

于是结构证书足够定义原process的exact finite p_t cone product，以及R009已审positive DD6 H_epsilon表达式。q_i>0、epsilon>0使cone/time/Fourier absolute integral `∫e^-epsilonψ/ψ`有限；log用right half-plane，重复/零knots用confluent derivatives，all-zero DD6=1/(6epsilon^6)。这个桥接不改变原64weight-one共同时间模型。

实际group images/phases和signed sums本批未enumerate或evaluate；可被数学上完整索引，不等于已经给出H数值certificate。后续计算仍须controlled cancellation、准确repeated-knot处理及closed face continuous enclosure。

## 8. 读取源码与实际135/212/629证据归属

本审只读run.py/check.py：all-catalog mask、29 sequential splits、rank7 child、完整zero-facet tests、global-orderpull、shared DAG一致性、topduplicate禁止、interior/ell/q_i/det/alpha/volume全部与证明接口相符。source中的assert不是本审执行证据；没有运行Terra fullchecker或任一rank/det job。

Terra receipt/report实际记载29层、135fullD cells、212top simplices、629encounteredfaces及volume1，并报告已完整检查保存的split/facet/pulling records；主集成人另行portable执行可把运行PASS与本报告的source hashes绑定。135/212/629不是由本审源码阅读推测或新计算得出，也不能用Sol作者receipt中0数值运行替它们背书。

该fullchecker除了volume还：验证R009复制ray witness及hash；重算每个level parent/child rank；逐encountered face检验全部36zero intersections；递归重建Pull并对比每cell top tuple；核每个simplex的signs/det/q_i/alpha/volumes和global uniqueness。因此它所检查的数据，配合本报告的数学完整性链，足以成为完整ordered-frequency triangulation证书。无需重跑R0092m，亦不需本审额外构造另一套cell generator。

本审已锁定全文、generator/checker和actual retained evidence hashes，核theory及numeric manifests、pre-run inputs。SHA元数据核验与科学rank/determinant复核严格区分；仅前者由本审实际运行。

## 9. 原G transfer与结论边界

与R009一样，H_epsilon≤G，连续full-face H下界可通过G的−6homogeneity给cube essential g；fixed-epsilon H本身不被误当homogeneous。原whole-space positive-part delta-source、一般L1/全部共同r>0比较仍用于G，不能用LH_epsilon=−p_epsilon替代它。

若未来有H^(7/6)全face integral上端，可用R001 whole-face short-time预算支付G−H，或复用独立full-G I预算。结构输出没有给出这些连续算术，也没有H、g、I或Cn evaluated bound。

因此本轮完成的提升是：R009 ray catalog在本轮split/facet/pulling证书与实际checker证据支持下，可成为完整ordered-frequency triangulation以及足够Fourier使用的全frequency cone decomposition。它不是新Green值或算子上界。

本审到此停止，无数学数值脚本、MC、MATLAB、子代理、commit/push，旧封存及其他工作者文件未改；不自动启动后续任务。主集成人保留最终证明和采用责任。
