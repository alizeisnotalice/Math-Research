# R070 LOCAL对R071的接入绑定

**局部调用一致。** 实际读取R071源：

`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round071/second_profile_cantelli/PROOF.md`

SHA-256：`0e36488654bfa1413c19335a3eb9cf7b900e16b25cf9b73b2b8f36efc7381a10`。

本文件只确认新H2局部接口在该源中的实际适用，不签完整Cantelli预算、计数、径向列、actual-history或C_n全链。

## 1. 参数与closed/baseline接入

源s=3/2、q∈[2/3,.91]、r∈[1+2^(−19),1.365−2^(−19)]、x=.91L、L≥8192、h=1/(256L)。它们满足R070的r<3/2、h≤2^(−20)及x∈[L/2,L]。源明确q_b=(r+h)/s、匹配区间[(r−h)/s,(r+h)/s]、active R≥1、ρ≥4r、εq_max≤2h/3，以及correct +ζρ gap。inactive R<1处不使用局部矩。

源§3的J2正是R070定义的capture phi²/ell积分，既不是J1平方，也不是逆权V²。闭式R/(3ρ)(Q_g³+Q_v³)、m>g时中央零profile、Fubini target均值均与R070一致。

R070 target mean是**完整n−1个free target坐标加一个确定pin的总均值**，不是n条free均值。n≥2000及h≤2^(−20)时已证其在[.00528,.005295]x；源使用更宽[.005275,.005295]，方向安全。源的额外pinned上包1/(Nq_min)较R070的1/(250Nq_min)保守；其明确保留(1−1/N) free比例，未丢失一项。

## 2. 有限h transport及target概率

R070新强导数<.1/ρ支持源§3的弱1/ρ，并按a.e./AC使用。源h在8192基点为2^(−21)，所以weak bridge 2xh≤10^(−6)x也成立；R070更强bridge xh/5无需源改写即可支撑同一matching区间。

源baseline事件ΣJ2∈[.00527,.0053]x经过transport进入[.005269,.005301]x，与LOCAL完全一致。target law是full signed target face、free coordinates独立、pin确定，未当作actual μ模型。R070的显式坏概率5,760,000/n已可支付该事件；源使用更宽的点界与自己的N(L)=L^64/64! gates也是合法上包，完整43-gate表不在本接入文件重跑。

源一阶profile采用其已接受的strong <.5/ρ及shift≤xh。在此wide-h参数下h<10^(−5)，baseline .04421x到matching .04422x预算仍成立。它不是将旧tiny h阈值无核对搬来。

## 3. Source-face variance接口的准确方向

源§5在每个captured signed source face的full inverse law下使用free means=J1_i/J_i、raw second moments=J2_i/J_i，pin为确定Q_pin。这与LOCAL概率空间相同，cap删除在其后另付。

R071只需要upper variance：

Var(V|face j)≤Σ_(i≠j)J2_i/J_i
≤Σ_all J2_i/(1−eJ)
≤.005301x/(1−2^(−20))<.00531x。

这里丢掉negative mean-square项、把非负source-pin profile加回去，均为正确的upper方向。R071没有调用R070的lower variance，也不需要在这个upper步骤显式减均值平方；不能因此误说它声称了variance等于raw second sum。

R070更强lower接口则已严格去任意source pin并减free均值平方，其附加guards n≥2000、x/n≤.001、J误差≤1/1024也由此源参数满足。但该lower接口不是R071全链的依赖，不应与其他BE候选混为一项。

源把条件方差上包用于每face分别作概率估计，之后混合概率；没有用mixture variance等于条件variance和，也没有在D-cap条件化后保留独立性。Cantelli的具体retention预算及完整operator结论留给R071整链审计。

## 4. 冻结身份与范围

实际核对R070三件文件在首次发布后保持字节不变：

* `PROOF.md`：`4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8`。
* `verify.py`：`0f28087a2600ac84d7dfd0f1cf8b35439097be4251e0e8ab8e7ff2e5c08857ee`。
* `verification.json`：`8d1e1d722b2a773bda7a42db4be259ea1c0d8b8dbf1b8448c5016273742be97c`。

这次是源读取及应用范围/身份核对，不是新计算run；未重跑R070新程序、作者43 gates、旧baseline或生产grid。R072/R073等其他tuple不能仅凭本文件继承R071完整阈值，须核自身guard；R070 §7只是已明确重推的LOCAL finite范围。
