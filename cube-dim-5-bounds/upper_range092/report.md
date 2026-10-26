# upper_range092：新chain ratio完整零H重range/mean/chord

**完整56139节点、14个有限程序全部完成并通过独立读回；零H、零新anchor、零DD/stencil、零old exp/asinh重算。** 091是新数学ratio机制，本轮没有重复旧点核基线。

最终真实完整正面F上的截时积分证书可简写为

\[
 26069759/1000000000000\le\int_F H_{1/512}^{5/4}\le885871/15625000000,\qquad
 \boxed{I=\int_FG^{5/4}\le6010267/100000000000.}
\]

内部细上端UH≈5.669574350580058e−5、I≈6.01026698256809e−5。**I上包严格优于088完整cap70301737/10^12；配合已合格g=1/10000，cost上端≈44.93722703599，准确32门未过，正式C5仍32。** g105未消费。LH只标记截时H积分下证，不将upper-budget表达式的下端称作实际I或C5下界。

## 1. 冻结新机制和输入

源091 manifest6589e35086bd5082c66c4e94219b0f7c80c393e735f5011eb6b2a6d24d22665b、ratio2931e5d34b71d93eb11dcdb52d5ea0921f9fdbfa0e5e814a60f4ff64eed0de3d、report0989e72688f0235aa5086cf8844ae8cb77facd8ebaf5d71ac7c1323b79dc81e8，全部source文件实际核hash。源088 manifestbd861990f4e623bda5697a6e266b3c6f8b8e61b0b7bb5844d9aff7f59bc5d39d保持字节锁定，复制初始state SHAe764f8ce6fd4f8f1c98813dd8afb8974e069d8fbf18c81100255c410b3acecf5。

同一body16未单位化方向、alpha1、逆Fourier(2π)^−5、eps=1/512、q=5/4、同anchor/真实box/73端点全部保留。085已认证的同几何M2=E||u||∞²、μ和088旧mean直接复用，不重复旧moment基线。source012已接受g100的小快照SHA0b5251bbad3fb507e1965d3ee7eb3fbab6526f7019880aed34a907a9d5c6beae；g105不属于本输入。

088没有逐band保存oldRU，因此按主明确补正：**每band仅用091新RU，之后与保存的旧aggregate range/mean/self积分做交集；不重新计算任何old exp/asinh。** 这是机制登记的一部分，不依赖新RU的192-bit向上端在数值上小于某个未保存oldRU。

对真实h∞和同一band左端a，x=h∞/a：

\[
R_{new}(x)=\left(\sqrt{1+x^2/16}+x/4\right)^2
\left(\sqrt{1+x^2/256}+x/16\right)^8.
\]

091全域face/segment/positive-time证明支付其使用，无144、D、sqrt5或16时间尺度改动。全部节点/73端点实际只有89个不同x；ratio_registration先冻结有限89输入，再由091 bits192产生catalog的整数根与幂向外端。不是从89案例推断全域定理，定理来自091已锁解析证明。

prepare逐原raw文件核SHA，只把原73端点及point复制到外部endpoint-only packet，不再让production/readback读取无关stencil。节点表冻结每个packet SHA及original raw SHA。预处理实测60.2643秒，主要是旧raw IO；没有生成H或weight。

## 2. 新range、mean、凸弦与严格交集

每band原质量[mL,mU]准确由保存H端点相减（H_1支付无穷尾），全部严格正。对新RU：

  Lnew=ΣmL/RU，Unew=ΣmU RU；
  Mnew=ΣmU min(RU,1+||μ||∞/a+(25/48)M2 RU/a²)。

保存旧088 aggregate均有效，故rangeL=max(oldL,Lnew)向下，rangeU=min(oldU,Unew,3)向上，mean=min(oldM,Mnew,rangeU)向上，随后断言0<L≤M≤U。全部量化到2^−160，错误不会靠强行调整M或吞异常通过。

z^(5/4)在[L,U]上的递增凸弦支付真实mean，给积分上端。volume·L^(5/4)给合法积分下端。selfIU=min(oldselfIU,weighted newchord upper)，selfIL=max(oldselfIL,weighted newlower power)。本阶段允许lower提高的唯一理由是新 **真实全box range lower**，没有从mean或upper recipe制造lower。56138个upper与全部56139个lower严格改善；1个upper保留旧证据，仍完整保存。

## 3. 完整F与层级账本

全部geometry、anchor、320条parent展开边、55819个physical leaves均不变。每个节点的self是完整orbit的weighted积分，weights/体积与082相同。数值层级继续effU=min(selfU,Σchild effU)、effL=max(selfL,Σchild effL)，同一parent full orbit上的合法交集。未消费节点付旧088 self，故每个4096前缀仍是完整F证书；UH不增、LH不降，体积始终1。

所有14个prefix及最终真实层级由独立递归路径核对。原source的物理frontier不变，未创建child、未执行空间refine、未对H_eps使用同eps齐次。

|完整阶段|UH上端显示|I上端显示|合格g100 cost上端|32门|
|---|---:|---:|---:|---|
|088输入|6.678388883851e−5|7.030173690973e−5|52.56280830|未过|
|092首4096|5.868512887644e−5|6.211509377581e−5|46.44186489|未过|
|092全56139|5.669574350580e−5|6.010266982568e−5|44.93722704|未过|

E076完整范数外端固定，B=UH^(4/5)+E，I≤B_U^(5/4)，无064tube。精确32测试5^5 B_U^5<32^4g^5为false；并实际用整数幂核(B_U−E)^5≥UH^4、I_U^4≥B_U^5、cost_U^4≥(5B_U/g)^5。显示小数不承担证书。

剩余最大physical weighted gap是1141960779405953454237411570073918300617/45671926166590716193865151022383844364247891968，node6:62,61,61,61；top1024仅保存为remaining_priority供主选择下一有限空间批，没有自动派发。

## 4. 实际运行与完整独立读回

14个独立有限程序：前13个4096、最后2891；每个程序先冻结selection/parent state/receipt/code/input hash，12worker按节点原子checkpoint，close/join后退出。无wall alarm/timeout，无无限controller。所有成功节点只读已有结果不重跑；本轮真实fault数0、无续算，不声称测试过故障恢复。

首主PID79284（workers79298–79309）、末主PID80626（workers80650–80662，准确列表在receipt）。所有worker现已不存在。receipt wall合计41.996723秒，单程序2.075482–3.512291秒；该口径到pool join，最终cert/state落盘在其后，不与完整calendar wall混淆。

独立读回PID81026、workers81039–81050，17.812451秒且全部join。它对89个primary192根逐项验平方与组合幂向外，并用 **√(16+x²)、√(256+x²)的224-bit重标度根** 独立界定新R，不调用old ratio或H。全56139节点采用另一R外端、逐项288-bit有理向外lower累加（必要时exact Fraction fallback，本次0次），准确rational upper/mean、224-bit整数四次根凸弦，所有新lower用四次幂直接比较。全部4,098,147个band质量严格正，全部aggregate交集及self min/max正确。

最后独立递归重建14个全域prefix，校验完整覆盖、lower/upper单调和所有B/I/cost幂。可携只读复核另实测约3.27086秒通过全部56139凸弦、89根及14prefix，不读取外raw packet。

## 5. 源码、可携复核与最小归档

global registration4043674b8fb0987ad2d90a194035550fb740408d28b815a936b37991cc22156d；engine755bc2e77d435a87bdd00e9e48a500c72cd23921d56137a5a16a11bf385f051c；nodes9751ea75ef452b112a1420adfaf7efd6e5542f352088a905007194d6f20f46bc；catalogb33176d60b4c7942262169f14cf3bff7db40ad255527e5661831a7f9cc642e74。完整文件身份以manifest为准。

`python -B verify_portable.py`默认只读，使用本包及相邻冻结082/088包，不读取旧raw。作者--record结果在portable_verification.json。`python -B replay_one.py --rank R --output NEW_DIRECTORY`只在新目录纯重包络，需已冻结endpoint-only packet，不调用H；本轮未为了测试入口重复成功节点。

外部完整最终state：`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_range092/state.json.gz`，SHA7cf145d50a08b07a081e5a8faea6d9d827a5617757275aa478b3a7090d7292f6。节点表给出每个endpoint packet路径与SHA。原raw和stencil不复制进Git；每块量化证书、有限登记、receipt、R/root catalog、层级/样本/优先序、源码与报告/manifest足以复核或确定性回放。

旧准备plan/source_snapshot已移至外scratch archive_preparation；当前mechanism明确取代“可选per-band oldRU min”准备文字。所有旧079/080/082/085/088/091源最终核hash不变，无MATLAB、提交、暂存或推送。g105未用。下一阶段只准备095方案，须主在092/094门后确认才派新H。
