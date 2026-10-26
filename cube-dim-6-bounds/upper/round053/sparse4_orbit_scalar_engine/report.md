# C053 sparse4 orbit／scalar 预备交付

本阶段只完成新frame `body32+(1/80)sparse4` 的结构适配、轨道准确计数、16项形式标量预检和一次Arb512 scalar core编译。**完整frame kernel调用0，物理正T代入0，G/H/K数值0。** 未运行152方向的新真实S2资格检查，未继承body32旧资格。有限部分的物理解释仍待B052理论接受；本包没有代替F050完整fan独审。

当前交付是可审查的预备包，不是一个已经发布的full point evaluator。编译binary仅能运行ABI模式；full外层driver及其新注册须在后续单独release内完成。详尽适配、数学公式、成本与拒绝条件见 `adaptation_plan.md`。

## 输入与完整结构检查

预登记input SHA字段：`6e877703e420f9534987fb669b14ebf0b9267d38fa305f38a5fe212ebb3bd22d`。
只读输入050 manifest SHA：
`891160dc9e98571218bcb56717b374a6c4e81f05db58c5fe340c387d3544ae51`。
完整fan gzip SHA：`2405ccad9b2b0c344167539e769e3c64e3275e6e4ead208fbedc3167be47c948`；解压内容SHA：`522308c63bf203f4b5d6170ada14a81f3e824d0b014fe1d50aa78bb5d22795d6`。

实际读取全部478 primitive rays、13229 simplices。逐ray重建152真方向的加权绝对形式和，核对其真实q；逐simplex以准确整数Bareiss行列式核对原件。对于q=p/d，逐列同时改成`R=d r,Q=p`，每个node比值及每个完整det/product(q)权重准确一致。

|结构范围|准确最大值或计数|
|---|---:|
|完整ray／simplex记录|478／13229|
|逐列缩放d|40|
|整数R坐标绝对值|520|
|整数rate Q|20773|
|primitive determinant|32|
|scaled determinant|2048000000|
|拟议点的安全node numerator／denominator界|156000／1038650|

全部符合int64；六个Q的乘积仍须使用大整数或逐次Arb除法，不许做int64积。factor24与det也不能误用int32。完整转换数据及完整权重流的SHA分别为
`bbe9274e5b235d5c7f94b8c226ca124bd2c1f778ba012a124440cdf711cd47b3`、
`77f1584527cf64989a5cfdd62aeea08e51dc0cf7883462bda31cf789d44390f2`。
没有重跑050 active-system、fan覆盖或旧21cone baseline；这只是全部条目的适配检查。

## 轨道和形式核验

独立枚举46080个准确群作用，拟议点分子`(50,50,50,45,27,0)`、分母50，得到3840个像，每像恰12次；global-negative代表1920个，每代表权重24。组合计数及奇偶证明与实际Counter一致。因此未来实际point-orbit/cone数为25399680，逻辑项为609592320，不能按旧21cone模式运行，也不能在权24外再乘12或2。

H的DD5在全局负号下为偶；整个epsilon-Laurent DD5及每一系数同样为偶。零坐标的符号稳定子已在12中，负号配对无self-pair。本点不是零点，故该1920配对成立；不能将其盲用于其他点。

实际16项FORMAL checks为8组含零／重复／混号配置各两项：H在符号T下的导数与独立实resolvent乘积一致；G的Laurent向量（π/2常数系数及epsilon^-1,-3,-5）由三角DD与独立Hermite多项式权重准确一致，且全向量负号奇偶通过。全部残差准确0，未代入`T=32/825`或其他物理时间。全零配置虽有formal constant0，却保留非零负幂；本包没有把它称为有限G。

结构检查0.39823295900714584秒；16项形式批0.08207712500006892秒。两批各运行一次；没有追加点或参数。

## 编译与待放行边界

一次编译、0修复、退出0。新 `joint_scalar_core.cpp` 复用C010的H normalized-jet算法，并加入明确命名的FORMAL Laurent constant scalar routine。它不是旧fan driver。仅运行ABI模式，结果为FLINT3.6.0、ulong/fmpz8bytes、arb_struct48bytes、Arb512；没有调用其G/H scalar入口。形式预检验证数学算法，不能伪称已完成编译binary的新T数值交叉。

新source SHA：`3f2daa865ee5ccc218a9ee003f314eded13054e119b2f62014c20a05a7da0c48`。
新binary SHA：`2dd40f70c307f01b9dbe4a1092ae6defafaee165d3e0eff09e9c9b95512b1d24`。
compile、private rpath patch、匹配header/runtime哈希保存在 `build_record.json`。没有环境安装、旧源修改或系统路径变化；scratch binary及headers不入归档。

下一物理full必须全部满足：F050完整fan独审PASS；B052有限部分定理接受；针对152方向／11476pairs另行准确S2 preflight严格通过；scalar与外层driver集成验收；主代理单独release和有限注册。`s2_preflight_plan.py` 已提供准确Gram拒绝器设计，但CLI明确拒绝当前阶段执行。任一真S2 pair距离为0、Gram异常、遗漏pair、overflow、不可分辨准确nodes或hash不一致都须整体停止，不能跳项或回退body32资格。

拟议唯一T仅以精确式`32/825`保留，按marginal33归一，无优化。本批不自动接受并发出现的依赖文件，也不自动启动full。未来单点`K=G−H`仅是必要条件probe，绝不等于全域mT、g、I或Cn界。

## 只读封存

`verify_saved.py`检查保存哈希、完整结构收据、准确orbit计数、16份形式收据及Hermite权重、ABI记录，不重跑群枚举／13k determinants／形式H-G／binary或S2检查，不写文件：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round053/sparse4_orbit_scalar_engine/verify_saved.py
```

所有文件仅写本新目录。无kernel、旧baseline、扫描、grid、MATLAB、扩招、stage/commit/push。主保留后续证明、release和集成责任。
