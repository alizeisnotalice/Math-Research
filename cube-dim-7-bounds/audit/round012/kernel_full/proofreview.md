# R012 full-generic独立证据/接口审计

2026-10-27。独立clone起点 `75eedb49b0e95dff748081d658d5a9f53a593974`，仅写 `audit/round012/kernel_full/`。已完整读新wrapper generator、checker、registration、result、receipt及报告；R011原数学/实现已另审，今核其冻结字节与完整执行接口。本审没有新DD、phase或whole-kernel计算。

**结论：接受固定exact-dyadic index2710、epsilon1/1024的完整真实H_epsilon点证书。** 独立Fraction求和5040个saved signed chunk intervals，再以768bit转换端点，证明

`0 < L_ind ≤ H_epsilon(x) ≤ U_ind`，

reported full-H interval包含[L_ind,U_ind]，且其宽及独立宽均严格小于1e−12。完整独立有理端点见results.json；值约3.77829151969660618696332123e−6仅作阅读展示。没有global g、新I或Cn界。

## 1. 原数学与冻结kernel

R011已审DD6正号、exact-repeat Hermite、real/complex等价、B^6及alpha-q尺度、完整R010135cells/212simplices和B7覆盖。R012 `kernel.py`逐字等于原R011 `run.py`，manifest及原文件SHA均相同；copied cells、rays、64prefix inputs、testcases和R010manifest也与原pilot的锁定项相同。没有改变q_i、alpha、common denominator、ray scale或模型。

wrapper只import kernel；其 `if __name__=='__main__'` 未触发，因而不会偷偷重跑R011 scalar/symmetry/pilot。每个new task只调用 `kernel.generic_chunk(j,256)`；该函数保留212signed simplex请求、exact integer knots、原phase-local caches和256bit Ball算术。无promotion、fallback、phase clipping或kernel重放路径。

每个generic chunk在原函数内已乘 `2/(2π)^7`，opposite half-group fibre为1。wrapper最后只加这些signed已归一intervals，没有再乘2、group order、pi或phase count。DD时间积分已准确支付epsilon至infinity，没有另一个数值frequency/time tail待补。

## 2. 固定physical point与lex全覆盖

完整point_floathex与原R011 testcases、R003 `truncated_numeric/results.json:minima[0]` 逐项相同，lowest_sobol_index准确2710，epsilon准确1/1024。minchecker将这些binary64 hex恢复为exact dyadic Fractions并保存在results，不用decimal输入替换。七绝对值不同且非零、各coordinate正，因此fibre1和first-sign+ half-representative规则适用。

lex permutation_rank0..5039，每个64 lex sign patterns，phase_id64j+s，构成322560half phases。wrapper condense在每任务返回时检查64个phase IDs准确连续；保存的terminal ranges和IDs覆盖所有5040chunks、无重号缺号。既有R011 header的phase记录保留在copied旧results中，新5036只保留condensed chunk certificate；新chunk内部phase执行依赖冻结generic_chunk及condense断言，本审不声称从未保存的phase balls独立重建它们。

R011原results文件逐字复制；前4 normalized interval payload（含准确endpoints/display）与其原chunk0..3完全相同，counts由原64phase records准确加和匹配。R012终态row因condense去掉phase_records、CPU改null并加source标签，不是整个旧record逐字复制；科学interval本身准确复用，未重算256旧phases。剩余5036个source=new，每个一次，没有隐形旧pilot第二次加入。

## 3. 准确signed求和及768bit端点

本审逐row将dyadic (mantissa,exponent)转Fraction，准确求

`L_ind=Σ5040 chunk_lower`，`U_ind=Σ5040 chunk_upper`。

其中2526个chunk intervals严格负，2514个严格正，没有跨0者。全部signed terms保留，未abs或clip。前256phase partial比full result大并不矛盾，缺失signed chunks可抵消；本轮如今已完整相加，才是H。

保存的256bit reported interval满足准确关系 `0<L_report≤L_ind≤U_ind≤U_report`。原reported宽≤1e−12，独立准确宽是

`1003292887117513 / 30354201441027016733116592294117482916287606860189680019559568902170379456331382784`。

768bit Arb只用于将两个exact signed endpoint sums向外转换；它们本身是可准确表示的dyadic，转换无损，768端点恰等于原Fraction sums。这消除final sum rounding，未将每个chunk的256bit arithmetic升级成768bit，亦未重算π、DD或任何phase。

本次强检查是准确包含而非overlap，因此直接验证保存reported full-H包住全部chunk certificate之和。真实H属于这个和仍依赖被冻结generator的每个chunk外包、完整fan/phase数学；本审没有把最终端点相加当作第二套whole-kernel实现。

## 4. 计数实测/重建分类

counts来自原kernel在每phase实际记录的cache逻辑分支；wrapper将它们condense成chunk总数，本审再准确加和。它们表示各phase cache context内的操作/unique计数之和，不是跨所有phases的global distinct keys。

|counter|全值证书，含复用|R011已执行复用|R012新执行|
|---|---:|---:|---:|
|top DD requests|68,382,720|54,272|68,328,448|
|unique sub-DD cache misses|659,175,346|517,440|658,657,906|
|log evaluations|14,837,760|11,776|14,825,984|
|cache hits|697,882,546|548,160|697,334,386|

212×64×5040也是top requests的理论identity check；其它总数由saved实际counters重建。本审没执行这些操作。Hermite、zero-ray-knot和repeated-top-node counters本point全0；这些来自被执行的完整generic编译/phase逻辑，不是凭“generic”名称推定。无epsilon0 G运行或额外physical point。

## 5. CPU、并发和内存的实际范围

receipt记录5036new chunk task CPU约898.730397秒、parent CPU3.332712秒、body wall133.678526秒。原4pilot CPU未再次计入R012新成本。按saved task记录重建并结合parent，平均约6.748 core-equivalents；这是定时器范围内的工作记录，不是全机器CPU积分，child import/setup等范围也不全由chunk timer覆盖。

本审只解析保存的resources日志，核26个worker阶段sample，各sample可见15个spawn workers，sum ps CPU从549.1%到713.7%，sum worker RSS354992到491968KiB。看到15进程不等于15核满载；报告没有理由声称full saturation。RSS可能重复shared pages，单worker高水位40,550,400bytes不是系统或全部workers的总峰值。机器共享负载/采样限制应继续保留。

## 6. 可携带checker、实际本审运行和边界

本审minchecker.py为CLI stdout-only，不import kernel，不执行author checker，支持--numeric-dir/--pilot-dir/--theory-dir/--point-source；默认合并后的项目相对路径。prehash的已知绝对theory/pilot路径通过这些目录参数重定位，不要求原电脑路径存在。

一次运行exit0，0failure、0repair，审计脚本0.13776762498309836秒；Python3.15.0、python-flint0.9.0，唯一Arb工作是上述768bit端点转换。所有source/pre-run/manifest和prior-pilot映射通过，前后hash不变。results保存完整independent endpoints、coverage/reuse、counts和resource bookkeeping。

这验证了真实截断kernel的一个固定point interval。H_epsilon≤G，低H不能给G上端或排除原模型；一个H point lower也不是整个face essential g。没有新I、g或Cn结论，没有MC、MATLAB、commit/push、新phase或whole-kernel replay。任务到此停止，主集成人保留最终证明与采用责任。
