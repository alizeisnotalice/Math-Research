# A046：sparse4 首父 cone 的完整准确细分

**首父 cone 的局部 fan 通过；未生成整个框架 fan。** 对 beta=1/80 的真实152列框架，独立确认首 cone 仅有5个 crossing forms，索引2、10、26、50、82。完整32个闭 sign patterns 已处理：12个五维、20个低维边界，无空 pattern。结果为12条 primitive rays、14个全维子 cone，准确归一体积和1。

在父 cone 参数 u 中，五个 cuts 为

    u1−u3, u1−u3−u4, u1−u4, u1+u2−u4, u2−u4。

它们均不含 u5,u6，因此本父 cone 的非平凡细分是 first4 参数三维截面与公共边 conv(e5,e6) 的 join；这没有被推广到其他父 cone。

实际只做首 cone 的720个新 form/ray 整数 dot products，没有重跑旧2160基线。所有 sign patterns 共用相同边界矩阵，因此只需分类462个准确6×6系统，而非重复14784次：237个奇异，6个解在父 simplex 外，余解去重为12个顶点。全体 vertices 的各 pattern 可行性及 affine ranks 全部准确验证。20个低维 pattern 各给一个包含它的 full-cell witness，不删除零面。

pulling triangulation 直接使用11个半空间导出的准确 facets、统一全局顶点顺序和递归归纳；保存66个 face proof nodes。完整覆盖来自闭半空间 partition 和 pulling 归纳，并非只凭体积相等。另对12 full cells 的66个非空共有面逐一比较诱导 triangulations，全部一致。

14个子 cone 全部核对152个真实 forms 的弱符号、真实线性 symbol、primitive ray rate、比例尺度及物理行列式。没有使用旧 q。精确 rates 如下；完整 slice coordinates、active-set witnesses 和 determinants 在 local_fan.json。

|ID|primitive physical ray|q_true152|
|---:|---|---:|
|0|(1,1,1,1,1,1)|249/4|
|1|(1,1,1,1,1,0)|249/4|
|2|(1,1,1,1,0,0)|999/20|
|3|(2,2,2,1,1,0)|399/4|
|4|(4,4,3,3,2,0)|3989/20|
|5|(3,3,2,2,2,0)|2991/20|
|6|(8,7,6,5,4,0)|14959/40|
|7|(5,4,4,3,2,0)|4489/20|
|8|(6,5,4,4,3,0)|1097/4|
|9|(3,2,2,2,1,0)|4989/40|
|10|(4,3,3,2,2,0)|3491/20|
|11|(2,1,1,1,1,0)|1497/20|

每个子 simplex 的归一体积为其六个 u 顶点的 determinant 绝对值，14项精确和1。参数 u 截面 Euclidean5-volume 是 √6/120；物理 ξ 截面法向量 R^−T1=(0,−1,1,1,0,0)，故物理体积是 √3/120。写入后的一个体积标签混淆已作元数据修正，并由独立保存证书 checker 验证；没有改变任何 vertices、facets、simplices、rates 或构造算术。原始运行和修正后的 source/fan hashes 均在 results.json/run_receipt.json 区分保存。

登记后实际构造一次：PID48776、exit0、0.0262069702秒，Python3.15/python-flint0.9.0，所有线性系统和 determinants 用 FLINT 有理数。随后独立读取保存对象，使用 SymPy 准确 rank/determinant 重建全部66个 proof-node 关系、14个子 simplex、真实rates、物理截面归一化，PID48924、exit0、0.1154339314秒；没有重新枚举 active sets 或产生新 sign patterns。未用 scipy/floatHull。

限制严格守住：32 distinct patterns、462 active系统、14<100000输出 simplex；0 kernel/native/G/H/K/p、MATLAB、MC、优化、全 fan、扩招、commit/push。B044 report/manifest/frame 三个预期哈希吻合，所有旧输入末次 hash 不变。所有句柄已退出。

最低验收对象已经完整交付：一个父 cone 的合法新局部 fan 与准确体积。其余20个父 cone、跨父 cone 的全域整合、global multiplicities、zero/repeated Fourier nodes 和热核数值接口仍未证。不能将本局部结果交给旧21-cone引擎冒充全fan，也没有继承旧 tube/G 值或证明任何 Cn 改善。批次至此停止，主保留最终证明责任。

保存证书的无新枚举复核：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper/round046/sparse4_first_cone_refinement/verify_certificate.py
```
