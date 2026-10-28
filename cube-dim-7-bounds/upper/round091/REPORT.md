# R091：一般配比的显式 1.1 覆盖上界

FULL ACCEPT（内部严格结果）。所有整数维数n满足自然log n≥2^33时，

    Cn≤H91=(1+A_n(t)/tau)(tau+B_n(t)−gamma n)<T−11n/10，
    T=nlogn−nloglogn，t=n/[(1−sigma)logn]，tau=n exp(−dlogn)。

固定准确参数包括s=4257987173/10^8，sigma=92808239/10^8，
d=46403803/50000000，gamma=373319131/10^8；其余完整分数见exact_registration。
一般L1、Lebesgue、中心轴平行共同半径及全部正尺度保持。

## 构造与全尾认证

使用已接受R090的任意固定s局部几何、五事件、weightedcap和自由tau完整H成本。
新阈值逐项支付D_s=4s²+6s+1、J1'≤5s/rho、目标s^4矩、source真variance2s²x
及pin1、e依s的计数、source>4(s+1)和同history位移。实际d_mu由log定义，
准确区间仅用作预算。旧s≤5常数、.95 cap底数、sqrt tau都未移用。
独审用N12/T8与range-reduced log20给出基点和每一全尾负幂/单调论证，
包括新tau的L²和完整H交叉项的L因子，不以单点PASS替代全域结论。

## 实际运行与准确来源

一次注册本机15worker有限探路，22,554评价；13区间有候选、2区间本批未找到。
wall约.764秒，5次采样最高合计CPU986.4%、RSS1106208KiB；采样非完整峰。
完整保存字段按预登记8位十进制ROUND_HALF_EVEN定义一个新有理tuple，
不是从展示小数恢复旧证书。原summary和取整规则保留。

首个预登记门槛2^33由本机PID33459实际58项Fraction PASS，因此未运行备用2^34。
独审另注册N12/T8/log20，PID34233，7组PASS，未复跑作者58或15worker。
主读取已归档完整成本区间，注册一次2项新的准确比较PASS。所有原始身份、
来源hash、失败分类和运行范围分别保存。有限浮点搜索不是最优性证明。

## 完整比较

旧H89≥B89−gamma89 n>T−[c89+20000/(4987L)]n。
准确原件给c89<.901，本轮L≥2^33上H89−H91>.1989n。
旧固定H88类似给T−[c88+40/(11L)]n，故真实共同域有H88−H91>.1969n；
该比较须满足H88未知eventual门槛。旧成本下包不是Cn下界。
本轮1.1的门槛2^33与R089的.9/2^27按域取小；R090的更强1.181387…
limsup仍有效、尚无显式门槛。未宣称新固定tuple支配整个R090参数族。

公开比较沿用已核文献口径，前导nlogn系数不变。本轮是内部显式认证，
未新检索文献、未称公开发表或全球最优。

## 交付与后续

全尾独审在audit/round091/threshold/review.md。必要原始注册、receipt、summary、
最小复核代码与来源hash归档；可再生15份worker细目留原独立clone，不重复上传。
原manifest保留完整来源索引。程序保留原路径身份：复核请在副本运行，必要时
将配置中的原proof_path、SOURCE、BASE指向同hash归档文件；不要改历史receipt。
exact_verify.py的单点参数为8589934592，独审与主比较check.py各有独立范围。

R092局部单侧引理已主审，条件缩放模型没有新Cn界。R093 Bernstein及新宽度
构造正在推导，R094/R095较低门槛均仍为浮点候选。论文不改，下界暂停，七维专项停止。
