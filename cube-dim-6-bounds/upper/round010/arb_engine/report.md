# C010：FLINT/Arb512 编译核引擎与单一旧点交叉验证

结果：引擎可用，唯一既有完整点的512-bit区间严格位于旧Arb384证书内。最终同一binary暖启动外墙1.26252秒，对旧A007记录6.87921秒提速5.4488倍；内核1.25710秒，提速5.4723倍。首次冷启动外墙仍有明显代价，本批只认证这一点的吞吐，没有声称所有节点类型都有统一5倍速度。

这是新实现验证，不是重跑旧研究基线计数学进展；没有新kernel数学界或C6界。6个不同标量DD锚全部通过，完整群共跑3次，始终同一旧空间点、同一T=1/250；0个新点，已按上限停止。没有MC、MATLAB、下界工作、扩招、提交或推送，未改其他clone/旧目录。

## 1. 本机依赖与可复现编译

实际版本：python-flint0.9.0，wheel带FLINT3.6.0、GMP6.3.0、MPFR4.2.2；Apple clang21.0.0，arm64。FLINT动态库来自现有venv的 `site-packages/flint/.dylibs/libflint.24.0.dylib`，未安装或替换任何环境。

本机无已安装FLINT头文件。仅在本轮私有scratch取得 [官方FLINT v3.6.0源](https://github.com/flintlib/flint/tree/v3.6.0) 的匹配头，源archive SHA256为4307a504622702bf0be6d8969791f7d7ff378645cf2ae3bb5a7a2b56653d97f1；GMP使用本机匹配6.3.0头。曾取官方MPFR4.2.2头核版本，但最小编译依赖不需要它；没有将无用头加入交付archive。

flint.h与flint-config.h从官方模板在私有目录生成：64bit、小端、long limbs，关闭asm/gc/blas等可选内联优化，启用线程/TLS声明。它不是伪称wheel完整build-config已获取；只使用匹配版本的generic C API。C++编译期检查sizeof(ulong)=8、sizeof(fmpz)=8、sizeof(arb_struct)=48；启动还检查flint_version严格为3.6.0。6锚及整点跨实现验证检查了实际ABI/算术链，没有用binary64或平台long double替代Arb。

`headers_min.tar.gz`只有实际compiler dependency列出的13个头，约44KB；每个header与三个runtime dylib的hash在dependencies.json。LGPL/GPL许可文本保存在licenses目录。完整176头暂存不在归档白名单。

复现编译：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/upper/round010/arb_engine/build_engine.py
```

脚本解包最小头，核header及本机dylib的SHA，再执行/usr/bin/clang++ -O3 -std=c++20，直接链接已有FLINT。现有wheel的LC_ID含构建机/DLC绝对路径，脚本仅用install_name_tool修补**本轮binary**到@rpath，不修改dylib、系统PATH或配置。完整命令在build_record.json。最小archive独立重编译已成功，binary SHA与最后实际核验binary完全一致；重编译没有再调用kernel。

## 2. 精确输入、群项与归一化

输入fan副本 `fan_integer_input.txt` 与主clone原件SHA相同：7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812。保持全部21 cells、720坐标permutations、64个符号，逻辑总数967680。

每cone的权为det/∏qj，节点为 b_j=(Σyi rji)/(den qj)，最后除(2π)^6，与A007 run_moments.py一致。节点分子/分母先在整数域约分，排序与相等判定均使用精确有理数及__int128交叉乘积。参与模型值的运算全为Arb512；double只用于chrono的时间输出。

实现只接受节点整数可放入带符号64bit、交叉比较可放入128bit的输入；超界立即报错。当前fan/点严格在范围内。CLI的T_NUM、T_DEN、X_DEN前提均为正整数，坐标分子为相应范围内整数；新T=0.048可直接用6/125，但本批未运行该T。

流式版本使用A007已有的精确全局正负成对：对f(x)=atan(x/T)，五阶DD满足f[−b0,…,−b5]=f[b0,…,b5]（f奇、阶数5）。所以每对(mask,mask xor63)计算一次、权重严格乘2，实际DD483840次，保留967680个逻辑群项。它不是删掉负号项。所有permutations、21cells、confluent/zero节点均保留。

## 3. Hermite DD及重复/零节点处理

新引擎不使用浮点阈值跳过近节点。节点准确相等时，以真实Taylor jet替代除差；否则执行Arb除差，若精度不足以把两个准确不同节点分开（差球含0）就停止，不给伪精确结果。

标量jet用独立于旧Acb幂公式的实递推。记D=T²+x²，g(h)=T/[D+2xh+h²]=Σgm hm，则g0=T/D，g1=−2xg0/D，gm=−(2xg_(m−1)+g_(m−2))/D，atan((x+h)/T)的k阶Taylor系数为g_(k−1)/k。此式也包含x=0与全部节点重合。

6个登记锚包括：普通不同节点含0、一对重复0、6个全0、6个全1/3、三对重复含0、两组不同重数有理节点。C++512结果的完整dyadic端点全部包含独立旧A007/Acb640参考的完整向外dyadic端点。相同6案例的参考端点已补全序列化，没有添加新锚。精确参考见anchor_validation.json；anchor_native.json保留C++端点。第二版只改变群driver；canonical/Node/DD函数文本hash与首版一致，因此未重新扩充锚批。

存档复核曾尝试拿640参考的100位展示串重造区间作包含比较，因展示舍入使其区间变宽而失败；已改为保存同6个参考的完整向外dyadic端点，原生结果未放宽。这次格式错误没有作为任何数值证据接受。

缓存键是约分后的精确节点值，每个点/T内只缓存12324个不同节点的atan/必要jets。512-bit算术缓存不混用其它T或空间点。早期整张DD键表版本需要约124MB；最终流式版约11MB。

## 4. 可靠外包端点

所有新数值输出调用 [官方arb_get_interval_fmpz_2exp](https://github.com/flintlib/flint/blob/v3.6.0/doc/source/arb.rst)：它把Arb所代表的整个区间准确写成整数端点[a,b]乘2^e，不从短radius展示反造区间。输出三个字符串 lower_mantissa、upper_mantissa、common_binary_exponent。

最终H_(1/250)完整端点在full_native_v3.json，公共指数−523。用于阅读的值约为2.24792268996421778043582023348519848438597×10^−5；该展示不作为端点。verify_saved.py用完整dyadic端点复核6锚和三次整点记录，不调用任何kernel。

旧A007现有证书按其保存的Arb384字符串原样读取作交叉参考，没有把它改称新引擎内部区间。所有512-bit新整点区间都**严格包含在该旧证书区间内**。流式版两次运行的完整dyadic端点逐位一致。

## 5. 唯一旧完整点与有限性能结果

点分母1000000007，分子

`(1000000007,998836847,892738545,807464540,467291939,96968636)`。

T=1/250，旧A007记录完整H证书耗时6.87921229202766秒。每次都严格核逻辑群数967680；该既有generic点本身zero/confluent群项均0，特殊节点由前述6锚验证，没有假称整点测试覆盖了full-group退化点。

|新实现运行|外墙秒|内核秒|旧/新外墙|RSS|
|---|---:|---:|---:|---:|
|v1全部群键聚合|1.54207|1.50618|4.46103|约124MB|
|v2流式新binary冷启动|1.67599|1.30479|4.10456|约11MB|
|同v2 binary暖启动|1.26252|1.25710|5.44880|约11MB|

首版未达5x，已保留失败benchmark与压缩源文件；流式版内核本身超过5x，冷启动进程开销使其首个外墙仍未达。最后暖启动无重编译，达到登记的5x吞吐门槛。只有本机同一旧点的三次有界实现验证，不能推断其它点、全confluent点或并行压力下的统一时间保证。由于3次full-group上限已到，未再测试、未增加新点。

最后可直接使用的CLI（工作目录为本audit clone）：

```sh
cube-dim-6-bounds/upper/round010/arb_engine/scratch/arb_kernel full \
  cube-dim-6-bounds/upper/round010/arb_engine/fan_integer_input.txt \
  1 250 1000000007 \
  1000000007 998836847 892738545 807464540 467291939 96968636
```

这是复现命令，不表示本报告又执行一次。父任务可复用引擎做其另行登记的T/点批。本轮未提供梯度、带限差Hε−HT或连续区间x-box evaluator，只交付任意合法正有理T、准确有理点的H_T；需要band时可对另行授权的两个H证书作区间减法。

## 6. 验收、归档与停止

编译和静态加载exit0；6锚PASS；3次整点exit0且严格区间吻合。最后运行工具chunk490a69，外墙1.262519500014605秒；没有活跃进程。build_record.json对应从最小13头archive重建的同一已验证binary。verify_saved.py最终实际PASS，0次kernel调用。

最小必要交付为C++源码、build/reference/只读复核脚本、匹配头archive与许可证、fan输入、依赖版本/hash、登记、6锚完整端点、3次benchmark完整端点和验证记录，以及本报告/manifest。私有scratch中的binary可直接使用，旧全量头和编译临时材料不列入归档白名单。未改系统path、未重装库、未写其他clone；C009及更早目录保持冻结。

本轮源代码、二进制、fan、旧证书、官方archive及各header/runtime dylib的SHA均在dependencies.json、validation_results_v3.json和manifest.json。没有新的数学上/下界，未登记为C6突破。
