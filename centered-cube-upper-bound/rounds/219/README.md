# 第219轮：三来源共同重叠与全位移几何衰减

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮建立同尺度、任意独立双软度、任意来源平移的一维共同重叠定理。记a=(1−s)(1−σ)、A(d)=(1−|d|)+、V=h*φ、C=h*max(h,φ)，则任意长度1窗口上的半密度积分满足 N²≤aA+(1−a)V≤C。证明使用原指数混合的共同参数单调正相关与三角不等式，不假设极值窗口居中。

进一步，c=1−φ(0)>0给全位移界 N²/C≤1−c min(a,1−a)，固定同尺度张量得到exp[−nc min(a,1−a)/2]衰减。远位移|d|≥1时N²/C≤1−a。a=0或1的退化完整保留：同尺度J腿虽有常数一界，却不能据此取得正指数衰减。

一般尺度部分精确保留原硬捕获门1_{R≤L}1_{Q_R(x)}(v)，接收几何是两个坐标乘积之差；没有提前按坐标求和。DG2–8给固定物理参数重叠及原good门的合法点态回代。异尺度、自适应soft/scale、全硬来源质量与两份后验分母、bad端联合score预算仍未支付；本轮没有登记新的原正交通付款。

## 验证与复现

- 完整独审通过，726行LaTeX文本保留，不编译PDF。主代理核对43个标签及35处引用、环境配对和指示函数记号。
- 本机Python/mpmath：36个同尺度配置、72个扩尺度三目标腿配置，各两精度/求积阶数；未拟合一般阶数。强式检查另复用36个积分值、10个解析端点及1260个Fraction抽象代数网格。
- 主代理从最终checkpoint执行driver，再迁移脚本与绑定输入到新临时目录。15项manifest哈希运行前后通过，三份完整JSON逐字节一致；运行时间分别17.833、3.138、0.060秒。有限求积不是有向区间或一般L1证书。
- 最终回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round219_20261016_driver_replay/receipt.json`。

## 精简归档

[research.zip](research.zip)：16成员，104325字节；SHA256 `79391c0ed887267c6b5c430f1d5692de8ee67135a9b9f2b8a73dae1b61df5543`。包含详细LaTeX、原几何推导、解析新引理、独审、三套脚本/JSON、数值笔记、任务合同、报告和manifest。

展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round219_20261016_final_checkpoint`；独立ZIP：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round219_20261016_research.zip`。逐成员与工作文件字节一致核验后移除本轮工作目录。仓库只新增本README和ZIP，另更新项目README；未改历史轮次或其他项目。

解压后运行 `python reproduce.py --verify-only`；完整复现用 `python reproduce.py --output-dir <新的包外目录>`，需mpmath1.3.0。scalar脚本需要相邻triple_probe.json，driver已处理迁移。提交与推送是否成功以Git及远端确认记录为准，不由归档存在推断。
