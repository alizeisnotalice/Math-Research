# 第229轮：自适应共同窗口与局部门区间证书

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮新证AG5/AG11：原三测试腿在任意实际连续L(x)、R(x)、s(x)选择下，接收积分≤24e^4 N_n c_j^(2n)a^(−n)E_a(y,v)，N_n=O(√n log n)。点态向上节点覆盖后才应用固定QG；物理M^(−n)几何求和避免另付线性n尺度节点数。上端截断、两朝向、第三尺度与soft费用均保留，精细版保同一hard来源和捕获窗口的两乘积差。

AG13回代保β、ν四来源、q、j、a^(−n)、双后验包络及clock长度。粗价中soft与clock根阶仍相乘，ν总质量/(q a^n)及包络质量也未付。下一步是共同窗口与这些真实来源权的联合预算，不能把接收界当整体root闭合。来源交换辅助稿仅给有价身份，无新付款。

## 独立审查与实际计算

新桥、完整费用和局部平滑定理已独审通过。本机Fraction对n2/16/64/512、比例带/截断/两朝向/三腿执行61,066项断言。另有作者稿内878项精确检查。24组几何菜单按初探、压力与扩大范围进行40/70dps计算，只是有限探路，未拟合维数阶。

独立区间审查发现原n512候选遗漏J与K0容量且实际失败；保留失败区间。确定缩减mv至原值1/4、mz至精确1e−150后，37项明确列出的局部门在N256、256bits有向区间中全部严格通过。独立积分采用单调矩形与带显式二阶导误差的梯形包络交集，不依赖探索脚本的expint实现。原版本J/K0失败、第一修复K0失败、最终修复成功均在证书中。

主稿证明严格有限局部数据可平滑为正L1输入并保正Lebesgue邻域；这不构造continuum FIRST、原Z、forest/cut/rank或history。修复后H/q仅约1.79e−42，点态代理约1e−136，不能声称极值或全局付款。区间算术保留所定位mpmath实现的信任边界，非证明助手认证。

最终包外驱动完整运行三脚本，三个完整JSON均逐字节一致，15项manifest哈希前后通过。运行Python3.14、mpmath1.3.0；约1.46秒（grid）、0.38秒（有限探路）、0.24秒（区间证书）。详细回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round229_20261016_final_replay/receipt.json`。没有可归属本轮的实际扣额数据，未估算模型成本比例。

## 紧凑归档与复现

[research.zip](research.zip)：16成员，120,591字节，SHA256 `9de85ae157a3fe8dedf1a339a53cf6093b04bfff49615a95b314ef9ad5fca92f`。新增2个仓库文件，另更新项目索引；无PDF或临时字节码。

展开检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round229_20261016_final_checkpoint`；独立ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round229_20261016_research.zip`。逐成员核验work、检查点和ZIP一致后才移除本轮work。

解压后 `python3 reproduce.py --verify-only` 验哈希；完整运行需已记录的mpmath1.3.0依赖：

```sh
PYTHONPATH=/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor python3 reproduce.py --output-dir <新建包外目录>
```

驱动复制合同及脚本到包外再计算，不以读取旧JSON冒充复现。逐字节一致在上述运行时与库源码上核验，不承诺不同环境逐字节一致。详细LaTeX保留，不编译PDF。
