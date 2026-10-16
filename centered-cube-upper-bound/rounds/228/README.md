# 第228轮：真实捕获几何与共同硬窗口

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮证明全参数点态几何二分：H0≥a时强捕获余支排除闭Hhat立方体，因而排除统一a核；H0<a时远时钟长度受cφs−(u+0.5logell)/n控制，其中cφ<2，明确低软度范围严格为空。所有截断状态与等号保留，高软度剩余支未删除。

保留原hard来源时得到共同窗口判零及固定参数积分的两乘积差公式；它不能无代价替代自适应参数积分。真实来源包络及接收覆盖的联合根阶预算仍未证。阶段分层核查没有新增付款，完整费用与缺口保留。

## 实际核验

独立审查先审命题，再读证明，核实量词、截断、边界与费用。几何脚本4组共同窗口、3组边界例，共34项断言通过。空间探路使用真实核、独立pilot冻结q与新seed测试；n=8/16/32，每模型2048个直接样本及2048个移位样本。补入原模型hard因子后联合正hard事件均零命中；仅有限采样证据，不证明事件不存在。移位采样低ESS作为失败诊断，不用于更强结论。无维数阶拟合、无FULL认证、无PDF编译。

最终包外完整复跑：两个完整JSON逐字节一致，12项manifest哈希运行前后通过。Python 3.14，mpmath 1.3.0；几何约0.97秒，空间约30.31秒。回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round228_20261016_final_replay/receipt.json`。

## 精简归档与复现

[research.zip](research.zip)：13个成员，64,829字节，SHA256 `6e86d15cc6a304e2f24c9ffe4678da811e11eb4a9de2df80549c5610b9e77743`。新增2个仓库文件，另更新项目索引。

展开检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round228_20261016_final_checkpoint`。独立ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round228_20261016_research.zip`。work/检查点/ZIP逐成员一致后才移除本轮work；临时字节码不入包。

解压后可用 `python3 reproduce.py --verify-only` 验哈希。完整复现需要mpmath 1.3.0，在本机实际使用：

```sh
PYTHONPATH=/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor python3 reproduce.py --output-dir <新建包外目录>
```

普通未配置依赖的Python不能直接运行数值脚本。浮点JSON逐字节一致是在上述同运行时核验，不宣称跨平台逐字节一致。数值证据与一般解析证明分别标注。
