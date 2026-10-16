# 第236轮：实际软参数几何与来源质量交换

**一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。** 本轮完成原far正余项的一条更紧全尺度估计，不把来源预算或其它分支视为已付。

保精确逐坐标后验KL余量，可用原接收点实际选择的软参数处理共同来源。原两槽FIRST与固定物理节点的hard核质量一列界给出真正的Lebesgue接收空间支付：

\[
\int S_{act}(y)\,d\beta(y)\le
\frac{\|\nu\|}{C_0}\left(\frac{\|\beta\|}{qa^n}\right)^\theta.
\]

随后使用同一完整原后验和增长重门，配合独立旧QG/Gini尾价，结算全部物理节点和比例带。令

\[
N=\lceil\pi\sqrt n\rceil,\quad
R_V=1+\log_+\frac{\|\beta_V\|}{qa^n}+\log_+K+\log(eN),
\]

得到

\[
\mathcal P_{far}(U)\le12\Lambda e^4(9e^{3/2}+2\sqrt e)
\mathbb E\sum_V\|\nu_V\|R_V^3.
\]

新界没有外乘N；菜单数只进入对数。原FIRST、来源、真实接收测度、时钟和全部维度/尺度费用均保留。关键是先固定hard物理节点再积分；不能把adaptive hard列免费归一。最终选择theta由outer来源一次确定，不能逐接收点优化。

**剩余缺口：** 尚未证明上述nu加权来源对数矩的根阶输入预算；Lambda、来源历史及其它分支还需完整回代。裸核两槽算子界也不能直接替代带逆权和后验的clock估计。以上不是一般弱型定理完成。

## 核验与复现

[research.zip](research.zip)保存详细LaTeX文本、两份辅助证明、独审、精简复现代码/结果、manifest；不编译PDF。本轮只新增本README与ZIP，更新项目索引。

本机实际运行：12个真实一维L1模型、12个一般维数后验恒等式、30组全尺度算术（498143个节点）；三阶段探路和压力/扩展复核。真实一维模型不是原large-n/FULL实例，普通求积加密差不是区间证书，不拟合n阶。另有可复现的120项两槽标量/有限算子诊断，独立计数。独审检查量词、端点、常数和同母量回代，独立重算两完整JSON及120项脚本。

解压后运行 `python3 reproduce.py --verify-only`；有NumPy2.3.5和mpmath1.3.0的Python环境下，运行 `python3 reproduce.py --output-dir <新建包外目录>`。驱动将代码与合同复制到包外，重算两完整JSON并逐字节比较，执行内嵌120项且核对输出，manifest运行前后均检查。具体本机Python3.12.14命令见包内numeric_notes.md。

使用F05核对两侧真实底测度/核界，使用E03核对历史质量流条件；没有据技能名称假借未验证的Schur/Bellman定理。没有可归属本轮的实际额度扣减数据。

## 冻结回执

包外完整重放通过：两份JSON逐字节一致（24配置约1.57秒、30配置约0.16秒），内嵌120项输出一致，12项manifest前后通过。ZIP共13成员、46,673字节，SHA256 `567c799363620655911e46bec2f447d1b3b9b4fa37b1a71e71b0046f2e78dd57`。

检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round236_20261017_final_checkpoint`；ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round236_20261017_research.zip`；回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round236_20261017_final_replay/receipt.json`。已逐成员比对工作文件、检查点、ZIP和备份后移除本轮临时work；无PDF。Git pack约366.94 MiB，本轮新增包远低于单文件限制。
