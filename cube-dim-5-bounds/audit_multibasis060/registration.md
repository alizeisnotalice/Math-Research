# audit_multibasis060：多符号基81点判据独立资格

预登记有限批`audit_multibasis060_integer_001`：只回读056选定的162个原split父盒，不重跑H、梯度、DD、log/exp/Gamma或旧043数值基线。不触052控制器、队列、STOP或缓存，不提交、不MATLAB、不改旧锁及其他代理文件。120秒软预算，先armed报告PID，按盒检查时间，每16盒外置checkpoint `/tmp/audit_multibasis060_checkpoint.json`；停止后不自动重启。

输入锁：batch0002/state.json `b79c134fdcec283c120d464cf78f131c669c905ac5cdacda309a34f07a9fb29b`、batch0002/registration.json `585a5de78b8f159ca03f138013894041dba350d144c3217b7bbda09fd5c7c9e9`、046/arb_cover.py `df3fb4165aea29e4aa3d1fb7ecacdb9cc5ad4a531bc03e3cc43a0820671fe028`、043/report.md `a298ca56ba980fce4fe2542929e03befb97b528afd3d390ed9d15651da68322a`、047/report.md `0b680c1bd7e4aea30c027fe32cc3e51436f7dda6dd8b74ddf521c2d74bb8ca16`。056最终manifest `0f3261b09dabc9557903547c0670125392e2bf106a8f4e7c3d952005b3978ba4`；report `7aee0b3a6f685f5e722f2abd95cbf303357ef7ef4d4dfe8659c0b34e487c5303`；results `564d454c0d0f2631843c1fcdff12f1b34596ba4697c177d1c6a2a6d69effa395`。数学回读前补记最终锁，默认入口仅核锁源，显式--run才执行。

理论审核：正Cauchy混合、可积支配两次微分、归一score方差给全空间Hessian下界；对每个合法diag(1,s)N分别积分，再逐点最大。共同物理gradient无需符号变换。tail曲率min_s Q_s=||v||²/4−||v||1²/36，每闭象限凹；全局凹性不是假设。anchor位于闭box内保证16个闭正交胞元的全部顶点包含于每坐标{lo−anchor,0,hi−anchor}的81点集合，处理所有零坐标和闭边界。

独立准确路径：仅标准库，无作者函数/SymPy/Arb。整数验证固定N与给定6N^-1乘积、16 signed basis及tail Gram。对每个盒将全部offset、Delta、gradient端点清共同分母到整数；由16闭象限胞元的角点并集生成候选。每点用16个signed inverse的五维平方和取最小曲率，再与l1闭式交叉，避免调用作者psi或沿用其Fraction循环。准确重建旧16点与新≤81点最小margin、worst offset及162回读结果；梯度ratio用正H端点交叉乘法验证。不重新求Delta的log，依赖冻结046/052球包围与Delta方向，报告共享算术边界。最多162*81*16个整数曲率评估；无新空间函数评估。

只得到这162个闭父盒的新资格；不改变覆盖树、不推出全域g或C5界。精简白名单：本登记、check.py、result.json、report.md（若需保存逐盒准确margin纳入result，不另复制旧原件）。
