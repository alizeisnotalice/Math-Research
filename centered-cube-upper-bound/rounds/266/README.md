# 第266轮：固定来源几何与原后验联合估计

一般 n^(1/2+o(1)) 弱型上界尚未闭合。

本轮证明两移动核的固定 Δ_b 缺口；与原 good−bad 正交通及原后验直接配合，完整辅助后验/clock 行价为显式 O(log n)，去掉旧初始 heavy 比值，不新增删除费用。J 腿另获统一固定窗口半密度损失及完整辅助槽 L^(3/2) 估计。真实 dx、维度、soft网格、比例带、β、叶/种子费用保留。

仍缺外层来源平方/分数次方预算及真实接收重数结算，不能把本轮行估计称为一般闭合。详见包内 G13–15。

[研究与复现包](research.zip)：74,117 字节，8个成员。ZIP SHA256 `818c1c51df8ea7080be134c226a674c42e09550155d533ebc7d14e00ba7caff4`。

- `proof.tex.txt`：完整详细 LaTeX 文本，含原公式槽位、解析证明、常数与剩余费用；不编译 PDF。
- `report.md`、`audit.md`、`source_reciprocity.md`：成果、独审与原来源关系核查。
- `numeric_checks.py/json`、`reproduce.py`、`manifest.json`：必要数值及完整复现。

独审 PASS 绑定 proof SHA256 `e427ea2f3a1071706ccd81b21c7010c0f286aa47a00dabd13eb64867d8bfde97`。主代理实际重跑12组真实核与9组有限后验模型，JSON逐字节一致；359字段40/70位最大差约3.44e−41。另有精确有限恒等式及J标量检查，明确不作为FULL或一般阶数证明。

解压后用装有mpmath的Python运行 `python3 reproduce.py`；`--check-only`仅检查保存产物。详细环境、范围、失败修正和未完成费用见报告。本轮只增加本索引和必要压缩包，并更新项目总览。
