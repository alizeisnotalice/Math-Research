# 152：150 的唯一实际有限批独立审计

**73 个新 factor、64 个原根、4672 个缓存 band 质量乘法全部通过独立算术审计。** 不同路径精确匹配生产 150 的所有 factor、逐根 mean/gain、每个完整 UH prefix 和最终 UH。64 根的 mean 均严格改善，总 gain 严格正；本轮仅确认完整 H512 积分上包的改进，没有合成新的 G 范数、I 或 Cn。

唯一新审计 PID **73479**，exit **0**；checker 时间 **0.08927041600691155 秒**，完整 subprocess 时间 **0.12113116699038073 秒**，已退出。输出位于 `run001/`。外授权 `primary_authorization152.json` SHA 为 `a0840b85e061fc93e7029a6be86913c50f2b8b576ae5a4ac90b614900be55e0e`。准确外层执行计划 SHA 为 `d0d85e1898779740afaef188755c1ec7d32c47600dd2de69da07afcf74714657`，绑定 SHA 为 `7890a763c07a2ff86914add104c3e562c4194da556fd333912ff11f9469b4a3c`，原 base plan 为 `4db99fefe092b2cb838a0df5104749d899b8548127ff7c83fd66e4b72b0892e5`。

实际核得

`old_UH = 32424400796227698610199181381141800612723929 / 730750818665451459101842416358141509827966271488`，

`new_UH = 64847899372736011555499717104202734173681319 / 1461501637330902918203684832716283019655932542976`，

`gain = 902219719385664898645658080867051766539 / 1461501637330902918203684832716283019655932542976 > 0`。

算术核心未导入作者 consumer、Fraction、R 或 power root。时间权由 z³、z⁴ 的积分构造 barycentric 系数，与生产 theta 实现不同；64 个 mean 将原 H dyadic 分母统一成整数网格，再乘 factor 的 192 位整数并一次确定 160 位 ceil，与生产 Fraction 累加不同。原 112 的 224 位端点 power 未重算。用于安全减量的是 exact φ 弦斜率下界 1/16；真实旧 active arm、无 child root、orbit、体积和原 packet 身份已在准备来源门独立实核。完整理论及其适用条件保留在 `theory_report.md`，没有改写。

本次准确消费 150 的实际 PID 69275、exit0 输出，execution manifest SHA 为 `fc47797bb77735bbec86e30aff3e1efc7644b26438e55dba855db497a2ac58d4`，原 results SHA 为 `79428e616f3d96715a05413d2f7b13e023e99ddf532a235f837fb616966319d8`。101 个静态输入和 16 个真实输出/执行证据输入在运行前后全部 hash 相同，3 个代码锁不变；独审逐项比较所有新数据，没有以生产 passed 标签代替数值核验。

原准备 `manifest.json`（`f454bec0df48ba225169724ad193f740d90cfaba666554832cf107cb1eef6935`）和全部准备成员保留原字节；新实际资格入口是 **`execution_manifest.json`**。154、145、141 冻结件也未修改。完整有限结果、注册、退出、源锁终检及明确白名单见新 execution manifest。没有重启、追加 case、H/R/log/gradient 求值、旧 power 或旧基线重跑；无 Astra、MATLAB、下界、提交或推送。主代理保留最后接纳及任何后续常数合成责任。
