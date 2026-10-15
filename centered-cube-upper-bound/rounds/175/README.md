# 第175轮：共同捕获储备与几何来源对包络

> 第180轮来源审计修订：ν≤β 的原始继承尚未定位，故 β−ν 为正残余及据此构成正 subprob 的部分暂属条件结论；同 α 的固定 ν 储备恒等式不依赖该排序。涉及 Wν→Wβ 的后续预算也需额外质量比较。详见 [第180轮](../180/README.md)，未修改本轮冻结归档。

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮独审通过：同alpha的bad来源捕获储备与正残余联合估计消去低alpha分母，按实际K0/wt比值切出的来源对获得根阶预算。精确标量包络扩大部分可付区域，未付低捕获补集完整保留。

双储备三点核保留真实outer/inner密度及j费用，但尚未证明全局列预算，不登记为新付款。

建设性数值：172个有界候选中，面上修正产生50个有限菜单下活跃R_eta，30个Gamma本身超过根阶阈值；n4/n16两例通过N128/256区间认证，并都在未付低捕获补集。actual Lplus定义已修正重跑。两脚本迁移复现，全部非runtime字段一致；14个有理包络极值×1001网格及7反解精确检查通过。原子/有限future不是一般L1、continuumFIRST或history认证。

[research.zip](research.zip)：87,624字节、17成员，SHA256 `1a9f32004a8b24de2706b565dab872cb65b8843b4ec5ddd694f8fa676699bf14`。含详细LaTeX、理论独审、数值及完整费用。未编译PDF；原工作已仓库外逐字节备份。仅新增本README和ZIP，并更新项目索引。

解压后 `python reproduce.py --verify-only` 核验字节；`python reproduce.py --output-dir <空目录>` 重放。依赖冻结本机核和vendor路径，非独立便携包。详见report.md、parent_verification.md及numerical_review.md。
