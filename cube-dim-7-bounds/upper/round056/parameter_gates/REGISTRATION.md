# R056 单一proposed参数包注册

对象仅为一个由主协调者指定的新proposed refined inverse tuple，不是搜索：d=3/10，σ=7/100，χ=7/200，L0=2^20，ε0=2^-24；P0=65275/65536，εJ=εN=εRet=εS=εT=εD=1/1024，εA=1/64，M=21/8。target rlo=43/40+2^-20，rhi=49/40−2^-20；source q∈[31/40,37/40]。

固定运算：target log用alternating Taylor4 lower；pin log用Taylor3 upper；source exponential用Taylor8、τ exponential用Taylor16、scalar exponentials用Taylor8；baseline使用n≥L^8/8!，其他geometry可用n≥L²/2。只作一次候选完整Fraction/SymPy gate重建，不扫描相邻参数、不重跑旧grid。

停止/报告条件：任一domain、target概率、校准、moment、cap、τ或完整scalar gate失败，必须原样报告失败，不自动换tuple或把小数当准确证书；craw lower若不大于1/40，应报告真实lower并另选安全c。所有无限尾必须由解析单调性或全域Taylor上包闭合。仅constant/gate包通过不能替代尚未接受的新generic/operator证明。

原始输入来源是本轮协调消息，当前几何作者的新full theorem尚在构造。旧inverse-coordinate闭式只作公式来源；它不自动许可旧operator参数域跨到d=3/10。作者新源落盘后其实际hash应另行记录；不得伪造一个尚不存在的源hash。

所有源文件仅写本目录。无MATLAB、生产扫描、旧网格、Git或他人源修改。
