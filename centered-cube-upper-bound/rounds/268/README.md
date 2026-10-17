# 第268轮：真实硬来源比例带折扣与单预算饱和分配

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。

保原实际硬来源，在来源核心x∈Q_L(y)、R≥L的分支证明δ_j=min(1,(1+d/a)^nγ^(−2nj))折扣。原唯一j的完整后验响应费用成为Σδ_j^p=O(√n)。用正确未归一硬门tildeη而非错误加强旧条件η；保全部dx、原来源与后验。

这让同三腿来源容量从R267的λa^n/√J0提高到λa^n，仍获低soft、固定原合法参数下的O(√n log³(n+2))来源价。另证明可逐步增加来源、最终只使用一次接收预算的饱和算法，并给精确有限步rate。只登记旧已付域外的新来源增量。

剩余包括该核心分支的饱和邻域补源、核心外及R<L方向，均未付；不把容量缺额矩当交通，不改变原signed母账。

[研究与复现包](research.zip)：73,921字节，9成员。SHA256 `965f510a6e20a6ab434c7494a5f019186ced269cc3be260f9149c0178a555fd8`。

- `proof.tex.txt`、`report.md`、`audit.md`：详细LaTeX文本、成果和独立先审及终审；不编译PDF。
- `source_checks.py/json`、`ratio_checks.py/json`、`reproduce.py`、`manifest.json`：必要复现。

主稿SHA256 `973f9d4f05ceb5627daa91fba62dc2c814de97df2d92455317bfe7af5299f364` 已独审PASS。主代理完整重跑两套程序，JSON逐字节一致：11模型（含3连续L1来源例），1252精确窗口检查；24比例费用、14直接求和、44有理几何放置；另198点态rate和45缺额矩检查。50/90位最大差约2.502e-47。有限来源例不具FULL资格，高精度并非区间证书。

解压后装有mpmath的Python运行 `python3 reproduce.py`；`--check-only`只核哈希及rate。全部展开文件已核验备份于仓库外；本轮仅新增本README和ZIP并更新项目总览。
