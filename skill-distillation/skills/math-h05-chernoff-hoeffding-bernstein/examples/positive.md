# 当前正例：独立有界和 Hoeffding

Claim: `H05-HOEFFDING-INDEPENDENT-SUM-01`（完整 MGF 证明见 evidence map 中 H05 的本地推导项；从安装后的 Skill 使用 `references/portable-audit-access.md` 解析该审计记录）。

取 `X_1,…,X_20` 独立且服从 Bernoulli(1/2)，令 `S=Σ_i(X_i−E X_i)`、`t=5`。每个区间长度为 `1`，故 `D²=Σ_i 1²=20`。精确尾概率是

`P(S≥5)=P(Bin(20,1/2)≥15)=21700/2^20≈0.0206947326660156`，

Hoeffding 上界为

`exp(−2·5²/20)=exp(−2.5)=0.0820849986238988`.

当前 case `H05-Hoeffding-independent-positive` 保留 runner 的 double-precision 数值近似 `0.02069473266601564`；精确有理值为 `21700/2^20`。有限计算检验该输入，不替代上面的通用 MGF 证明。
