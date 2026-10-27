# 第049轮：端点原算子门槛量化到log n≥2^52

已完整独审ACCEPT：对每个整数n满足自然L=logn≥2^52，Cn≤F(A_n(t),B_n(t)−2^-40n)<T−2^-41n，T=nL−nlogL，t=n/[(1−2^-43)L]。这里是全部机制条件闭合后的真实operator域，且适用于一般L1与全部正共同半径，不是R045比较的2^240标量门槛。

固定d2^-24、eps0 2^-20，沿用R045已证soft b=min(1,theta/G)与同history增量。n=e^L≥L8/8!给baseline log和的mean≤8![24d³/L6+4d/L7]、variance≤48·8!d³/L6。基点使mean≤1/32、variance≤9/8192，故log坏率≤1/8；Kw、targetcap坏率各≤1/8，保留5/8每面而得到A质量≥nd/16。pinned坐标只入均值，未当free或漏去。q0即使越source端点h也保完整正域，individual J、a.e/AC全q桥、inverse-row概率Markov均明确量化。

原sourcecap/k≤2^20Ld²e^-Ld/64，由Taylor4化为24·2^44/(L³d²)并全尾递减；基点<1/36。tau/k用Taylor64在4d²L0=64给2^74·64!/2^384<3/128，尾是L^-125/2。完整scalar的三个正项在基点由Taylor16给2^-103+16!/2^144+2^27·16!/2^128<2^-42，尾按L^-2,L^-16,L^-31/2递减；准确付delta2^-40及buffer2^-42，得到严格负2^-41n。没有把有限盒外推为全域。

beta/kappa、h/d、zeta t/h、source/target分离、参考R/g/扩展q0、n及nd、实际finite-a≥8n4的center guard与ca缩放全部量化。b≤1保真实同一target一次资格支付tau；原witness有限覆盖、紧水平集与内正则性给一般L1/全尺度，未启动其他history、support假设或选择器极限。

与旧完整cost比较：R042在自己L≥2^108域有Phi42≥T−n(4/L+2^-31L^-1/3)，该域全尾loss≤2^-106+2^-67<2^-41，所以本轮完整Phi严格低于Phi42，未只比展示负项。更旧小门槛R034等仍按各自有效域取min，本轮不宣称所有维数都排序。所有公开/内部比較边界保留，内部成果不称发表。

作者实际37项Fraction gates与独审15项不同Fraction gates均过，后者另逐项检查全解析尾及全算子桥；主也重建基点与所有guard，未重跑旧基线/网格/LP/下界/n7/MATLAB。准确源、最小代码及实际receipt在本轮source/audit目录；原R045未改。R048的30页核心稿在本轮前冻结，未纳入此显式域；最新合稿另轮推进。继续一般覆盖，旧暂停项目不重启。
