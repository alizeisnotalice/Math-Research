# 51个声明节点的语义依赖索引

绑定初始数学快照205596954c5a3265532b8dfbda0354e20561dd68；四份数学源和PDF在本轮保持原字节。索引来自三份实际重建，root补齐8个原索引未列的incoming边，并纠正一般上界支撑数不受下界晶格“偶数平衡”限制的转写。它不是自动生成的数学证明或外部专家裁决。

精确假设、量词、误差、实际对象、证明位置和逻辑依赖逐节点保存于[JSON](063-dependencies.json)。全部51项已有映射；TeX语法引用单独保存，不作为逻辑证明边。

| 声明标签 | 源位置 | 重建节点 | 解析记录位置 |
|---|---|---|---|
| thm:endpoint | paper/main.tex:127 | M4, M6, C1, C2, C3, S9, S10, S11, E0, E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13, L0, L1, L2, L3, L4 | reviews/063-foundation-reconstruction.md §1 M4 and §5 C1; source paper/main.tex:127-136; reviews/063-eighth-reconstruction.md §§4-9 and node table §4; source paper/main.tex:137-185; linked theorem paper/sections/eighth.tex:9-25; foundation proof_record.md §5 C2-C3 + eighth reconstruction §9; source paper/main.tex:151-164 |
| thm:reduction | paper/main.tex:198 | M5 | foundation §1 M5; source paper/main.tex:198-205 |
| lem:defect | paper/main.tex:228 | M2 | foundation proof_record §1 M2; source paper/main.tex:228-234 |
| lem:partition | paper/main.tex:241 | M2 | foundation proof_record §1 M2; source paper/main.tex:241-248 |
| lem:bip | paper/main.tex:263 | M4, B0 | foundation §1 M4; sixth rebuild §2 B0; source paper/main.tex:263-265 |
| lem:finitegreedy | paper/main.tex:273 | M3, B0 | foundation §1 M3; sixth rebuild §2 B0; source paper/main.tex:273-283 |
| prop:upper | paper/main.tex:295 | M4, B0 | foundation §1 M4; sixth rebuild §2 B0; source paper/main.tex:295-301 |
| lem:parameter | paper/main.tex:314 | C3 | foundation §5 C3; source paper/main.tex:314-329 |
| prop:lower | paper/main.tex:415 | M6, G0 | foundation §1 M6; sixth rebuild §13 G0; source paper/main.tex:415-421 |
| cor:uniform | paper/main.tex:480 | C3 | foundation §5 C3; source paper/main.tex:480-495 |
| cor:eighth-multi | paper/main.tex:509 | C3 | foundation §5 C3; source paper/main.tex:509-516 |
| lem:global-small-gap | paper/sections/global.tex:52 | G1 | foundation §2 G1; source paper/sections/global.tex:52-65 |
| lem:global-linear-radius | paper/sections/global.tex:123 | G2 | foundation §2 G2; source paper/sections/global.tex:123-130 |
| lem:global-face-capacity | paper/sections/global.tex:169 | G3, B1 | foundation §2 G3; sixth rebuild §3 B1; source paper/sections/global.tex:169-203 |
| lem:global-face-neighbor | paper/sections/global.tex:431 | G4, B1 | foundation §2 G4; sixth rebuild §3 B1; source paper/sections/global.tex:431-441 |
| lem:global-all-source-face | paper/sections/global.tex:479 | G5, B1 | foundation §2 G5; sixth rebuild §3 B1; source paper/sections/global.tex:479-496 |
| lem:global-full-graph | paper/sections/global.tex:550 | G6, B2 | foundation §3 G6; sixth rebuild §3-4 B2; source paper/sections/global.tex:550-569 |
| lem:global-dyadic-small | paper/sections/global.tex:614 | G7, B3 | foundation §3 G7; sixth rebuild §3-4 B3; source paper/sections/global.tex:614-623 |
| lem:global-free-face-strips | paper/sections/global.tex:692 | G8, B3 | foundation §3 G8; sixth rebuild §3-4 B3; source paper/sections/global.tex:692-706 |
| lem:global-small-source-recovery | paper/sections/global.tex:741 | G8, B3 | foundation §3 G8; sixth rebuild §3-4 B3; source paper/sections/global.tex:741-753 |
| lem:global-common-defect | paper/sections/global.tex:821 | G8, B3 | foundation §3 G8; sixth rebuild §3-4 B3; source paper/sections/global.tex:821-837 |
| lem:global-augmented-capacity | paper/sections/global.tex:862 | G9 | foundation §3 G9; source paper/sections/global.tex:862-876 |
| lem:global-exterior-union | paper/sections/global.tex:991 | G10 | foundation §4 G10; source paper/sections/global.tex:991-1007 |
| lem:global-interior-recovery | paper/sections/global.tex:1121 | G11 | foundation §4 G11; source paper/sections/global.tex:1121-1133 |
| prop:global-coarse | paper/sections/global.tex:1226 | G12 | foundation §4 G12; source paper/sections/global.tex:1226-1229 |
| sharp-critical | paper/sections/sharp.tex:77 | S1, S2 | sixth rebuild §4 S1-S2; source paper/sections/sharp.tex:77-112 |
| sharp-input-tiling | paper/sections/sharp.tex:247 | S3 | sixth rebuild §5 S3; source paper/sections/sharp.tex:247-259 |
| sharp-angle-integrability | paper/sections/sharp.tex:315 | S4 | sixth rebuild §6 S4; source paper/sections/sharp.tex:315-325 |
| sharp-incidence | paper/sections/sharp.tex:363 | S5 | sixth rebuild §7 S5; source paper/sections/sharp.tex:363-368 |
| sharp-fresh-cdf | paper/sections/sharp.tex:453 | S6 | sixth rebuild §8 S6; source paper/sections/sharp.tex:453-461 |
| sharp-source-price | paper/sections/sharp.tex:483 | S7, S8, S9 | sixth rebuild §§8-10 S7-S9; source paper/sections/sharp.tex:483-496 |
| sharp-upper | paper/sections/sharp.tex:646 | S9, S10, S11 | sixth rebuild §§10-12 S9-S11; source paper/sections/sharp.tex:646-652 |
| eighth:theorem | paper/sections/eighth.tex:9 | eighth:theorem, eighth:exact-expansion | reviews/063-eighth-reconstruction.md §§2.1, 9-10; node table §4; dependency_graph.mmd (node table §4 lines 160-181); source paper/sections/eighth.tex:9-25 |
| eighth:priors | paper/sections/eighth.tex:43 | eighth:priors, eighth:positive-prior, eighth:mass-prior | reviews/063-eighth-reconstruction.md §2.1 steps 1-3; §4 table; foundation §4 G10-G12; sixth §§4,11,14 (node table §4 lines 160-181); source paper/sections/eighth.tex:43-57 |
| eighth:window | paper/sections/eighth.tex:133 | eighth:window, eighth:void, eighth:root-transfer | reviews/063-eighth-reconstruction.md §4.1 and §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:133-146 |
| eighth:soft-capacity | paper/sections/eighth.tex:188 | eighth:soft-capacity, eighth:soft-weight | reviews/063-eighth-reconstruction.md §4.1 (node table §4 lines 160-181); source paper/sections/eighth.tex:188-193 |
| eighth:good-slots | paper/sections/eighth.tex:333 | eighth:good-slots, eighth:bad-root, eighth:good-high | reviews/063-eighth-reconstruction.md §4.2 and §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:333-342 |
| eighth:high-net | paper/sections/eighth.tex:417 | eighth:high-net, eighth:high-payment | reviews/063-eighth-reconstruction.md §4.2 and §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:417-423 |
| eighth:joint-price | paper/sections/eighth.tex:642 | eighth:joint-price, eighth:joint | reviews/063-eighth-reconstruction.md §4.2 node table; detailed derivation §§5,7 (node table §4 lines 160-181); source paper/sections/eighth.tex:642-654 |
| eighth:bootstrap-upper | paper/sections/eighth.tex:784 | eighth:bootstrap-upper | reviews/063-eighth-reconstruction.md §5; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:784-786 |
| eighth:exact-stability | paper/sections/eighth.tex:926 | eighth:exact-stability, eighth:exact-stable-budget | reviews/063-eighth-reconstruction.md §5; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:926-944 |
| eighth:exact-high-gross | paper/sections/eighth.tex:989 | eighth:exact-high-gross | reviews/063-eighth-reconstruction.md §4.2 tail of derivation; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:989-995 |
| eighth:exact-axial | paper/sections/eighth.tex:1046 | eighth:exact-axial, eighth:exact-height-nec | reviews/063-eighth-reconstruction.md §6; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1046-1062 |
| eighth:exact-iteration | paper/sections/eighth.tex:1146 | eighth:exact-iteration, eighth:exact-second-cost | reviews/063-eighth-reconstruction.md §6; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1146-1163 |
| eighth:exact-normal-price | paper/sections/eighth.tex:1221 | eighth:exact-normal-price, eighth:exact-normal-integral, eighth:exact-model-price | reviews/063-eighth-reconstruction.md §7; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1221-1242 |
| eighth:exact-complement | paper/sections/eighth.tex:1395 | eighth:exact-complement, eighth:exact-complement-price | reviews/063-eighth-reconstruction.md §7; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1395-1408 |
| eighth:exact-small-box | paper/sections/eighth.tex:1553 | eighth:exact-small-box | reviews/063-eighth-reconstruction.md §8; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1553-1568 |
| eighth:exact-small-source | paper/sections/eighth.tex:1608 | eighth:exact-small-source, eighth:exact-inside | reviews/063-eighth-reconstruction.md §8; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1608-1617 |
| eighth:upper | paper/sections/eighth.tex:1689 | eighth:upper | reviews/063-eighth-reconstruction.md §9; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1689-1695 |
| eighth:macro-lower | paper/sections/eighth.tex:1811 | eighth:macro-lower, eighth:finite-lower | reviews/063-eighth-reconstruction.md §10; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:1811-1817 |
| eighth:lower | paper/sections/eighth.tex:2035 | eighth:lower | reviews/063-eighth-reconstruction.md §10; §4 table (node table §4 lines 160-181); source paper/sections/eighth.tex:2035-2041 |

基础语义图和参数次序见基础重建nodes（其内容已在报告逐项展开）；六阶S1–S11见六阶报告§1。八阶G0–G4→E0–E13和独立L0–L4下界见顶层报告及[图](063-top-level.mmd)。当前图终点是N8精确八阶，不自动供应新N9、真实任意阶或模型外的极值。
