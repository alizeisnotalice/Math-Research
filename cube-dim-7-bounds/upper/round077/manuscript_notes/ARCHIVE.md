# Manuscript-only archive whitelist

Owned paths: paper/main.tex, paper/main.pdf, paper/BUILD.json, paper/MANIFEST.json, paper/README.md, and the explicit files listed in whitelist.json under upper/round077/manuscript_notes. Historical paper attachments and all unrelated untracked directories are excluded. No commit or push by the author.

Frozen baseline main.tex SHA256 8135adfb51470228ed18371c6c5f1a22caac0b35523b0cd49c50b13cfb9b386b, published commit f26db46db9e320aaffca74807ad25a815caad237, blob d4a09760b285ebc9fc39cbfd87d90aff2a7aed18; main.pdf SHA256 e0cc7e246ea8e8a3b0bc0f00873a9d5c857f4d802f082c58ac5a8c5fe75108a3, blob ad8abde8c6f3f159858a0d9f70bfd4e35fc1b666. Recovery was verified by reading the published commit blobs. The baseline files were not overwritten in the primary clone.

No duplicate old PDF or old TeX is archived. revision.patch and the published blobs provide exact recovery. Supplement contains only original exact candidate/certificate, the two independent exact review scripts and receipts, review and source identity, and the three local mathematical dependencies.

Actual author work:14 grouped Fraction transcription checks PASS; final two XeLaTeX passes; Ghostscript all37-page render and visual inspection. No search, optimization, production, shared MATLAB, lower-bound, n7, commit or push.

Independent full manuscript review remains the integrating agent's responsibility.

主整合去重：supplement_hashes.json 将13份本地作者对照副本逐一映射到主clone现存已跟踪原件的仓库相对路径。这13份字节副本排除归档白名单；主整合直接保留原R073作者与audit、R070二阶剖面、R064相对J、R066一阶剖面路径和准确hash。论文所指补充程序可从这些原路径读取；本稿的verify_transcription.py另为便携式新转写程序，原运行身份仍单独保存。

实际排版警告：3项Underfull hbox（正文1、文献2），全页图审可读；Overfull/undefined/missing均无。
