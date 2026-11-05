# R004 便携发布去重补签

2026-11-06。**PASS：仅元数据与同字节日志去重，没有数学或执行变化。** 本补签保留此前portable-review.md及portable-review-binding.json原字节，不替换原37文件签审记录。

已将当前scratch/round004-public-stage/round004/{weights,jet}逐文件对照原37文件inventory：当前36文件，只有三处差异——删除weights/portable-transfer-finite.log，修改weights/README.md的共用stdout说明，修改weights/SHA256SUMS的README哈希并删去冗余日志行。无新增文件，其余34文件哈希全部相同。

Terra源包中的transfer_finite.log与portable-transfer-finite.log均为397字节，逐字相同，SHA256均为9b970bc3ba5ac576b3225d1c820f405c9d6ddde8a468e140e359f58ba62b2fe1；发布保留的weights/transfer_finite.log也与二者逐字相同。一次保留该相同stdout不丢失任何数值输出。两个原源文件仍在本地未修改，原运行身份由已有receipt及root的两源到一发布文件绑定保留。

weights README准确说明原始与移位有限窗stdout同字节合并为transfer_finite.log，移位WL输出仍单独保留。SHA256SUMS除该README新身份与冗余行删除外未变，当前全部行均重新按字节核实。没有执行源码、输入、数值、其他日志或receipt变化，没有运行数值程序，也没有重审数学。

最终36文件inventory、三项精确差异与旧签审身份另见portable-cleanup-checks.json。root后续只需添加本补签、总source-bindings与manifest。
