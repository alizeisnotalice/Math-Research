# 第43轮运行前独立静态审查

- 审读最终 `extend_piecewise.py` SHA `4bc087786a59772d8cc559aab411ada8446836a2609e67724ee0bf3d39f35c33`；PLAN SHA `d2780c37bcc5f0eadfab334dd899f1fa33420b08ba5f6d0dadb90c631ae19448`；wrapper SHA `54bccde8567bc1f97e5acf7de5e55ceccfd5fae26e42d1a937d0adfff4b72221`。结论：GO，无静态阻断。
- 旧362整叶及615/661双来源强制绑定最终模型；663追加已核四分片补证，初始363树叶、366证明区域、216 pending。历史partial仅元数据，不要求存在。
- 新盒仅一次AS和6中心Hint+2近角Hab；T1未先扣D、未乘c。新单侧修正为same-sign零、opposite-sign Muv·δuδv，AS仅扣原ES。4piece合计64行全正才接纳树节点；失败父piece不混入接受前沿。
- 仍以原GRID根深度为准，depth≤4，首depth4失败即停；1500新盒、12000 Hint/真实积分、6000piece/96000行及300秒外层hard明确。组件与调用账本分列，completed后guard如实处理，不把有限返回改称中断调用。
- 新保存piece使用已pin的round23纯Fraction限制代码；这不是独立代数核验。另编标准库幂基转换/代入/回Bernstein核验器，检查原树frontier及叶内piece两层分区，不运行新模型或积分。
- 最后checkpoint同步树叶及证明区域计数已局部复读。外部wrapper仅清理自有进程组，拒覆盖已有证据；当前审查未执行数学。
