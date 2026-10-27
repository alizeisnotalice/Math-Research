# 168：158复用source与完整域best的独立预审

本目录只准备不同整数路径和来源门，不执行新corner/q/log。165与157等原冻结件不改。158真实输出已按原math manifest冻结；实际待独审范围是19661个新16角点case（314576角点）、1833条first-PASS后的skip，合计固定21494对、1240 receivers、4769 selected sources。准备计数不等于168已取得通过结论。

## 复用资格与允许负margin110

每个source的EPS64/EPS128、point、q153、Delta115和四个physical directed loggradient必须逐字段回绑已独审157的155原q event。168不调用153、155或158的数学代码，不计算q，不再证明或求值公共log，只匹配原合格log文件SHA及157已通过的upper字段。新q可以小于256，Delta115可为负，anchor可在receiver外。

每个新source仍只用一个anchor、一个epsilon和整组梯度测试全16闭角点。165式正margin断言不适用于此任务：新的source对当前receiver可能给负margin110，仍是合法但较弱的整盒下包。168整数helper没有`margin110>=old>0`限制，负值照常比对并计入已测集合。

令Ulog为已资格155公共上端。这里恢复110必须逐source成立

`saved_Delta115 + Ulog = saved_Delta110`。

于是同一source的两条曲线严格只差Ulog，故`m110=m115+Ulog`、最坏角点相同。这不是对任意115下界加一个log上端来反推110；若只有任意`m115≤log(H/g115)`，加Ulog（大于真实log比）一般不能推出110下界。这里安全性来自已经保存并绑定的Delta恒等式。

168先清除整证所有项至一个共同正整数分母，再用16个整数分子确定最小值和字典次序最坏角点，区别于生产Fraction逐项求和；不导入157带旧margin正性断言的helper。m115严格正才PASS115；其余包括零与负都只是该证未通过，不能当作G或C5反例。

## best、first-PASS与失败source

每receiver以原155正margin110作为fallback。每条实际evaluated完整单source证书都参与

`best110=max(old155110, all actually tested whole-box newsource margin110)`。

失败115或负110的source也必须进入这个max流程。只比较整盒最小值，不能先对各corner切换source再声称整盒通过。相等时保留原先证书，与生产顺序一致；first-PASS115后可以停止此receiver，余下预备对只能记skip。best只覆盖实际测试集合，不声称达到所有预备source的最大值。

全部固定pair必须精确消费一次，receiver按child_index%15保持单worker；hash-chain/seq/checkpoint/receiver_best文件及每条skip的既有first-PASS均独立核验。每次更新best即使115仍失败也保留改善，因此不允许只处理PASS标签。未重新运行KDTree或浮点候选排序：158的float只负责已冻结选择清单，数学结论不依赖候选最优。

## 原完整域与新全局最小值

只替换1240个旧pending receivers，其他109961条155已过115叶保持原margin与来源。168把155压缩margin投影逐项匹配其已资格原leaf events，并核receiver的原149闭box/node/index；新source旧155已测同box失败的exact-newcore pair必须被排除，不能重放旧case。

每个best110≥原155正值，全部111201叶和原52360根闭覆盖不变。各source只使用H64或H128，因H128≥H64，整盒best可统一为H128在F上的下界。选取max允许不同整盒source，绝不混梯度或corner；在完整旧域上再取准确min，核witness/source种类及pending115。目标115的partial状态与完整110域严格分开。

后续若主用正m做P8/down192接纳，应另登记；168只核生产保存的linear候选与原完整域。固定epsilon的H不向Q作齐次外推，只有G≥H128后沿029的G齐次性处理。成本表达式区间下端不能作为C5下界。本轮不计算P8、新g或Cn。

## 有限执行合同

`sources168.py`准备核4769 source全部alias原155 event/body/raw-record hash、Delta110 identity来源、已资格log、固定pairs、source index、1240 receiver闭几何与全111201原margin投影。`integer168.py`只做新source/box的整数角点与best操作。`audit168.py`准备15显式worker，仅对真实evaluated事件调用角点helper；skip只核元数据。Q核心原字节复用已审整数/gcd路径。

实际执行必须在准确plan、真实158完整输出binding经主审授权之后开始，先检查全部源，再做任何新角点。输入SHA变化、非法skip/漏pair/重复receiver、source复制不一致、best遗漏失败case、min/witness/pending不符都保存失败并停止，不自动重试、改source或增加case。0新q/log/H/gradient、0旧corner/log基线，原log仅继承hash与已经通过的证书，不重新eval。

当前实际行为仅为源码/schema读取、byte hash及AST解析；0新角点。旧165/157/152/154冻结件不改，无Astra、MATLAB、下界、提交或推送。主保留最终接纳责任。
