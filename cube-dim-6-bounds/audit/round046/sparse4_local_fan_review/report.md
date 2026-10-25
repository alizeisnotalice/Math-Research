# F046 sparse4 首父 cone 独立审计：PASS

A046 的首父 cone 局部 fan 可接受：32个闭pattern完整覆盖，12个full cells与20个非空低维区域全部保留，12条primitive rays形成14个full simplices。覆盖通过完整半空间顶点与facet/pulling归纳证明，独立检查没有只依赖体积和1。

真实新frame的 Gram 为(2561/80)I，六个坐标边际rate均为33；这是152个共同时间Cauchy参数的正保守齐次半群。首父cone被五个新form切开，旧q-only更新不合法。本审计仅接受这个父cone，未接受全fan或任何核值、tube迁移、I/g/Cn新界。

## 独立方法与实际结果

独立review.py不执行或导入作者writer/checker。唯一完整462 active系统枚举先消去u6，再用标准库Fraction做5×5 RREF，区别于作者6×6 FLINT求解。结果：237 singular、6 outside、12 distinct vertices，逐项与冻结fan相符。每个系统的输入活动集合、分类及非奇异解完整保存在results.json。

围绕这一次枚举仅做12个结构关系组：frame/Gram/边际、父cone与物理截面、首cone crossing、全部active系统、全部closed patterns、准确facets、bottom-up pulling、共享面、primitive rays/rates、152forms与线性slope、signed determinants与体积、global负号范围。

独立facet算法先取包含关系极大的proper zero faces，再核其维数，未使用作者的facet维数筛选作为唯一来源。全部66个proof DAG节点、12个根的递归准确吻合。全部66对共享face比较了所有极大诱导simplices，包含低维tie；结果一致。每个非full pattern都有经完整顶点核验的full-cell闭包包含见证。

首cone只重算720个new-form/ray整数dot，未重跑旧2160基线。12条新ray的真实152forms共1824个dot用于rates与符号；整体负ray只核同一偶symbol对称性，未生成新空间/时间点。14个子锥分别核全部152forms、真实slope与带符号primitive determinant，未沿用旧q。

## 体积修正与权重

参数u截面体积为√6/120；物理xi截面体积为√3/120。独立五条物理边的Gram determinant=3，直接验证最终标签；物理法向ell=(0,−1,1,1,0,0)，面积Jacobian=1/√2。相对体积仍为|det U_child|，14项准确和1。此截面不是psi=1截面。

原始记录与最终元数据hash已有区分。本审计验证最终对象及真实Jacobian，未将历史hash声明本身当作算术证明。primitive determinants的±号只表示orientation，积分使用绝对值；局部对象未添加global负号因子2。未来若合并相反锥，只能在全fan覆盖和共轭关系均确认后取2Re。

## 收据与范围

PID 49425，实际独立计算0.047149秒，源锁至完成数学核验130.523346秒，小于15分钟。一次462枚举、12结构组、0数学失败或重试、0 kernel/native/G/H/grid/MATLAB/作者writer/全fan/旧2160/扩招/push。初始只读查找B044不存在的theorem.md返回缺文件，已登记；随后读取其冻结report，未运行计算或writer。该文件查找不是数学失败或重试。

所有源hash在执行前后核对，B044完整白名单及A046旧输入均核对；不修改任何源。独立代码、完整active ledger、共享面证据、determinants、proof与收据只写本目录。

## 规范入口

- A046 manifest: e0b943a42323ddb0f9cf2c937141bb252ee17c93e55927fd3de343cba6f74ca7。
- A046 local_fan: bae26cd08d31dac644e5ae0e77483ae912de0f687a03d6bc3d7b98a4a35a78b6。
- A046 theorem: aebca68e55d6720191c91cea4a02efe0df8b58bdec363542d12177b22ba11251。
- B044 manifest: b1a705ea2f4a06052a837642cfde6469907729baccba1cb58fc175aad5ad6306。
- B044 frame: 67a12cd67961baf5ee554c138c82d6b4b627701c66dc691cb6f8b2a823c9be1c。

完整证明见proof_review.md，独立代码与实际范围见review.py、results.json、execution_receipt.json；本目录由manifest.json封存。其余20父cone及全域数值接口仍未证，本批停止，父代理保留证明与集成责任。
