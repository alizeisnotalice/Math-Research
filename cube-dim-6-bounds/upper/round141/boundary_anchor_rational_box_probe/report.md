# C129终态残差与C141有限代数准备

C129已正常终止：detached receipt为exit0，当前controller计数新started=reaped=18257、active=0、live owned process表为空；累计C127+C129 reservations为21384，其中包含C127的3127。全325根终态为290根closed、35根 `NATIVE_BUDGET_EXHAUSTED`。续接55根中终态20根closed，但6根在C129开始前已闭，本次新增闭根准确为14，全森林从276→290；不能把20误写成新增。35个残差根均用完其原累计512 native预算，停止原因是预算，不是中心/真实G的反例。

只读终态树共40537节点、20106 split节点和20431终端叶；终端叶分为10165 passed、10231 queued、35 blocked，无empty叶。剩余10266叶中10262叶没有已存prefix，只有4叶保存了一个band。已核根状态与queue/blocked、叶分区及closed-coverage列表、child-parent引用、终态checkpoint计数一致；终态index/journal/counter之间hash链也一致。完整全域覆盖仍为否，尤其不能从290 closed根推出global g或新C6。这里没有重算旧核/旧完整盒证书；已存数学证书的独立重审由另一Sol承担，本轮结构核对不替代它。

控制逻辑给出比“继续加空间节点”更明确的改进依据：20091次非恢复几何split都有已存prefix中心lower高于g0=1/65000，而完整box lower仍低于g0；其中18252次只用了首band、1835次用了前两band、4次用了前三band。此前控制器一旦中心过g0便先尝试完整盒界，盒界不过就split，未先尝试追加下一个正时间带；另有962次split的已存中心lower尚未达到更强gstar。这支持以后由“完整box证书缺口”触发补已注册短band，而非只按中心阈值决定继续时间窗、随后重复盲分空间。它只是机制诊断，不预言某个未知短band的收益。

R137在真实boundary anchor c=(1,1,1,3/8,1/4,1/8)的五个正带jet均已128位accepted。本准备已先核原prepared pins及五份raw/receipt/jet/intent/process/stderr hash链，五semantic不会重发。首band中心值仅以约7.0164e-9余量超过gstar，而五band中心累计约2.693248e-5；点值改善本身并不是完整box证明。C141因此先只复用这五jet，对终态原未决叶作零native、零exp的完整盒代数探测，优先检验现存数据在空间证书上能否兑现短band收益。

固定选择全部10266残差叶中离该anchor最坏frame32二次距离最小的128叶，以root/node ID破同分，来自3根：C121.r11011110110的63叶、C121.r11011101111的57叶及C121.r11011110111的8叶。选择是纯几何排序，未预先计算新盒下界。固定odd集合为{1,3,5,7,9,11,13,15}；每个原robust角值f exp(r)用f·max(0,P1(r),…,P15(r))作精确Fraction下包。32角规则来自原robust logT的凹性/原T的角点控制，绝不假定Podd凹；每band先min32角，再对五band求和，每band用各自a。零下端安全但不能称严格正下端。详证见PROOF.md，独立理论审阅由142承担。

C141预算固定为128盒×5band×32角=20480 case，共享Taylor递推至15阶最多307200次精确更新；单CPU、nice≥10、无时间或CPU上限。g0与gstar分别验收，lower未过只记未决，真实g反例/实际G上界计数均为0，因为本法不计算用于反驳的upper。没有新源、native、exp、checker或细分；任何hash/schema/预算不一致即停，不自动重跑。当前保持UNRUN等待主与独审gate。若局部gstar通过，也不能把原来仅g0-paid的叶自动提升为global gstar；若代数法未决，未来再由明确盒缺口注册有限短窗/中心预算，不从这一弱下界失败直接判定真实g或原指数证书家族失败。
