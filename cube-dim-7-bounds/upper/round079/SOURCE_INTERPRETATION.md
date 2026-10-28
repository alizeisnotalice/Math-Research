# R079 原始登记字段解释

certificate_registration.json 基于 R073 登记模板生成，保留了不适用于本候选的 original_model_parameters、planned_gamma、planned_final_c、旧Cantelli retention与source提示文字。原文件的hash身份保持不变。这些字段不作为R079参数或理论依据。

实际执行的certificate.py只取 chosen_parameters，与certificate_amendment.json中的Ret、gamma、mL、price、c、slack；更强c=17/100由随后单独登记stronger_registration.json和stronger_guard.py验证。应按这三个具名输入和实际代码解析确认对象，不把模板文字当原件来源证明。

当前固定参数为L0=4096、s=3/2、qmin=2/3、qmax=89/100、sigma=11/100、chi=51/1000、D=1536，原误差容许量不变；最终待审gamma=7166/25000、price=5827/50000、c=17/100、slack=1/10000。

独立审核必须重建此tuple的全尾和算子桥。weighted_cap子目录是另一个局部引理，不属于本轮53项旧cap模型，不能用于追溯替代该次执行证据。R080另外登记和运行新模型。
