# McDiarmid upper及observer replay预审

结论：`round041_mcd/PROOF.md`和`round041_replay/upper_math.py`、`reconstruct.py`的数学/预运行接口通过，未发现必须修正项。父线随后完成同协议重建；本审已核对三stream hash及24个旧payload数组逐位一致，见 `replay_checks.json`。本人没有启动replay或MC。

先运行自有n3、8状态小例：48个Fraction正多项式坐标ratio检查、10个数值log-MGF、4个sum-tail检查通过。随后读完整MCD正文及生产probe；其6个case共84个hidden状态，独立两点law用于检验任意独立law的McD，未混用uniform-bin的variance/3。另以Fraction核对r0及C预算，均通过。

## 数学链

- 非空正multilinear项替换一个坐标的ratio由该坐标所有class核ratio的最大值控制；空count段、零prior项单独删去，不取log0。macro核每类均严格正，非空项log有限。
- 振幅v及真实均值下包l给min≥l−v，故max/min≤1/(1−v/l)。旧quarter-cap r=max(v/l)²/4给q=2√r，log坐标振幅≤−log(1−q)。这里必须保留**旧quarter r**；把它直接替换为r/3会少付因子。
- `Σc_i²≤[4(1+r0)/(1−2√r0)²]Σlog(1+r_i)`正确，分别用了−log界及r≤(1+r0)log(1+r)。r0=10^-6、C0=1000001/249001<81/20=4.05，固定C不按收益选择。
- source的r0预算约4.60163×10^-7<10^-6，第二个均值误差e及L≤2的4e均已支付。该预算仍以真实均值误差前提成立为条件；解析中点误差+经验guard小于10^-8不认证完整FP，正文已保留。
- Doob增量的conditional range长度为c_i，不是2c_i；自含Hoeffding证明给t²Σc_i²/8，优化后exp(−2gap²/S)。Jensen `E logF≤log EF≤log u`方向正确。
- λ>Σu时按λu_b/Σu分阈值，逐项union不需表间独立；Smax控制各S_b后可用R exp(−2log(λ/Σu)²/Smax)。空项和S0分别处理，λ≤Σu返回1保守合法，不能把负log gap平方。

## 实际upper实现

fixed0使用table0的旧logproduct；joint/prior使用各表最大旧logproduct乘C，R项union。每个L先形成同一L的表和均值，再取旧Markov/Cantelli、variance/3 Cantelli、McD的min，最后对233个共同尺度作union。旧Markov保底保留原source contribution cap，没有用lower active条件给upper设帽。

均值为0、variance为0或λ等于均值上包的分支均保守；strict equality处返回1上包虽然可松，但不制造错误小upper。`upper_math.micro_checks`及本审推导一致。

## Observer replay

reconstruct只包装原feature/bounds：feature读最大r、不改返回值；bounds先调用原函数，新数据通过stack/copy取得，不改stats或返回out；没有新增调用原Generator。BASE.HERE只改输出目录，导入时固定的SOURCE仍指40source。旧40文件不写入，新base_replay与side NPZ均放41目录。

最终必须精确核对count、sampled-anchor、source-statistics三stream hash，并逐数组比较全部旧payload；相同seed本身不是验收。side NPZ在完整验收前已经写出，后续新upper CI必须要求最终`RECONSTRUCTION_RECEIPT`验收通过及文件hash对应，不能仅因side文件存在就使用。

保存的logmean、old logvariance、max/table0 logproduct足以离线重算新upper。该过程是同一40数据的协议重放，不是新独立确认；新CI按1008预留预算、同原source/BE/完整FP限制报告。
