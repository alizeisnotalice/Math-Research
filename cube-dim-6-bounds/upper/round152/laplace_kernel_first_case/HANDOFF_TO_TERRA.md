# C152首例交接（Sol未执行）

round152起初没有既有目录。这里只完成独立正Gamma推导、一般
二维接口与未执行代码；没有运行prepare、预检、Arb或新数值。
下界150保持STATUS_USER_PAUSED，未作任何继续工作。

请把本目录源文件复制到Terra独立clone相应upper/round152目录。
唯一首例是147真实cone0/pair(16,17)，h1/128、lambda30；不生成
或重复149的density。先由Terra运行prepare_for_terra.py，冻结
terminal147原型/next证据、独审和130源文件，生成登记/源hash。
源metadata准备有准确有理检查，Sol没有执行它。

主及独审批准后，parent_release.json须含round="152"、authorized、
main_approved、independent_review_approved均true，以及登记中的
case_sha256、proof_sha256、general_interface_sha256、frozen_sources_sha256
和registration.json自身的registration_sha256。首次执行前不得
有started/certificate。单worker、nice10、线程1；四个环境变量
OMP_NUM_THREADS、OPENBLAS_NUM_THREADS、MKL_NUM_THREADS、
VECLIB_MAXIMUM_THREADS均为1。无时间/CPU cap，不自动重试。

只运行single_case.py一次。它用Arb512整体包围I30表达式，给
900I30平方和(46080/pi^3)*alpha0的本pair贡献。验收为严格正、
概率上界及固定10^-80宽度门。0逐点k、0quadrature、0旧6D网格、
0MC/native、0density再生成。保存完整child exit/reap墙钟/CPU、
stdout/stderr及certificate。代码的lambda0/2分支仅静审不另跑。

结果仅一个真实正贡献，不是完整beta2上端，不更新M/C6。一般
9961模式的后续signed密度积分与全域1e-9误差建议已写在接口文档；
是否可执行、选哪一个非可分离案例，等待149完整metadata和本次
实际反馈再由主登记，不能自动扩大。
