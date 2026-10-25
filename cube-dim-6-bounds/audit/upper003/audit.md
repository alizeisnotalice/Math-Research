# upper003 direct-kernel 修补：独立审计

## 结论

**部分接受，且严格限于诊断。** 新的 `direct_kernel_laurent.cpp` 对零和重复
节点的 *形式 Laurent jet* 写法正确；旧 `direct_kernel.cpp` 对含零/重复节点
的删项输出仍为 `STRUCTURAL_FAILURE`，完全不能使用。新的 C++ Laurent 输出是
53-bit 浮点诊断，而不是 principal-part 精确相消证书，不能用于 (G) 的严格
点值、(g)、face 积分、成本或 (C_6)。

独立 exact fmpq 单点重建确认：对一个已登记的全 generic 有理点，direct
formula 的整数节点分类和 (operatorname{DD}_5(operatorname{sgn})) 结果
正确。它只校准该一个点的数值接口；不验证其余 31 点，也不支付有限点到连续
face/cube 的缺口。

## 数学检查：four-vector Laurent recurrence

新程序使用

\[
(\pi,\epsilon^{-1},\epsilon^{-3},\epsilon^{-5})
\]

四个坐标。对于 (F_\epsilon(b)=\arctan(b/\epsilon))，其 `jet` 恰为：

\[
 F_\epsilon(b)\to \tfrac\pi2\operatorname{sgn}b\quad(b\ne0),\qquad
 F_\epsilon'(0)=\epsilon^{-1},\quad
 {F_\epsilon^{(3)}(0)\over3!}=-{1\over3\epsilon^3},\quad
 {F_\epsilon^{(5)}(0)\over5!}={1\over5\epsilon^5}.
\]

All other displayed limiting jet coordinates are zero. Thus the source's
`p=+/-0.5`, `e1=1`, `e3=-1/3`, and `e5=1/5` are correctly normalized. After
exact sorting, `eq(x[i],x[i-o])` is exactly the confluent Newton condition:
equality means the intervening sorted block is the same node, so substituting
the order-(o) derivative jet is legitimate. This avoids the invalid former
operation of discarding a cone term just because it contains a zero or a
repeat.

There is an essential limit qualification: at a nonzero node the omitted
derivatives and the nonconstant part of (F_\epsilon) are (O(\epsilon)).
The four-vector program therefore describes a Laurent *limit* of the complete
finite cone sum; it is not a finite-(\epsilon) evaluator and cannot itself
prove a tail or a uniform bound.

## Actual independent checks

1. I compiled and ran `direct_kernel_laurent.cpp` afresh on all four v1
   points. Every non-timing field of each JSON row exactly equals the saved
   `laurent_v1_results.json`; elapsed times differ normally. The fresh run
   took 0.0682, 0.0680, 0.0681, and 0.0686 seconds per point.
2. I independently reconstructed the already-fixed first high-denominator
   generic point with python-flint `fmpq`, using all 21 base cells and all
   (46080) signed permutations: exactly 967,680 terms, no imports from the
   producer's DD/aggregation code, and no term skipping. It produced zero
   zero-nodes and zero repeated pairs, and exact four-jet negative coordinates
   (0,0,0). Its 512-bit enclosure is

   \[
   [6.3984213567624446970612568299985983373603382701502085420992921376686675008584333
   \cdot10^{-5}\ \mathbin{\pm}\ 1.50\cdot10^{-85}].
   \]

   The C++ generic diagnostic at that point is
   (6.39842135673982659186\cdot10^{-5}), with actual relative discrepancy
   (3.5349508793\cdot10^{-12}). This is a calibration observation only,
   not a rounding enclosure.

## Integer and floating-point boundary

For the independently rebuilt generic point, every raw C++ node numerator
has absolute value at most (8{,}327{,}839{,}898); the largest multiplication
summand is (5{,}000{,}000{,}035), and the largest absolute l1 bound before
signed accumulation is also (8{,}327{,}839{,}898). These are below
`LLONG_MAX`. The largest raw denominator is (160{,}000{,}001{,}120), and
the largest actual exact-comparison product — the source's `n*d` in
`__int128` — is (1{,}202{,}548{,}648{,}737{,}840{,}482{,}240), far below
(2^{127}-1). The smallest adjacent exact node gap encountered is
(667/36{,}000{,}000{,}252), so the generic point is genuinely nondegenerate
under the source's integer comparison.

On this host, compilation reports `LDBL_MANT_DIG=53` and
`sizeof(long double)=8`: source `long double` is IEEE-double precision, despite
printing 21 digits. Kahan summation cannot bound the divided-difference
subtractions or the final cancellation. The generic cancellation ratio is
(5.5977\cdot10^5), while the 32-point report records a maximum
(1.4048\cdot10^6). Hence these long-double values are valid only as
short-feedback numerical diagnostics. In particular, the C++ Laurent output's
displayed zero negative-power sums do **not** establish exact cancellation.

## Explicitly excluded evidence

The original denominator-257 first four outputs are labelled
`STRUCTURAL_FAILURE`; they skipped cone terms with zero/repeated nodes. They
are excluded from this audit and from every numerical comparison. The repaired
Laurent program does not skip them, but its floating sums remain nonrigorous.
Neither the 32-point sample mean, its displayed formal cost, nor a sampled
minimum can be called (I), (g), a cost certificate, or a (C_6) bound.

## Reproduction and hashes

Independent exact command (actual 44.18 seconds, exit 0):

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  /Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/upper003/independent_generic_fourjet_fmpq.py
```

| artifact | SHA-256 |
|---|---|
| repaired Laurent source | `f6713a81144f17f7e636706b70d30048f21e3b47a7485c29d33e6b05c52b54d4` |
| locked integer fan | `7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812` |
| fixed generic points | `c6505d85ad7566cbf24886ec0b7dac932909bfa7fd9863615361ea3579b0736f` |
| generic C++ diagnostic | `ac6b4a7567edbe6da6bf5fc63f1272d3bd94a70e328de4f8b6c6feaff29e7e95` |
| independent fmpq program | `17920407aeaefc235e6433d0c51e92c637f4cdc9c7a076daa6b37796ae1ea32c` |
| independent fmpq output | `4b12a936a4296943b1d4a3e770c7166b7b94399e1b54fa8905352904e9a03bb0` |

主集成将输入路径改为项目相对位置；本次主clone复跑记录见同目录结果与rounds/003运行说明，集成脚本hash与独立原脚本hash分别记录。
