# upper002 anchors：独立有限截时审计

## 验收结论

**接受（有限范围）**：在已锁定的 21 个 base fan simplex 以及理论文件给出的
六节点截时恒等式下，两个实际运行的锚点的 Arb 区间实现通过独立复算。该结论只
认证所列有限 (H_\epsilon) 值；它**不**认证完整 Fourier cone 覆盖、
\(\epsilon\downarrow0\) 的全局势、连续 face 下界、\(g_0\)、\(I_0\)，或任何
\(C_6\) 上界。

本审计使用的理论接口是

\[
 H_\epsilon(x)=(2\pi)^{-6}\sum_R\alpha_R
 \operatorname{DD}_5[\arctan(\cdot/\epsilon);b_{R,1},\ldots,b_{R,6}],
 \qquad \epsilon>0,
\]

即 `kernel/002_kernel_interface.md` 的 (A)。本次并未把它替换成
\(\epsilon=0\) 的符号函数公式。

## 实际独立运行

运行命令（2026-10-25，本地，exit 0，约 2.7 秒）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  /Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/upper002/anchors/independent_anchor_audit.py
```

运行环境：Python 3.15.0、python-flint 0.9.0、SymPy 1.14.0、Arb 512 bit。
审计代码没有 import 作者的 `truncated_anchor_batch.py`，也没有调用其
`aggregate_nodes` 或 `hermite_dd`。

| case | point | \(\epsilon\) | exact orbit / stabilizer | direct cone terms | merged node tuples | result |
|---|---|---:|---:|---:|---:|---|
| A | \((1,(3/4)^4,0)\) | \(1/1000\) | 960 / 48 | 20,160 | 20,076 | Reproduces the producer's non-pilot interval |
| B | \((1,(1/2)^2,0^3)\) | \(1/1500\) | 480 / 96 | 10,080 | 8,686 | Fresh-epsilon independent finite evaluation |

The direct-minus-merged Arb ball contains zero in both cases. For A, the
direct result is

\[
[3.48079723832451724433456244140170679164286169372536827982159274195282507457497633615859601
\cdot10^{-5}\ \mathbin{\pm}\ 3.05\cdot10^{-95}],
\]

and its difference from the producer's stored, display-rounded 70-digit ball
contains zero. For B, the fresh value is

\[
[0.000117420000785534475132830169970009858797661865665626884726522123324283638636643863193780833
\ \mathbin{\pm}\ 1.38\cdot10^{-94}].
\]

The complete enclosing balls (including both residual balls) are preserved in
`independent_anchor_audit_results.json`; no conclusion here is made from the
displayed decimal midpoints.

## What was checked

1. **Group orbit/stabilizer.** The code enumerates all (2^6 6!=46080)
   signed permutations into an exact `Counter`, rather than deriving a
   stabilizer from a set. Every image in each selected orbit has the stated
   uniform multiplicity and orbit-size times stabilizer is exactly 46,080.
2. **Fan normalization.** The 21 base-cell coefficients sum to exactly
   `full_exp_minus_symbol_integral / 46080`; after orbit multiplicities, the
   exact aggregate coefficient sum is precisely the JSON's full integral.
   This group factor is essential: equating the base sum itself to the full
   quantity would be incorrect.
3. **Hermite rule.** All nodes are exact `Fraction`s and sorted before a
   confluent Newton table is formed; hence a zero denominator occurs exactly
   when the contiguous Newton block is repeated. The code uses
   (F_\epsilon^{(k)}(b)/k!\) there. SymPy algebraically verified all five
   hard-coded derivatives of (F_\epsilon(b)=\arctan(b/\epsilon)).
4. **Real part/conjugation.** All (b_j=x\cdot r_j/q_j), coefficients, and
   group multiplicities are exact real rationals. The evaluator implements
   the stated *real* divided-difference identity for
   \(\int_\epsilon^\infty\Re\prod_j(t-i b_j)^{-1}dt\); it does not form a
   complex approximation and then discard unmatched imaginary terms.
5. **Truncation and rounding.** Both epsilons are positive exact rationals;
   the calculation remains a finite-ε integral to infinity. Equality,
   orbit multiplicity, sorting, and cancellation happen before conversion to
   Arb. Thereafter all operations, including `atan` and \(\pi\), use Arb balls
   at 512 bits. Thus the saved output is an enclosure under the installed
   Arb library's ball-arithmetic contract, not a binary floating-point
   sample.

## Qualification and remaining gaps

The producer's `verify_anchor_aggregation.py` imports the producer module and
tests only its first pilot anchor, so it is useful as a self-check but is not
an independent implementation. The present script is independent at the
formula/group/aggregation level, although it necessarily shares the locked
fan JSON. The latter was already subject to a separate fan audit; this anchor
audit does not repeat its (4096)-word/chamber certification.

Neither a finite collection of anchor intervals nor direct/merged agreement
controls untested points. In particular, it supplies no lower bound for
global \(g\), no face integration certificate, and no weak-type constant.

## Exact input and output hashes

| artifact | SHA-256 |
|---|---|
| locked fan JSON | `18d8aa578f84283d00e1c9c30522fb14cc49812092afd5716f2cc9ae68ba2798` |
| producer batch source | `6d5fe3fb515a7f31efa2b9af440170964c832f28d24c65a5040fd959274ee6e9` |
| producer results | `f5c91cfc7fdaae5b9650970085bd0664250eba380bf6c1e4b113d1617f27b6a8` |
| independent audit script | `7eeb848df39c46126c4ae18090a2e8a140bfedb7862602991a884ea247702b80` |
| independent output | `50f195e30ddd061e6fdfb93a18f0c50c0707af8867678074215ff48d0fba62e9` |

主clone将来源路径改为项目相对路径并实际复跑退出0。集成版源码SHA及结果见同目录*_results.json；上述原SHA保留作为独立代理原始运行来源。
