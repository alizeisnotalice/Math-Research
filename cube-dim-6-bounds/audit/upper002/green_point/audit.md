# upper002 green point：primary Laurent 四 jet 独立审计

## 结论

**通过（一个点，条件性）**：在锁定 21-cell fan 及已写明的有限-cone 截时
DD5 恒等式条件下，主点

\[
x_*=(1,3/4,3/4,3/4,3/4,0)
\]

的完整 signed-permutation cone 和的三个 Laurent 负幂系数以准确有理数全为零。
因此该有限和的 Green 极限存在，且

\[
G(x_*)=
\frac{142703764082773109166538358568744627214353397671239027600712621842743987660709129826323249}
{13362678374297216728025812287452466802901590672721726385717434629016827885771180000000000\,\pi^5}.
\]

所以它只给 cube infimum 的一个**点值上界**：\(g\le G(x_*)\)。它不是
\(g\) 的下界或等式，且并未认证 \(I\)、模型成本或 \(C_6\)。

## 独立复算

实际命令（2026-10-25，本地，exit 0，1.8 秒）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  /Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/upper002/green_point/independent_laurent_primary_audit.py
```

该程序不 import `laurent_points.py`，不调用其 DD 或 aggregation。它以一个独立的
四分量 Laurent jet
\((\epsilon^{-5},\epsilon^{-3},\epsilon^{-1},\pi)\) 进行 confluent Newton
递推，并以 `Counter` 实际枚举全部 46,080 个 signed-permutation actions。

| check | exact result |
|---|---:|
| distinct orbit / stabilizer | 960 / 48 |
| all actions | 46,080 |
| 21-cell direct pre-aggregation terms | 20,160 |
| exact aggregated node tuples / tuples containing zero | 20,076 / 4,408 |
| \(\epsilon^{-5},\epsilon^{-3},\epsilon^{-1}\) coefficients | 0, 0, 0 |

The recomputed 512-bit Arb enclosure is

\[
[3.48973580695852672672944665488126173233092438526283402746678659700462196148970765822527279
\cdot10^{-5}\ \mathbin{\pm}\ 4.51\cdot10^{-95}],
\]

and its difference from the producer's saved 384-bit ball contains zero.

## Why the cancellation implies a finite limit

For (F_\epsilon(b)=\arctan(b/\epsilon)), at a zero node the only nonzero
jets relevant through DD5 are exactly

\[
F'(0)=\epsilon^{-1},\quad F^{(3)}(0)/3!=-1/(3\epsilon^3),\quad
F^{(5)}(0)/5!=1/(5\epsilon^5).
\]

At fixed nonzero nodes, the zeroth jet is
\((\pi/2)\operatorname{sgn}(b)+O(\epsilon)\), and every positive-order jet
is \(O(\epsilon)\). The exact finite sum has all three negative powers
cancelled before the limit is taken. Its remaining discarded part is therefore
a finite (O(\epsilon)) sum. No zero-node principal part is deleted cone by
cone.

The base 21 coefficients total `full_exp_minus_symbol_integral / 46080`; exact
orbit multiplicities restore the full JSON normalization. This is checked in
the independent program.

## Scope limit

This does not reprove that the 21 base cells are a complete fan, nor the
Fourier/DD5 kernel interface; it uses both as locked inputs. It does not
control any neighboring point, the cube minimum, the face integral, a tail
uniformly in a region, or a weak-type bound. The exclusion-gate document is
correct that a pointwise (g\)-upper bound alone cannot exclude a model cost:
a face lower-bound mechanism remains necessary.

## Hashes

| artifact | SHA-256 |
|---|---|
| locked fan JSON | `18d8aa578f84283d00e1c9c30522fb14cc49812092afd5716f2cc9ae68ba2798` |
| producer source | `ebf204a8ae1da93edb2ed785ea0bf64e61029f38b7d6e9dc98a6c8e7865396c2` |
| producer results | `7a11ef2029b7d202d042c434b523e5e70bd0a0b43a1a3cf703c413cd6f74d884` |
| independent script | `2151f9b010b2cc8c000196701622b9ce9159b23a10b37133364cabb9a0ad3de7` |

主clone将来源路径改为项目相对路径并实际复跑退出0。集成版源码SHA及结果见同目录*_results.json；上述原SHA保留作为独立代理原始运行来源。
