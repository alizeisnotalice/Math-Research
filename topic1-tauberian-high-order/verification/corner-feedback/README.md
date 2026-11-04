# R26 computation archive: corner feedback

This minimal package retains the symbolic derivative/integral checks and the finite rational arithmetic checks only. The separate `boundary_checks.py/json/log` files are excluded; those checks are optional auxiliary evidence and are not a dependency of the retained scripts or the reviewed argument.

## Symbolic scope

`corner_symbolic.wl/log` records a Wolfram 15.0 run. For the capacity term
\[
 a((R+s)^n-R^n)-[C(s-h)-m]_+,
\]
the kink is at `h+m/C`. Before the kink the derivative is `n a (R+s)^(n-1)>0`. After the kink it is `n a (R+s)^(n-1)-C`; the claimed strict decrease on `0<=s<=S` uses the explicit condition
\[
 C>n a(R+S)^{n-1},
\]
which makes that derivative strictly negative throughout the post-kink portion. The general-n symbolic output leaves this needed inequality in the derivative condition; monotonicity of the power and the stated condition supply the implication. Thus the maximum is at the kink if it lies in `[0,S]`, and otherwise at `S` (with the given nonnegative parameters).

The two-arm integral identities and the recovery coefficient checks are symbolic algebraic checks for the supplied bounds. They do not verify geometric inclusion of the arms or disjointness/counting of source layers.

## Finite checks and scope

`corner_finite.py/json/log` checks the supplied analytic upper/lower budget bounds by exact rational arithmetic for `n=2,3`, `u=2^-30`, and `C=10^6`. The checks establish positivity of the computed bound slack and listed guards under those inputs. In particular, these are not actual measurements of the underlying sets: the `n=2,3` budget checks apply the given analytic estimates, rather than computing the measure of the actual bad set. The package therefore records a candidate one-sided lower construction only; it has not matched a corresponding upper bound and does not establish the true next asymptotic coefficient.

## Reproduction

From this directory, run:

```sh
wolframscript -file corner_symbolic.wl
python3 corner_finite.py --j 30 --C 1000000
```

The Python script writes `corner_finite.json` in the current directory, so reproduce in a disposable copy or preserve the archived output first. The Wolfram command emits the symbolic checks to standard output; compare with `corner_symbolic.log`.

## Whitelist and SHA-256

- `corner_symbolic.wl` — 2332 bytes — SHA-256 `45fd4c26c719fce17bbb4a0e09443002aca1a64982ce6e6a5dc945a351b966de`
- `corner_symbolic.log` — 877 bytes — SHA-256 `2c55427bf01d78d60a1806dcc9acaac44ac6652595d6a559ea66ca909d244f54`
- `corner_finite.py` — 2114 bytes — SHA-256 `069ac80df882c35ad95c1c9ec6a47db664a2f677aac5243603d579ed4a29c83a`
- `corner_finite.json` — 2346 bytes — SHA-256 `a2ab039588a4aec7e32e6df089a560c05a2483368fa0402530bc8a4e85c382a3`
- `corner_finite.log` — 93 bytes — SHA-256 `45573fb09b10924cc2ac8015ee6db7b50c6f18fef3f17962ce5d69996ce6d25d`
