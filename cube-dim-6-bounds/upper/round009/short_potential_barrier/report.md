# A009: a uniform positive short-potential barrier near all of S2

For the existing **alpha=1 body32 Cauchy** model and T=4/125, let
delta(x)=dist(x,S2), with Euclidean distance to the union of all 496 body
pair planes. The batch proves the simple whole-cube lower barrier

    K_(4/125)(x) >= (1/5000000000)*(1/delta(x)-40000)_+,

with K_(4/125)=+infinity on S2. In particular,

    ||x||_infinity<=1 and delta(x)<=1/80000
       => K_(4/125)(x)>=1/125000>1/128000=T/4096.

The exact excess of the simple lower bound over the required threshold is
3/16000000. This includes S1 and every pair-plane intersection. It is a
continuous-region certificate, not finite-point evidence or an I estimate.

The remaining compact set

    {x in [-1,1]^6:dist(x,S2)>=1/80000}

is explicitly unproved. No global m_T or new C6 upper bound is claimed.
B/F retain ownership of the short-potential source/weak-type transfer; that
interface is not re-proved in this subtask.

## Main proof and constants

The self-contained proof is in `theorem.md`. It uses the weaker but proved
Rademacher estimate

    8sqrt(6)||xi||<=psi(xi)<=32||xi||,

obtained from EX^2, EX^4<=8||xi||^4/3 and Hölder. No sharp Khintchine theorem
is assumed. For a selected pair v,w, its orthogonal R4 quotient has symbol
equal to the original psi restricted to span(v,w)^perp. The quadratic cosine
lower bound therefore gives

    q_t(y)>=a t^-4-b||y||^2t^-6,
    a=3/(4pi^2 32^4),
    b=15/[8pi^2(8sqrt(6))^6],
    b/a=5/108.

Deleting one additional variable leaves coercivity 7sqrt(6) and hence the
four-dimensional peak bound C_rest t^-4, where

    C_rest=3/[4pi^2(7sqrt(6))^4].

For the remaining 30-variable sum Y, removing the event ||P_UY||>L costs
at most

    (C_tail/L)t^-3,
    C_tail=(1800sqrt(6)/pi)C_rest.

This follows from a union of the independent events
|Z_j|>L/(30sqrt(6)) and their positive quotient-density convolutions.
It is a statement about pushforward submeasures and their densities, not
conditioning on QY=y. The initially suggested stronger t^-2/L^2 tail was
not needed or used.

Take kappa=sqrt(5/27) and tau=min(T,aL/(4C_tail)). On the positive time
window kappa dist(x,U)<=t<=tau, the Fourier second moment costs at most
a t^-4/4 and the truncation costs at most a t^-4/4. Thus at least
a t^-4/2 of the quotient density remains. This window is only a positive
lower bound for the full K_T; it does not replace the kernel or alter its
source identity.

The two pair parameters satisfy

    |z_i|<=sqrt(6/D_k)(sqrt(6)+L),
    D_k=det Gram(v,w)=4k(6-k).

Their change of variables has Jacobian sqrt(D_k), so the inverse factor
1/sqrt(D_k) is retained. Combining the two Cauchy lower factors with the
retained quotient mass gives

    p_t(x)>=d_k(L)t^-2,
    d_k(L)=a/[2pi^2 sqrt(D_k)(tau^2+B_k(L)^2)^2].

Positive coarea/Tonelli proves the density inequality almost everywhere;
fixed-time continuity of the true p_t then extends it to the whole closed
ball/tube intersection, including the pair plane itself. Integrating yields

    K_T(x)>=d(L)[1/(kappa dist(x,U))-1/tau],

where d(L)=min_k d_k(L). The time integral diverges on the plane. No
continuity of singular G or K_T across S2 is assumed.

## Exact geometry and preservation of intersections

Only five canonical geometries were evaluated with exact SymPy matrices:

|Hamming k|Pair count|Gram determinant|Smallest Gram eigenvalue|
|---|---:|---:|---:|
|1|80|20|2|
|2|160|32|4|
|3|160|36|6|
|4|80|32|4|
|5|16|20|2|

The counts sum to496. Every pair is covered by signed-permutation symmetry.
For T=4/125, d_k is increasing with the determinant throughout 20<=D<=36,
as shown by the exact derivative calculation in the proof. Hence k=1 or5
is worst, and only one scalar coefficient per L is needed.

At a point near the finite union of planes, choose a closest plane and apply
its lower bound. There is no factor496 loss and no summation of lower
contributions. The same argument applies to every overlap and S1 line.

For implementation, a plane distance can be computed exactly from

    dist(x,span(v,w))^2
      =||x||^2-[6(b1^2+b2^2)-2s b1b2]/D,
    b1=v.x, b2=w.x, s=v.w, D=36-s^2.

This is a usable geometric certificate for the tube/complement split; this
batch did not enumerate a new spatial grid or evaluate its complement.

## Eight actual scalar budgets

The registered physical L values balance tau proportional to L against the
pair-factor loss proportional to (sqrt(6)+L)^(-4). Each row uses the same
proved kappa, and all resulting tau values are strictly below T.

|L|Maximum distance rho_star for the sharp threshold formula|
|---|---|
|1/4|7.30998509312e-6|
|1/2|1.15597482002e-5|
|5/8|1.25278512754e-5|
|3/4|1.28703128518e-5|
|7/8|1.27364641602e-5|
|1|1.22713263946e-5|
|3/2|9.15624016728e-6|
|2|6.28467019229e-6|

The best **registered** choice is L=3/4, not a claim of globally optimal L.
Its outward Arb384 values are near

    tau =1.10128002484511221411e-5,
    d   =8.70465308705056072095e-11,
    c=d/kappa=2.02277858643958811773e-10,
    C=kappa/tau=39075.57326053003637384.

The strict comparisons c>1/5000000000 and C<40000 give the simple barrier
in the opening statement. At distance 1/80000 the unrounded formula has
the narrow Arb384 enclosure

    [8.278105407091553513260630791682380144836990978175321941536978469604070370682592421927267345682429390e-6 +/- 4.55e-106].

The direct rational lower 1655621/200000000000 is also verified and retained;
the simpler canonical value 1/125000 follows from the globally usable rounded
barrier. Both are strictly above 1/128000.

## Actual verification and stop status

`check_barrier.py` completed with exit0, Python PID **24029**, script time
**0.057840000 seconds** (tool chunk206abd). It checked five exact Gram
geometries, the two symbolic coefficient identities and exactly eight
parameter budgets. There were zero body32 kernel calls, zero MC samples,
no MATLAB, no grid, no subagents, no commits or pushes.

`verify_barrier.py` independently reconstructs the Fourier constants using
Gamma/sphere normalization, the Cauchy-tail union coefficient, all eight
budget values and the two rounded rational barrier comparisons. It adds
zero new parameter cases. The checked file records the exact distance,
K_T lower bound, strict margin, remaining compact set and source hashes.

The retained artifacts are this report, the self-contained theorem, finite
registration, two scripts, input manifest, complete results and checked
results. Original body32 model and F009 interface/point-diagnostic sources
were read-only and hashed. Source-interface proof ownership remains B/F.
F reran its context-only point-diagnostic output during this work; the
manifest preserves both the executed context hash and the later observed
hash. No A barrier constant depends on that point table.
The batch is complete and stopped; the singular neighborhood can now be
removed with a genuine lower certificate, while the compact complement
remains the next unresolved global task.

Reproduce from this directory using the existing environment:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_barrier.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_barrier.py
