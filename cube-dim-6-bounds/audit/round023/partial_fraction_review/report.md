# F023 independent audit of the symbolic partial fractions and Green value

**PASS.** The partial-fraction identity, real fold, zero/repeated-pole
handling and local reduced-denominator argument are valid. The saved
density table has 359 positive pole parameters and true denominator degree
417 in z=T^2. Independently, one complete Hermite-DD5 Laurent pass through
all 40320 frozen tuples reconstructs exactly the same rational Gamma_G
as B023, with all negative epsilon powers cancelling exactly.

This review performed **one full symbolic kernel audit**, eight formal
scalar cases, and saved-table algebra. It performed zero native/numerical
kernel calls, zero positive-time evaluations, zero Sturm/root operations,
and no global numerator construction. It is not an all-time certificate,
a cube-wide lower bound, or a new C6 bound. No next-stage work is released
or started by this review.

## 1. Frozen common-time Fourier model and exact multiplicities

The physical model is the original body32 frame, at the new rational point
x0=(50,50,50,45,27,0)/50. For every stored six-node tuple B, the density
formula is

    p_T(x0)=(2pi)^-6 sum_B w_B Re product_(b in B)(T-i b)^-1.

This follows by the Laplace integral on a simplicial Fourier cone after
normalizing its ray coordinates by the exact positive rates. Its weight
is det(rays)/product rates. The density is a function of one common T;
there is no independent choice of time per coordinate or cone.

C022's source explicitly multiplies this cone weight by the point-orbit
multiplicity 12 before its first merge. Its subsequent global-negation
merge adds BOTH partner weights and counts. A negated tuple represents
the complex conjugate product with the same real part. Hence each frozen
w_B already includes the stabilizer and both global signs. Neither 12 nor
2 may be applied again. This review read that source logic without
re-enumerating its point orbit or old fan.

The independent pass uses each stored weight exactly once. It checks
40320 tuples, 80640 represented point-orbit/cone pairs, 967680 logical
full-group terms, and the full weight checksum

    520066141/48769479475200.

The checksum is accounting evidence, not a substitute for the kernel
identity. Zero-node tuples represent 88320 logical terms, repeated-node
tuples 201408; all are retained by the independent recurrence. The input
archive and decompressed content are pinned to

    562e3703286fc9a9712efdad30bbebb71433c4f3f15f867c985ee41b80c6f719
    7deb058449e11fa74d7e1a4658311cb738a38fbbb57ab03118b494138cd39390.

The original C022 structural baseline and B023 writer were not executed.

## 2. Complex local coefficients and the i phase

For a distinct signed node b of multiplicity m, set s=T-i b. Factoring
the other nodes gives

    product_(c!=b)(T-i c)^(-m_c)
      =i^(-(6-m)) product_(c!=b)(b-c+s/i)^(-m_c).

If a_k is the coefficient of u^k in the last product with u replacing
s/i, the pole term of order j=m-k has coefficient

    i^(j-6) a_k /(T-i b)^j.

The source computes a_k by exact negative-binomial series through degree
m-1, using distinct-node differences b-c only. Thus it does not divide
by zero at repeated nodes or at b=0. Matching every principal part and
vanishing at infinity proves the whole rational identity: the remainder
is a pole-free polynomial tending to zero. This also covers six identical
nodes. The actual input has maximum multiplicity four, but the formula
and formal tests retain the sixfold case.

All aggregation uses arbitrary-precision exact fmpq values. Inspection of
the source finds no positive-time substitution, binary64 coefficient,
zero-node deletion, or replacement of a confluent term by a generic one.

## 3. Real fold, zero poles, and genuine denominator degree

Let A_(b,j) be the complete signed aggregate before the i phase.
Conjugating the denominator inside the real part gives

    C_(b,j)=A_(b,j)+(-1)^j A_(-b,j), b>0.

The parity factor is necessary; the chosen global-negation representative
need not itself have a conjugate-symmetric complex table. At zero, odd j
is purely imaginary and has identically zero real part. The even real
channels are respectively A_(0,2)T^-2, -A_(0,4)T^-4, A_(0,6)T^-6.
All three coefficients are actually accumulated and found exactly zero.
The surviving nonzero A_(0,1) is harmless because it is purely imaginary;
its presence is not a discarded input node or an overwritten coefficient.

For b>0, z=T^2, the real numerator is

    N_(b,j)(z)=sum_(r=0)^floor(j/2)
       (-1)^(j-r-3) binom(j,2r) b^(j-2r) z^r,

over denominator (z+b^2)^j. At the pole,

    N_(b,j)(-b^2)=(-1)^(j+1) 2^(j-1)b^j !=0.

This follows directly because the signs of all even-binomial summands
agree after substitution and their binomial sum is 2^(j-1). Consequently
the highest nonzero C order at a fixed b cannot be cancelled by lower
orders. Different positive b have distinct squares and hence coprime
linear factors. No polynomial part is hidden: every real principal term
is proper. Once the zero real channels vanish, the exact reduced
denominator degree is therefore the sum of these local highest orders.

The independent saved-table check reconstructs the signed fold with Python
Fraction, expands every small numerator around z=-b^2, compares the saved
ordinary principal table, and verifies all 359 highest-pole witnesses.
It obtains 311 factors of order 1, 38 of order 2, and 10 of order 3:

    degree_z =311+2*38+3*10=417.

There are 417 nonzero folded entries and 417 ordinary real principal
entries. The raw degree 889 is not the reduced degree. A count of signed
A entries alone would also be insufficient. This saved-table verification
is not a second independent full density aggregation; confidence in the
complete density table additionally uses the proved local identity and
the audited exact aggregation source. The independent full computation
below specifically reconstructs Green and all zero Laurent principal parts.

## 4. Independent full Green calculation by Hermite divided differences

The full audit uses a different algorithm from the pole Taylor expansion.
For f_epsilon(b)=atan(b/epsilon), let D5 denote the ordinary confluent
divided difference on the six sorted nodes. The exact resolvent identity is

    D5[1/(T-i b)]=i^5 product_(j=1)^6(T-i b_j)^-1.

Multiplication by i^-5=-i and taking real parts yields

    Re product_j(T-i b_j)^-1 =D5[b/(T^2+b^2)].

Repeated nodes mean normalized derivatives f^(k)(b)/k!, not a finite
perturbation. The identity extends to repeated and zero nodes by continuity
for each T>0. Integrating from epsilon>0 to infinity gives

    H_epsilon(x0)=(2pi)^-6 sum_B w_B D5[atan(b/epsilon)].

Interchanging this finite divided-difference functional and the time
integral is valid: each evaluation or normalized node derivative of
b/(T^2+b^2) is integrable on [epsilon,infinity). This fixes the sign and
Fourier factor independently of B's complex partial coefficients.

For b!=0, the constant term of atan(b/epsilon) is sign(b)*pi/2, and all
higher normalized node derivatives have only positive epsilon powers.
At b=0 the value is zero but the odd derivatives have exact negative powers:

    f_epsilon^(k)(0)/k!
      =(-1)^((k-1)/2)/(k epsilon^k), k=1,3,5,

and the even derivatives are zero. Therefore each tuple's full principal
and constant expansion is exactly represented by four rational coefficients:
the coefficient of pi/2, and those of epsilon^-1, epsilon^-3, epsilon^-5.
The omitted remainder is O(epsilon). All divided-difference denominators
are fixed nonzero node differences, independent of epsilon.

The independent triangular recurrence divides neighboring rows by their
node differences. When all nodes in a subinterval coincide, it inserts the
correct normalized derivative above. In particular it never treats sign(0)=0
as permission to omit the divergent zero-node derivatives. The coefficients
are multiplied by the complete stored weights and summed exactly over all
40320 tuples once.

The complete result is

    coefficient_(epsilon^-1)=0,
    coefficient_(epsilon^-3)=0,
    coefficient_(epsilon^-5)=0,
    Gamma_DD5 = Gamma_B023 exactly as full rational numbers.

The corresponding B zero-channel relations are A_(0,2), -A_(0,4)/3,
A_(0,6)/5; the independent recurrence agrees with all three. All possible
negative powers through DD order five have been accounted for, not merely
the maximum zero multiplicity observed in the current data. Because the
tuple sum is finite, the O(epsilon) remainders remain O(epsilon). Taking
epsilon down to zero therefore proves

    G(x0)=(2pi)^-6*(pi/2)*Gamma_G=Gamma_G/(128*pi^5).

The complete rational Gamma_G and G*pi^5 are saved in this directory's
`results.json`. No low-precision decimal is a proof input.

## 5. Direct Green integration of the partial table agrees

After the exact zero-real-pole cancellation, every nonzero-pole real term
is absolutely integrable. For j>=2 the complex integral of
(T-i b)^-j is (-i b)^(1-j)/(j-1). Multiplying by i^(j-6) leaves a purely
imaginary value, so its real integral is zero. For j=1 use its real part
directly, avoiding a divergent complex logarithmic integral:

    Re[-i A/(T-i b)]=A b/(T^2+b^2),
    integral_0^infinity = (pi/2) sign(b) A.

Thus Gamma is sum_(b>0) C_(b,1). Independently integrating the ordinary
real principal table uses

    integral_0^infinity (T^2+b^2)^-j dT
      =(pi/2) binom(2j-2,j-1)/(4^(j-1)b^(2j-1)).

The saved-table algebra reconstructs this expression exactly and agrees
with Gamma and with the full Hermite-DD5 result. Integration is in dT,
not dz. The final denominator 128*pi^5 already includes all group weights;
there is no remaining stabilizer or global-negation multiplier.

## 6. Formal tests, actual run, limits, and reproducibility

Before the sole complete tuple pass, exactly the eight B023 registered node
multisets were checked. An independent degree-five Hermite interpolation
matrix reconstructs the four Laurent coefficients and matches the triangular
recurrence. The same interpolation matrix applied to the normalized node
derivatives of b/(T^2+b^2) is compared, as an exact rational function of the
indeterminate T, with Re product(T-i b)^-1. These checks include mixed signs,
single and double zeros, double/triple/quadruple repeats, and six zeros.
No numerical T is substituted in any case.

The actual run exited 0: PID **38086**, terminal session **17091** closed,
start chunk **84d19c**, final chunk **d5bc9e**. The complete symbolic pass
took **0.9492714589869138 seconds**; total script time, including eight
formal cases and saved-table checks, was **1.1692438749887515 seconds**.
The environment was Python 3.15.0, python-flint 0.9.0, SymPy 1.14.0.
The full-pass guard allowed one pass only, at most 900 seconds and 768 MiB;
the actual memory receipt is retained in `run_receipt.json`.

Run in a fresh review output directory with the same frozen input paths:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_symbolic.py

The script intentionally refuses a second complete pass when its run receipt
already exists. This review did not execute B's source, rerun C022's orbit
enumeration, or recheck a previous point's full baseline. It executed one
new complete symbolic audit at C022's x0, explicitly counted as such.

The exact Γ equality strengthens the complete symbolic evidence, but does
not by itself independently reprove every density coefficient. Those are
covered by the rational-identity/source review and saved-table checks as
stated above. The degree-512 feasibility gate is met by 417; it is not a
root certificate, an all-time q bound, or authorization to compute either.
No future stage, numerical minimum, or C6 improvement is included.

All source and archive hashes are pinned. Only this review directory was
written; C022/B023 and earlier F directories remained read-only. There was
no MATLAB, subagent, native engine, positive-time evaluation, Monte Carlo,
global large numerator, Sturm call, git staging, commit, or push.
