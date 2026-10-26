# A059 independent mathematical review

The fixed model is the original 32 normalized body directions v=(1,s2,...,s6), all with unit Cauchy scale t. Let Omega={(1,y):y in [0,1]^5}, using ordinary five-coordinate face Lebesgue measure, rho=1/64, q=6/5 and epsilon=1/1000. No symmetry multiplier is added in this review. Every source budget uses this same positive face and measure.

## 1. Why exactly 31 unsigned distances suffice on the positive face

For two distinct normalized body vectors v,w, put A={i:v_i=w_i} and B={i:v_i=-w_i}. The set A contains coordinate 1 and B is nonempty. Their plane is the orthogonal direct sum of the signed line supported on A and the signed line supported on B. Consequently, for any x>=0,

    ||Proj_span(v,w) x||^2
      = (sum_A v_i x_i)^2/|A| + (sum_B v_i x_i)^2/|B|
      <= (sum_A x_i)^2/|A| + (sum_B x_i)^2/|B|
      = ||Proj_UA x||^2,

where U_A=span(1_A,1_B). Thus dist(x,U_A)<=dist(x,span(v,w)). This inequality remains valid with zero coordinates, equal coordinates, or a vanishing signed sum.

Conversely U_A is itself the body-pair plane spanned by ones and 1_A-1_B. Taking minima proves both inequalities needed for

    dist(x,S2)=min_(A contains 1, A proper) dist(x,U_A), x>=0.

One-body spans and the origin are already contained in pair planes, so defining S2 with “at most two” directions does not alter this minimum. There are 2^5-1=31 such sets A, with counts binom(5,k-1)=(1,5,10,10,5) for k=1,...,5. The convention 1 in A selects one representative of each unordered partition.

This does not replace the full-space S2 by 31 planes. To see the actual full-space count, any normalized body vector r in span(v,w) can be written r=av+bw. The first coordinate gives a+b=1. A coordinate in B gives a-b=+1 or -1; hence r=v or r=w. Therefore a body-pair plane contains exactly its two normalized body directions. Distinct unordered pairs give distinct planes, totaling binom(32,2)=496. For a signed x, the projection-dominance inequality above is not asserted.

The same one-block argument gives |v.x|<=sum_i x_i on x>=0. Hence the all-ones line is nearest among all body lines, so the cap in the current theorem agrees exactly with the positive-face S1 cap paid in 042.

## 2. Exact priority partition, all boundaries included

Write d_A=dist(x,U_A), b_A=mean_B x and T_A={d_A<rho}. On Omega, 0<=b_A<=1, including singleton B. The three ranges

    [0,1/4], (1/4,3/4], (3/4,1]

are disjoint and exhaustive. Thus the raw unions E,D,H satisfy E union D union H = union_A T_A, even though the raw unions can overlap each other through different A witnesses.

Every U_A contains the ones line; therefore d_A<=dist(x,onesline), and the cap C is contained in every T_A. Now assign a point to the first true entry of C,E,D,H, and to F if none is true. This is precisely

    Omega_C=C,
    Omega_E=E minus C,
    Omega_D=D minus (C union E),
    Omega_H=H minus (C union E union D),
    Omega_F=Omega minus (C union E union D union H).

The five sets are measurable, disjoint and have union Omega point by point. There is no appeal to ignoring a null seam. In particular b_A=1/4 goes to that A's low range and b_A=3/4 to its middle range, subject only to an earlier global priority witness. A point with d_A=rho is outside that specific open tube; it is still assigned using another tube if one is closer, or otherwise goes to F. The cap boundary is handled by the same rule. Zero coordinates and ties remain ordinary points in this partition.

Using C subset union T_A and the distance equality from section 1,

    Omega_F=Omega minus union_A T_A
           ={x in Omega: min_A d_A>=rho}
           ={x in Omega: dist(x,S2)>=rho}.

For completeness the stated high-region coordinate bound is valid on each raw high tube. Decompose the A block as x_i=alpha+e_i, where sum_A e_i=0 and x1=1. Then x_i=1+e_i-e1>=1-sqrt(2)||e||>1-sqrt(2)rho. On B, x_j=b_A+z_j>b_A-rho>3/4-rho=47/64. Since 1-sqrt(2)/64>47/64, every coordinate is at least 47/64. This geometric containment provides no Green budget by itself.

## 3. Correct source-domain and time inclusions

The saved 042 theorem pays the complete t>0 Green function on C (equivalently its S1 cap) with the stated b_C. The saved 056 theorem and its final independent audit pay the same complete G on the union of E_A with 0<=mean_B<=1/4 and d_A<rho; this is the raw E of 059. Therefore restricting the 056 payment to Omega_E=E minus C is safe. No subtraction of an estimated overlap contribution is required.

The saved 051 theorem pays only S_epsilon=integral_0^epsilon p_t dt on the union of d_A<rho and 1/4<=b_A<=3/4. The new raw D has the smaller range 1/4<b_A<=3/4, and Omega_D is smaller still. Hence the old short-time norm upper b_D applies, but it supplies no norm bound for the complete G on Omega_D.

Write L_epsilon=integral_epsilon^infinity p_t dt. Positivity gives G=S_epsilon+L_epsilon on every point as extended nonnegative quantities. If a future proof gives ||L_epsilon||_Lq(Omega_D)<=L_D, then Minkowski on this same domain gives

    integral_(Omega_D) G^q <= (b_D+L_D)^q.

The temporal pieces share a spatial domain; their q-power integrals cannot simply be added. Conversely the spatial priority blocks are disjoint and use the same complete G, so their q-power integrals do add:

    I=integral_Omega G^q
      <= b_C^q+b_E^q+(b_D+L_D)^q+I_H+I_F,

once I_H and I_F denote valid upper bounds on the full-G integrals of their respective blocks. Without those upper bounds the decomposition is still an exact extended-integral identity, but no finite global bound follows. The same distinction applies to the remaining all-scale source-transfer assumptions, which this local ledger does not reprove or replace.

## 4. Meaning of the saved integer-root arithmetic

Let P be the saved rational upper enclosure for b_C^q+b_E^q. The saved paid-B enclosure is for P^(5/6), so it is a conservative norm budget for the already paid union Omega_C union Omega_E. It is not B for the whole face and not a lower estimate of the actual Green norm.

Condition on a future same-model proof g>=1/65000. For B=I^(5/6), the target 64 has critical Bcrit=4/48750, because

    (6 Bcrit)^6 = 64^5 (1/65000)^6.

Let J be the saved lower enclosure of Bcrit^(6/5), and R=J-P>0. If a future aggregate unpaid budget

    (b_D+L_D)^q+I_H+I_F < R,

then I<J<=Bcrit^q and the cost inequality is strict. The saved remaining-B allowance is a lower enclosure of R^(5/6), hence offers a conservative sufficient threshold for the same aggregate remainder. It is neither the ordinary difference Bcrit-paid_B nor a measurement/bound already established for the unpaid regions. If a non-strict future test is desired, its validity additionally uses strictness of the recorded lower root inequality J^5<Bcrit^6, which the saved-certificate checker explicitly inspects.

The root certificates are verified by rational integer powers only. No root search, old scalar budget evaluation, original 059 checker execution, or kernel call is needed. Their exact endpoints, not displayed decimals, carry the arithmetic proof.

The unconditional conclusion of the review is limited to the exact spatial account and valid partial-budget composition. Full I, global g and a new C6 bound remain unproved.
