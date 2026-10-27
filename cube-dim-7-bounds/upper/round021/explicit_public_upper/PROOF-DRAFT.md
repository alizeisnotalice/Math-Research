# Stein–Strömberg selection argument: candidate explicit finite-dimensional bound

Status: **proof draft frozen for independent audit; no numerical optimization or claimed new accepted bound**. Source: E. M. Stein and J.-O. Strömberg, *Behavior of maximal functions in R^n for large n*, Ark. Mat. 21 (1983), 259–269, DOI 10.1007/BF02384314; relevant complete argument §2, printed pp.260–263 (PDF indices1–4). Original article reproduction: https://archive.ymsc.tsinghua.edu.cn/pacm_download/116/7398-11512_2006_Article_BF02384314.pdf . Author names/title/year are publisher-original article text, not a secondary summary.

The article states only c n log n with unspecified absolute c in Theorem1 (p260). It fixes B*=nB, B**=(n+1)B, B***=(n+2)B. Below is an explicit quantitative derivation from its full selection proof, with a freely chosen dilation a>1 and an exact finite-sum layer-cake identity replacing the article's geometric partition of radius ranges. It is NOT a quotation of a published explicit numerical constant.

## Model and proposed formula

Let n>=1, let D be any open bounded convex body symmetric about0 in R^n, and let |.| denote Lebesgue measure. The Minkowski functional ||.||_D is a norm (D contains0 as an interior point). For x∈R^n,r>0 put B(x,r)=x+rD and

M_D f(x)=sup_{r>0} (|rD|)^(-1)∫_{x+rD}|f(y)|dy.

All radii r>0 are allowed, a common scale multiplies the entire fixed body, and f is any L1 function. The n7 application chooses D=(-1,1)^7. Candidate bound for every a>1:

C_n(a)=[1+((a+2)/a)^n] * [1+((a+1)/a)^n(1+n log a)].

No assumption of radial f, finite support, fixed scale window, differentiability, or discretized input is used.

## Finite selection lemma

Given any finite family {B_alpha=B(y_alpha,r_alpha)}, with all r_alpha>0, select B1,...,BN recursively. If B1,...,B(k-1) have been selected, define the *small-ball increments*

I_j=B_j \ (B1∪...∪B(j-1)).

These I_j are pairwise disjoint measurable sets and ∑|I_j|=|∪selected B_j|. They are NOT increments of the star-ball union. Put w_j=|I_j|/|aB_j|, where aB_j=B(y_j,a r_j) denotes dilation about its center (not multiplication of center coordinates). A remaining ball B_alpha is eligible when

∑_{j<k} w_j 1_{(a+1)B_j}(y_alpha) <=1.

Choose an eligible remaining ball of maximal radius, if any; otherwise stop. Because adding terms only increases the eligibility sum, radii of selected balls are nonincreasing. Zero increments are harmless. The procedure is finite.

### Covering conclusion

If an unselected B_alpha is rejected, then

∑_{j:r_j>=r_alpha} w_j 1_{(a+1)B_j}(y_alpha)>1.

Justification: otherwise it would still be eligible after all selected balls of radius>=r_alpha have been taken, so before any smaller selected ball the maximal-radius rule would select it. If no smaller selected ball is ever taken, termination also forces rejection by the same set of selected balls. Equal radii cause no difficulty because all are included in >=r_alpha.

For x∈B_alpha and a contributing j, convex norm triangle inequality gives

||x-y_j||_D < r_alpha+(a+1)r_j <=(a+2)r_j.

Thus ∑_j w_j 1_{(a+2)B_j}(x)>1 on the union of unselected balls. Integrating that pointwise bound yields

|∪unselected B_alpha| <= ∑_j w_j |(a+2)B_j|
= ((a+2)/a)^n ∑_j|I_j|.

Hence

|∪all B_alpha| <= [1+((a+2)/a)^n] ∑_j|I_j|.    (A)

This is p262, original formula(2.3) and its integration, keeping the exact volume factor instead of its e^2 upper estimate.

### Weighted star overlap conclusion

For every x∈R^n put W(x)=∑_j w_j 1_{aB_j}(x). If W(x)=0 there is nothing to prove. Otherwise choose among contributing selected balls one B_k of **smallest radius**, and translate/dilate so x=0 and r_k=1. (There is a minimum because the family is finite.) All contributing radii r_j>=1, and ||y_k||_D<a.

Split W(0)=I+II by r_j>=a and 1<=r_j<a.

If j contributes to I, then j<k, since r_j>=a>1=r_k, and

||y_k-y_j||_D < a+a r_j <=(a+1)r_j.

Therefore the selection condition at y_k gives I<=1. This is p263 original(2.5), with the range n<=r_j replaced by a<=r_j. The threshold a, rather than1, is essential.

For II, let J={j:0∈aB_j,1<=r_j<a}, set v_j=|I_j|, and define the cumulative raw increment volume

F(t)=∑_{j∈J:r_j<=t} v_j,  1<=t<=a.

Since each I_j⊂B_j and ||y_j||_D<a r_j, for r_j<=t we have B_j⊂B(0,(a+1)t). Disjointness of the SMALL-ball increments gives

F(t)<=|D|(a+1)^n t^n.    (B)

The exact identity r^(-n)=a^(-n)+n∫_r^a t^(-n-1)dt implies, by summing finitely many nonnegative terms,

∑_{j∈J} v_j r_j^(-n)
=a^(-n)F(a)+n∫_1^a F(t)t^(-n-1)dt
<=|D|(a+1)^n[1+n log a].

This needs no differentiability or continuity of F and no extra hypothesis. It is a finite sum identity plus an ordinary integral; it may equivalently be expressed as Stieltjes integration by parts, with the atom at1 included in the boundary term. Thus

II=1/(|D|a^n)∑_{j∈J}v_j r_j^(-n)
<=((a+1)/a)^n[1+n log a].

We obtain

W(x)<=1+((a+1)/a)^n[1+n log a].    (C)

The original p263 uses a bound on each geometric radius bin instead; (B) is its same disjointness observation before partitioning. The finite-sum identity is the quantitative sharpening introduced in this derivation.

## From the finite lemma to all scales and general L1

Take f∈L1, replace it by |f|, and fix lambda>0. M_D f is measurable, indeed lower semicontinuous: for each fixed r>0, convolution with |rD|^(-1)1_{rD} is continuous as a function of x (L1 translation continuity gives a uniform-in-x bound using the bounded kernel), and a supremum of continuous functions is lower semicontinuous. Thus E_lambda={M_D f>lambda} is open.

Let K⊂E_lambda be any compact set. For each x∈K choose an actual radius R_x>0 with ∫_{B(x,R_x)}|f|>lambda|R_xD| (strict supremum level). Put small ball B_x=B(x,R_x/a). It is open, contains x, and its star aB_x=B(x,R_x) has average>lambda. The open small balls cover K, so choose a finite subcover. Apply the finite selection lemma to this finite family.

For every selected B_j,

lambda|I_j| < [|I_j|/|aB_j|]∫_{aB_j}|f|

when |I_j|>0; for zero increments use <=. Hence Tonelli for the finite sum and (C) imply

lambda∑_j|I_j| <= ∫ |f(y)|W(y)dy
<= [1+((a+1)/a)^n(1+n log a)]||f||1.

Combine with (A) and K⊂∪all B_alpha to get

lambda |K| <= C_n(a)||f||1.

Lebesgue inner regularity on the open E_lambda gives the same inequality for E_lambda by taking the supremum over all compact K⊂E_lambda. This also proves its measure finite without presupposing it: the right side is finite and inner regularity applies to an open set even if its measure were initially infinite. There is no truncation in the radii or support; only the finite subcover for each compact level subset was used.

Open vs closed cube averages agree for L1 input because their boundaries have Lebesgue measure0. The conclusion therefore applies to the exact fixed centered cube model of c7.

## Rational, non-decimal checkable specialization awaiting audit

For a=64,n=7 the exact expression is

C_7(64)=[1+(33/32)^7][1+(65/64)^7(1+42 log2)].

The elementary binomial bound (1+x)^7<=1/(1-7x) for0<=x<1/7 gives

1+(33/32)^7 <=1+32/25=57/25,
(65/64)^7<=64/57.

Already log2<1 gives C_7(64)<(57/25)[1+(64/57)43]=2809/25<128, entirely rational after the elementary logarithm inequality. A sharper optional bound log2<7/10 follows from e^(7/10)>1+7/10+(7/10)^2/2+(7/10)^3/6>2; then C_7(64)<10013/125<81. These are algebraic proof candidates, no numerical evaluation/optimization has been run in this subtask. Parent must independently audit before adopting.
