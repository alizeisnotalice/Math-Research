# Main appendix: exactly one incident rank-three span at the new body32 point

Let V={(1,s2,...,s6):si in {-1,1}} and x=(1,1,1,1,1/2,0). Let S2 denote
the union of actual spans of at most two columns of V. This auxiliary lemma
does not certify g, a curvature constant, an integral, or C6.

First x is outside S2. If x=av+bw, its first coordinate gives a+b=1.
Its sixth coordinate cannot be plus or minus a+b, so it forces a=b=1/2.
But then its fifth coordinate belongs to {-1,0,1}, contradicting x5=1/2.
The same argument includes single-column representations by setting b=0.

Consider any independent triple of body columns through x, with unique
coefficients c1,c2,c3. They are all nonzero, because otherwise x is in S2.
The first row is (+++), hence c1+c2+c3=1. If any of rows2–4 differs from
(+++), comparing it with row1 gives a nonempty subset of coefficients
whose sum is zero. It cannot have one element (ci nonzero), or all three
(sum1). Relabel so that c1+c2=0 and c3=1. The sixth row gives
sigma1*c1-sigma2*c1+sigma3=0. Equal sigma1/sigma2 are impossible;
unequal signs force c1=1/2 or -1/2. Every row value is then an integer,
since +/-c1 +/-c2 is in {-1,0,1} and +/-c3 is an integer. This contradicts
x5=1/2. Thus rows1–4 are all (+++), so every column of the triple belongs to

    U={(a,a,a,a,b,c):a,b,c real}.

Its span is U by rank three. Conversely choose the three body columns with
first four entries1 and last pairs (++),(+-),(-+). The minor on rows1,5,6
has determinant -4; coefficients (1/4,1/2,1/4) give x. Therefore U is an
actual incident rank-three frame span, and the preceding argument proves
it is the only one. Its maximal in-plane column set consists of exactly
the four columns with last pairs in {-1,1}².

Use the orthonormal tangential coordinates y=(sum_(i=1)^4 xi/2,x5,x6).
The four in-plane columns become (2,+/-1,+/-1). Thus this identifies the
actual local in-plane frame, rather than merely a three-column subframe.
Quotient columns and any weighted moments must still be treated with their
full original labels as required by076.

There are also uniform rational geometric margins without enumerating all
pairs or triples. For independent r columns A, r=2 or3, the Gram determinant
D=det(A^T A) is a positive integer at most6^r by Hadamard. If x is not in
span A, the integer Gram determinant E=det([A,2x]^T[A,2x]) is positive, so

    dist(x,span A)^2=E/(4D)>=1/(4*6^r).

Consequently dist(x,S2)^2>=1/144, and the squared distance to every
nonincident rank-three span is at least1/864>1/1024. Dependent subsets
reduce to bases; single-column spans are already in some two-column span.
Hence delta>=1/12 and gamma>1/32 in the076 local construction. Taking
rho=1/512 is valid. The column-norm sum is32*sqrt6<80, so the conservative
fixed scalar cutoff R0=1/(512*81)=1/41472 has 2R0*S<2rho as needed.

This completes the incident/nonincident classification for this single new
point. It does not give margins for arbitrary spatial boxes, actual C4
amplitude bounds, full remainder budgets or a body32 semiconcavity constant.
In particular it does not promote the single leading negative norm to a
whole-box or C6 certificate. The new appendix awaits independent review.
