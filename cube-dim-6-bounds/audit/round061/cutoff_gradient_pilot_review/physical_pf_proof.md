# Independent physical time-pole reference for F061

The accepted fan is fixed by hash. In one cone and one signed permutation, put
alpha=|det R|/product q_j, b_j=x.r_j/q_j and d_jk=partial b_j/partial x_k.
Its band gradient contribution, before alpha/(2 pi)^6, is exactly

    Re integral_a^T i product_j(t-i b_j)^-1
                         sum_j d_jk/(t-i b_j) dt.              (1)

This is direct physical differentiation of the six-factor density product.
There are finitely many fan/group terms; a>0 uniformly bounds the denominators
away from zero, so time differentiation is justified. No divided difference,
copied-node table, finite difference, or Euler identity is used as a gradient
oracle.

Group equal rational b_j into distinct real beta_g with multiplicities m_g and
delta_gk=sum_(j in g) d_jk. This grouping retains labels through the delta sum;
unequal physical derivatives at an equal pole are not dropped or equated.
Let

    H_k(z)=product_g(z-beta_g)^(-m_g)
                     sum_g delta_gk/(z-beta_g).

For each g let Q_g(s)=product_(h!=g)(beta_g-beta_h+s)^(-m_h). The principal
part coefficients in H_k=sum_g sum_(r=1..m_g+1) A_gkr/(z-beta_g)^r are

    A_gkr = delta_gk Q_g[m_g+1-r]
       + sum_(h!=g) delta_hk
           sum_(j=0..m_g-r) Q_g[j] (-1)^(m_g-r-j)
                                 /(beta_g-beta_h)^(m_g-r-j+1), (2)

where the second term is zero for r=m_g+1. A_gkr are real rational numbers for
rational inputs. The Q series is obtained by multiplying

    (Delta+s)^(-m)=sum_l>=0 (-1)^l binomial(m+l-1,l)
                                      Delta^(-m-l) s^l.

Substituting z=t/i shows that the integrand in (1) has coefficients
i^(r-6) A_gkr at pole t=i beta_g. In particular the simple-pole coefficient is
purely imaginary (-i A_gk1), so no real logarithm is missing from its real
integral. The primitive contribution for each r is

    Phi_1(beta)=atan(beta/a)-atan(beta/T),
    Phi_r(beta)=Re{ i^(r-6) [(a-i beta)^(1-r)
                                -(T-i beta)^(1-r)]/(r-1) }, r>=2. (3)

Equation (3) holds at beta=0 without a limit or division by beta. The reference
implements its higher primitives by complex reciprocal powers, with real and
imaginary Arb arithmetic; it does not use the C016 rational jet recurrence.
The simple primitive uses two separate principal arctangents, not C016's atan
quotient. Both endpoints are positive, so their phases lie in (-pi/2,pi/2).

The complete cone gradient is sum_(g,r) A_gkr Phi_r(beta_g). The implementation
shares pole-series coefficients across six physical fields but evaluates each
field. It explicitly enumerates all 64 sign masks, all 720 permutations and
all 21 cones with alpha (not doubled alpha), providing independent validation
of the native sign-pair reduction and physical coordinate mapping.

Five new scalar fixtures exercise distinct poles, repeated zero/two clusters,
all zero, six equal nonzero poles, and multiplicities 3/2/1. A signed Hilbert
matrix supplies six linearly independent physical directions. The exact oracle
uses a separate SymPy apart decomposition of H_k and symbolic primitive checks;
the C++ program uses local series formula (2). This tests the code independently
of its coefficient construction, while the displayed analytic argument proves
the formula for all positive-window inputs.

At preparation only, each genuine case computes 128 prespecified tuples, not a
complete field. These normalized partial sums and exact tuple indices will be
reused in any later complete aggregation. Later work must skip those tuples,
so a full gradient covers 967680 tuples exactly once; it cannot merely overlap
intervals after silently rerunning the benchmark work.

The two old C016 full-gradient outputs remain unverified component by component
until the independent 6+6 complete sums have been compared with their exact
native dyadic intervals. Euler supplies only one scalar relation among six
components and is insufficient. Neither a partial benchmark nor the completed
cross-check, if later released, establishes global coverage or a new C6 bound.
