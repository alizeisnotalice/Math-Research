# Time bands to a continuous g lower certificate

Let c be on the positive face x1=1, epsilon=a_l, and
a_k=2^k/1000 (0<=k<=7), T=a_7. Write

    Delta_k = H_(a_k)(c)-H_(a_(k+1))(c) >= 0.

For any x with ||x-c||_infinity<=r, the previously audited real comparison
gives p_t(x)<=exp(r/t)p_t(c). Hence

    J2_epsilon(x) <= U2_l(r)
      := sum_(k=l)^6 a_k^-2 exp(r/a_k) upper(Delta_k)
          + T^-2 exp(r/T) upper(H_T(c)).

All terms are nonnegative and the tail is the full [T,infinity) tail. This
formula holds throughout the closed clipped box, not only at its center.
The center lower moment bound is sum a_(k+1)^-2 lower(Delta_k); its upper
bound is U2_l(0). Negative lower endpoints from Arb subtraction are replaced
by zero, using the exact positivity of the time-band integral.

The round006 upper curvature estimate therefore implies

    partial_ii H_epsilon(x) <= (4225/4224) U2_l(r).

For the g lower bound, a stronger opposite-sign bound is useful. In the
coarea representation, R=A^T/32 and ||Rh||_2^2=||h||_2^2/32. Since
(log f_t)'' >= -2/t^2, the product F obeys

    D_h^2 F/F = (D_h log F)^2+D_h^2 log F
             >= -||h||_2^2/(16t^2).

Positive coarea/time integration gives the full quadratic-form inequality

    D_h^2 H_epsilon(x) >= -||h||_2^2 J2_epsilon(x)/16.

The local derivative domination used in round006 justifies differentiating
under both integrals. Every segment c+s h in the clipped convex box remains
inside it. Integrating the second derivative with the positive Taylor kernel
1-s therefore gives

    H_epsilon(c+h) >= H_epsilon(c)+grad H_epsilon(c).h
                       - ||h||_2^2 U2_l(r)/32.

The implementation bounds the linear expression using the actual five
intervals [max(-r,-c_i),min(r,1-c_i)], with the first face coordinate fixed.
It bounds ||h||_2^2 by the sum of their squared maximum absolute endpoints.
Every acceptance is an Arb comparison of this lower enclosure with 1/60000.
There is no substitution J2<=epsilon^-2 H in this certificate.

For comparison, the same center's best original multiplicative radius among
the eight tau values is obtained from exp(-r/tau) H_tau(c)>1/60000. This
comparison uses the already computed time nodes, not a rerun of old boxes.
The time-band multiplicative lower bound is also computed as a separate
check, and is never confused with the Taylor certificate.

The gradient is evaluated by the exact real node derivative identity

    partial_(b_j) DD5 f(b_1,...,b_6)
       = DD6 f(b_1,...,b_6,b_j).

At distinct nodes ordinary real automatic differentiation of the triangular
divided-difference recurrence gives the same derivative. Repeated nodes use
the displayed identity and the true confluent atan derivatives. The signed
permutation contribution is differentiated with the exact node slopes
sign_j*r_j/q_j. Global sign-paired group terms have identical H and gradient
contributions after this chain rule, so processing half the group with
multiplicity two is an exact reduction.

These certificates cover only the listed local box and its symmetry images.
Their complement on the positive face remains unproved. Homogeneity extends
covered directions into the cube; the general L1/all-scale C6 interface and
the I certificate remain separate primary integration obligations.
