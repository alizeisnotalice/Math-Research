# Mixed three-critical local arithmetic

Run wolframscript -file cap_convolution.wl with WolframScript 15 ARM. The retained log is the actual successful local execution, not a numerical fit.

For ell(z)=Min[w,Max[0,z-w/2+b]] integrated over [-w/2,w/2], the three b-branches [0,w], [w,2w], and [2w,infinity) give b^2/2, w^2-(2w-b)^2/2, and w^2. The script checks exact branch integrals, endpoint continuity, monotonicity, and the supplied lower bound when tau w<=b<=w. The later positive-part difference and normalized resolved-target clipping inequalities have empty violation sets.

Only the local clipped target bound and positive-part inequality are needed in the new proof. The full convolution supports mechanism selection; it is not inserted as an unproved full-width price. These scalar checks do not prove common-source capture, set disjointness, all-input limits, or any missing tail estimate. Those steps are covered by the proof and independent review, whose conclusion is restricted to each fixed two-narrow-axis window.
