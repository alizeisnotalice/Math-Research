# CI4 continuous integral route registration

2026-10-25. Ownership is only this directory; independent clone, no shared ledger edits.

Model fixed: n=4, eight unnormalised body directions, Lebesgue, common r, all scales, general L1. Target I4 = integral_[0,1]^3 G(1,y,z,w)^(4/3) <= 19/100000. No tau or g search. Source Green and representation hashes are fixed by 011L registration.

Locked initial diagnostic: exact sign-hyperplane polytope partition of the ordered face tetrahedron 0<=w<=z<=y<=1. The 30 interior node-zero hyperplanes are derived from the locked 240 rays. Incrementally split convex polytopes with exact rational vertices/edges; stop if 4096 regions would be exceeded. Regions are only geometry diagnostics until continuous evaluation succeeds. Whole-simplex node intervals use their true affine vertex hull, not an axis box. Share upper_theory's sign-stable DD evaluation; do not duplicate that algebra. Generic point values or quadrature differences are not bounds.

Cap mechanism: audited round011_face_caps.md / audit030, nearest-corner Voronoi restriction, epsilon=1/10000000000. On [0,1]^3 the only body-line intersection is (1,1,1). A radius-epsilon ball cuts one octant, so its cap is at most E4(epsilon)/8. This fraction requires a new direct radial geometry check and does not inherit 030's full-ball count automatically. Whole signed face surface is 64 copies of the unit positive face octant integral; corner count matches 8 faces x8 corners. Ordered tetrahedron is one sixth of the positive octant, not one eighth or one forty-eighth of it.

After the initial finite diagnostic, at most one representation revision is allowed, and must respond to an explicit failure reason (sign-change/cross-denominator cancellation/interval excess). No mechanical grid expansion. If the integral is not enclosed, report every unpaid domain/error item and keep U4=16. Any exact local point evidence remains point evidence. Any evaluator use remains conditional on the full U11 representation audit.
