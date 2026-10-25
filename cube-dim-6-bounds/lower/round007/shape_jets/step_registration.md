# One selected finite step after the positive direction

The registered direction scan stopped on its third direction,
`h=(0,1,0,0,0,0), D=E=0`, with a strictly positive exact derivative.
No remaining directions are evaluated.

Select exactly one step `t=1/1000000`, giving
`a2_new=(1000001/1000000)*(u^3+u^4*v^3)/4`; all other spacings and all
unit masses stay as before. This is a candidate pending E's independent
rational inward and finite-L1 certification.

The D follow-up only checks endpoint order at this one step. Include
all unclipped lower/upper endpoints and moving clipping faces from the
full 382 baseline nonnegative-width raw boxes. Only original coordinate
2 moves. If every adjacent distinct lex-jet endpoint retains strict order
at the selected step, all orders remain valid for 0<t<=1/1000000 because
these endpoints are affine. Hence the witness union volume is exactly
affine on that interval and its ratio is R0+t*R1. This check does not run
another union algorithm or scan step sizes. Failure stops and is reported;
it does not silently choose a smaller step.

The arithmetic source, direction result, original certificate, registration,
and this check's source are hashed before execution. The first-order
certificate's baseline raw family is reconstructed, without extra union
calls. E owns the separate independent finite-volume and L1 acceptance.
