# Round011d: unique contact-guided mass exchange (2026-10-25)

Authorization: parent authorized this one finite direction after locking round011. Only new round011d_* files may change. No random centers, phase grid, extra starts, baseline union rerun, or additional optimization.

Input: round008_certificate_input.json SHA256 8fa317bdb9d2b4209eb18f8fb1bb0fb6bc95abd6861728113c02abec3ab236bf, all 639 keys retained. Target is the entire exact upper endpoint of main audit/022_results.json; snapshot/hash saved before selection.

Set delta=eps2=0.000385935067256 exactly. eps1(s)=eps1+delta*s; eps2(s)=eps2-delta*s. Thus W is exactly constant. P_j(s)=P_j*(1+a_j*s), phi_bj(s)=phi_bj+v_bj*s, b=1,2. One rational direction with 9 geometry velocities and fixed mass exchange.

Select three constraints, one per normal axis: old nonempty receiver lower face B versus upper face A, different boxes and different orbit winners, each raw face inside [P_j/1000,999P_j/1000], tangent overlap >=P_t/10000, nonzero radius-speed forcing. Rank by the rational upper bound for |L_B-U_A|/P_j and then original key IDs. No objective values enter selection. For each independent axis row c=(LC_B-HC_A,P_j*(1[b=1]-1[c=1]),P_j*(1[b=2]-1[c=2])), solve c.u=midpoint(r'_B+r'_A) in minimum Euclidean norm, quantizing velocities to denominator 10^12. Report the exact interval residual; do not call quantized compensation exact.

Before any objective query, bound local weights, P>1, phases in [0,1], all original raw atom winner branches, selected-face clipping and tangent overlap for the entire symmetric path interval. Start s0=1/4, shrink at most eight halvings; stop without union if guards fail. Persist selected constraints, rational direction, guards and step bank/hash before evaluation.

Fixed bank: {0,+s0/2,-s0/2,+s0,-s0}. Baseline 0 reuses audit022 and is not queried. At most four new complete joint unions, one lightweight process, no parallel worker. At each step regenerate all masses/common radii/raw max-min/clipping and use full denominator W. Family contains all old639 plus previously captured round011 patterns, their conservative periodic placements on the whole guarded interval, and existing singleton/pair/triple cluster library. This remains a finite witness family; no exhaustiveness claim about the maximal operator.

Stop at first floating value greater than the complete022 upper endpoint by 1e-10 and immediately rational/interval-certify that point, otherwise after the four points. Certification of the same point is not a new parameter query. No new direction after a failed or exhausted batch. MATLAB only via parent after a concrete direction. Exact certification must compare new lower endpoint to full022 upper, and include explicit index collar and same-mass finite-block/smoothing transfer. Numerical nonimprovement only concerns this direction/bank/family.

Pre-query lock: round011d_query_lock.json. Direction SHA256 f62eefeb0894f3f09f44449de9b6fd6e4aeff84a0cde6d6d31af9b2995c4f36c. Actual selected constraints: axis0 (lower330,upper415), axis1 (lower319,upper315), axis2 (lower310,upper413). s0=1/4; all639 raw atomic extrema and selected clipping/tangent guards PASS. Four nonzero steps exactly +1/8,-1/8,+1/4,-1/4. Family22627 finite keys. No baseline volume query.
