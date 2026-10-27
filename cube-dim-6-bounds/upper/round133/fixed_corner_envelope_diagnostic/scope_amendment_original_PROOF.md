# C133: fixed robust-family corner obstruction, no new exponentials

For a fixed time band W and a fixed C132 selected source center c, let T_(W,c)(x) be the exact robust tangent value from the original saved value/gradient interval rectangle and the A057 f-clamp formula. This is a certificate value and lower surrogate for F_W(x), not the true F_W(x). Every C132 original box corner v has a saved Arb enclosure [l_(W,c,v),u_(W,c,v)] for that value. Set

    L_W(v)=max_c l_(W,c,v), U_W(v)=max_c u_(W,c,v),
    L(v)=sum_W L_W(v), U(v)=sum_W U_W(v).

Then L(v)<=sum_W max_c T_(W,c)(v)<=U(v). The max lower and max upper may have different maximizing centers; this remains a valid enclosure. Missing bands are identically0 in this fixed family. Distinct centers of one band are never added.

For any later axis box B containing this corner v, each same-source robust complete-box certificate satisfies l_(W,c,B)<=inf_(x in B) T_(W,c)(x)<=T_(W,c)(v). Therefore sum_W max_c l_(W,c,B)<=sum_W max_c T_(W,c)(v)<=U(v). Thus if U(v)<target, this particular candidate-jet family and robust formula cannot certify target on any such child box. An exhaustive closed axis subdivision containing the original corner must include it in at least one child, so merely subdividing cannot cure this fixed-family obstruction there. Changing jets, adding centers/bands, exploiting correlations beyond the same independent interval rectangle, or using a different proof is outside the excluded family.

This is never an upper bound on F_W or G, a genuine Green-floor counterexample, or a C6 conclusion. The family surrogate may be strictly below the true function. Conversely, if all32 saved corner lower bounds are at least target but the original wholebox certificate failed, this only marks a finite local-subdivision candidate: max-of-center robust tangents has no asserted concavity, so interior coverage is not established by paid corners. Intervals that do not separate the target remain unresolved.

The unique registered batch is the fixed93 original boxes times32 corners, with targets1/63500 and1/65000. All calculations are exact Fraction/max/sum/comparisons on the existing C132 outward endpoints. No exponential, native kernel, checker, or new subdivision is run. Input files and program hash are registered before the unique start marker. Any structural/hash failure stops without a retry.
