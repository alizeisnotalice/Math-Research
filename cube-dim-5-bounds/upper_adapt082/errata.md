# 082 frozen mechanism wording clarification

mechanism.md says “eight fixed small negative jitters”; the frozen manager.py actually registers seven negative fallback vectors, in addition to midpoint and three positive geometric candidates (11 candidates in total). Candidate order/coordinates and exact280 qualification in the pre-production manager code and job diagnostics are authoritative. No kernel, geometry, registration, chosen point, numerical result or source byte is changed. This corrects documentation only and does not expand the finite candidate list.
