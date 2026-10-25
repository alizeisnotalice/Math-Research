# Round 003: certified second-spacing refinement

The retained construction is the round-001 one-atom rectangular lattice
with only its second spacing changed:

```
a_2 -> (10001/10000) a_2,
a_2=(2+u^5 v)/4,  u^6=2, v^6=3, u,v>0.
```

The author-side exact number-field reconstruction regenerated all 127
nonempty candidates, retained 96 under exact containment and inclusion,
and computed the complete six-dimensional union.  Its interval is saved in
`second_spacing_certificate_results.json`; in particular its lower endpoint
is greater than `2.0651669810405669`.

The independent audit then replaced the author arithmetic with strict
inward rational boxes and an integer union.  It passed with the simple
strict rational conclusion

```
C6 > 25814587263/12500000000 = 2.06516698104.
```

The independent result is the read-only main-project artifact
`audit/lower003/rational_inward_candidate_results.json`, whose producer is
`rational_inward_candidate.py`.  It explicitly compares its lower bound to
the full exact R6-001 upper endpoint, not its decimal display, and supplies
an actual finite smooth-L1 witness with strict threshold and positive bump
width.

The transfer argument in `second_spacing_transfer.md` gives the general
Lebesgue L1/all-scale passage for the same common-radius witness boxes.
Together, the author certificate, independent inward certificate, and
transfer proof establish this round's lower improvement.  It remains an
evaluation/refinement of the rectangular-lattice mechanism, not a claim
that the construction family or public literature has been optimized.
