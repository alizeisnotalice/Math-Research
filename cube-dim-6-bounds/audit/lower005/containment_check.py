#!/usr/bin/env python3
"""Rational-interval containment audit for the span-prune shell survivors."""
from fractions import Fraction as Q
from itertools import product
import hashlib
import json
from pathlib import Path

from span_prune import sixth_interval, spacing_intervals

N = 6
HERE = Path(__file__).resolve().parent


def box_intervals(ell, a):
    P = 1
    for e in ell:
        P *= e
    # sixth_interval(P) encloses 2r=P^(1/6), whereas box endpoints use
    # r=P^(1/6)/2.  Keeping this division explicit prevents the earlier
    # factor-two containment error.
    diameter_lo, diameter_hi = sixth_interval(P)
    rlo, rhi = diameter_lo/2, diameter_hi/2
    ans = []
    for e, (alo, ahi) in zip(ell, a):
        if e % 2:
            # min of two interval-enclosed upper endpoints
            half_lo, half_hi = alo/2, ahi/2
            rad_lo, rad_hi = rlo-(e-1)*ahi/2, rhi-(e-1)*alo/2
            hlo, hhi = min(half_lo, rad_lo), min(half_hi, rad_hi)
            hi_tag = "half" if half_hi <= rad_lo else "radius" if rad_hi <= half_lo else "unknown"
            ans.append((Q(0), Q(0), hlo, hhi, "zero", hi_tag))
        else:
            llo = max(Q(0), e*alo/2-rhi)
            lhi = max(Q(0), e*ahi/2-rlo)
            rad_lo, rad_hi = e*alo/2-rhi, e*ahi/2-rlo
            lo_tag = "zero" if rad_hi <= 0 else "radius" if rad_lo >= 0 else "unknown"
            ans.append((llo, lhi, alo/2, ahi/2, lo_tag, "half"))
    return ans


def certify_contained(new, old):
    # old.lo <= new.lo and new.hi <= old.hi, sufficient interval directions.
    return all((olo_hi <= nlo_lo or (olo_tag == nlo_tag and olo_tag in ("zero", "half"))) and
               (nhi_hi <= ohi_lo or (nhi_tag == ohi_tag and nhi_tag in ("zero", "half")))
               for (nlo_lo, nlo_hi, nhi_lo, nhi_hi, nlo_tag, nhi_tag),
                   (olo_lo, olo_hi, ohi_lo, ohi_hi, olo_tag, ohi_tag) in zip(new, old))


def certify_not_contained(new, old):
    # A single rigorously separated reverse inequality rejects containment.
    return any(olo_lo > nlo_hi or nhi_lo > ohi_hi
               for (nlo_lo, nlo_hi, nhi_lo, nhi_hi, _, _),
                   (olo_lo, olo_hi, ohi_lo, ohi_hi, _, _) in zip(new, old))


def main():
    a = spacing_intervals()
    prior = json.loads((HERE/'span_prune_results.json').read_text())
    old = [(e, box_intervals(e, a)) for e in product(range(1,5), repeat=N)]
    output = {"status": "strict_interval_containment_audit_corrected_radius", "shells": {}}
    for key in ('m5', 'm6'):
        rows = []
        for ell_list in prior[key]['new_shell_survivors']:
            ell = tuple(ell_list)
            new = box_intervals(ell, a)
            new_nonempty = all(nhi_lo > nlo_hi
                               for nlo_lo, nlo_hi, nhi_lo, nhi_hi, _, _ in new)
            assert new_nonempty
            contained = [list(e) for e,b in old if certify_contained(new,b)]
            all_rejected = all(certify_not_contained(new,b) for _,b in old)
            status = ("contained_in_old" if contained else
                      "not_contained_in_any_old" if all_rejected else
                      "interval_unresolved")
            rows.append({"ell": ell_list,
                         "contained_in_old_lte4": bool(contained),
                         "old_container_samples": contained[:3],
                         "new_canonical_half_box_certified_nonempty": new_nonempty,
                         "certified_not_contained_in_any_old_lte4": all_rejected,
                         "status": status})
        output['shells'][key] = rows
    output['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'containment_check_corrected_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__ == '__main__':
    main()
