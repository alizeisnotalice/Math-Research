#!/usr/bin/env python3
"""Three-band counterpart of band_search.py; finite model only."""
from __future__ import annotations
import json
from pathlib import Path
from scipy.optimize import differential_evolution
from band_model import BandGeometry

def geom(z, gamma):
    center, gap, outer, rc, ro = z
    return BandGeometry(((-center-gap-outer,-center-gap),(-center,center),(center+gap,center+gap+outer)),
                        (ro,rc,ro),gamma)

def run(gamma=.5):
    def score(z):
        try:
            answer=geom(z,gamma).optimize_heights(margin=1e-7)
            return -answer['J'] if answer['success'] else 1e3
        except (ValueError,FloatingPointError): return 1e3
    ans=differential_evolution(score,[(.02,3.),(0.,3.),(.02,8.),(.03,3.),(.03,3.)],
                               seed=320321,popsize=12,maxiter=80,polish=True,workers=1)
    g=geom(ans.x,gamma); inner=g.optimize_heights(margin=1e-7)
    return {'gamma':gamma,'best_vector':ans.x.tolist(),'geometry':{'bands':g.bands,'radii':g.radii},
            'inner':inner,'de_nfev':int(ans.nfev),'de_success':bool(ans.success),'de_message':str(ans.message)}

if __name__=='__main__':
    result=run();Path('band_search_three_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('J={:.12g} nfev={} success={}'.format(result['inner']['J'],result['de_nfev'],result['de_success']))
