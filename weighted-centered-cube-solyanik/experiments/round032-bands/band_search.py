#!/usr/bin/env python3
"""Search two/three disjoint band geometries; writes signed checkpoint JSON."""
from __future__ import annotations
import hashlib, json, platform, sys
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
from band_model import BandGeometry

SOURCE = hashlib.sha256(Path(__file__).read_bytes() + Path('band_model.py').read_bytes()).hexdigest()

def geometry_from_vector(z, gamma):
    # Two symmetric disjoint bands [gap/2,gap/2+length], reflected; common radius.
    gap, length, radius = z
    bands = ((-gap/2-length, -gap/2), (gap/2, gap/2+length))
    return BandGeometry(bands=bands, radii=(radius, radius), gamma=gamma)

def run(gamma=0.5):
    def score(z):
        try:
            g=geometry_from_vector(z,gamma)
            ans=g.optimize_heights(margin=1e-7)
            return -ans['J'] if ans['success'] else 1e3
        except (ValueError, FloatingPointError):
            return 1e3
    result=differential_evolution(score, [(0.0,4.0),(0.05,8.0),(0.05,4.0)],
                                  seed=320320, popsize=12, maxiter=90, polish=True, workers=1)
    g=geometry_from_vector(result.x,gamma); inner=g.optimize_heights(margin=1e-7)
    return {'source_sha256':SOURCE,'gamma':gamma,'best_vector':result.x.tolist(),
            'geometry':{'bands':g.bands,'radii':g.radii},'inner':inner,
            'de_success':bool(result.success),'de_message':str(result.message),'de_nfev':int(result.nfev),
            'python':sys.version,'platform':platform.platform()}

if __name__=='__main__':
    data=run(); Path('band_search_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print('J={:.12g} nfev={} success={}'.format(data['inner']['J'],data['de_nfev'],data['de_success']))
