"""Selected U02 candidates, two finite multiprecision hybrid expressions.

This validates candidate values. mp.quad is not rigorous interval quadrature.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp

CANDIDATES={3:'0.47716536955561023',4:'0.47095317569629574'}
ROWS=[]
for digits,split_text in ((40,'0.12'),(60,'0.20')):
    mp.mp.dps=digits
    split=mp.mpf(split_text)
    image_m=8
    spectral_k=40
    odds=[mp.mpf(2*k+1) for k in range(spectral_k)]
    eigen=[(mp.pi*k/2)**2 for k in odds]
    def survival(t):
        if not t:return mp.mpf(1)
        if t<split:
            return 1-2*mp.fsum([(-1)**k*mp.erfc((2*k+1)/(2*mp.sqrt(t))) for k in range(image_m+1)])
        return 4/mp.pi*mp.fsum([(-1)**k*mp.exp(-eigen[k]*t)/odds[k] for k in range(spectral_k)])
    for n,a_text in CANDIDATES.items():
        a=mp.mpf(a_text)
        coeff=[mp.cos(k*mp.pi*a/2) for k in odds]
        distances=[((a-4*m)**2,(a-(4*m+2))**2) for m in range(-image_m,image_m+1)]
        def heat(t):
            if not t:return mp.mpf(0)
            if t<split:
                return mp.fsum([mp.exp(-u/(4*t))-mp.exp(-v/(4*t)) for u,v in distances])/mp.sqrt(4*mp.pi*t)
            return mp.fsum([coeff[k]*mp.exp(-eigen[k]*t) for k in range(spectral_k)])
        mass=mp.quad(lambda t:survival(t)**n,[0,split/3,split,mp.mpf('.5'),1,mp.inf])
        green=mp.quad(lambda t:heat(t)**n,[0,split/3,split,mp.mpf('.5'),1,mp.inf])
        cost=mass/(2**n*a**n*green)
        ROWS.append({'n':n,'dps':digits,'split':split_text,'a':a_text,
                     'unit_mass':mp.nstr(mass,digits-5),'unit_corner_green':mp.nstr(green,digits-5),
                     'cost':mp.nstr(cost,digits-5)})
print(json.dumps({'status':'multiprecision numeric evidence, not interval',
                  'rows':ROWS,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
