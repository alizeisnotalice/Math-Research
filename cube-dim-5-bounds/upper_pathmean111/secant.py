"""Convex scalar-R secant improvement using exactly the existing catalog.

F_{k,d}(x)=exp(k asinh(x/d)) has F''>0 for k>=1,x>=0.
R=F_{2,4}F_{8,16} is therefore increasing and convex. No transcendental call.
"""
from fractions import Fraction as F
import pathmean as P


def secant_envelope(plan,catalog,old_mean_factor=None):
    base=P.evaluate_plan(plan,catalog,old_mean_factor)
    if plan['profile']['h']==0:
        return {**base,'secant_J_upper':F(0),'positive_endpoint_weights':[]}
    upper=F(0);records=[]
    for cell in plan['cells']:
        l,r,w=cell['left'],cell['right'],cell['weight']
        pieces=[p for p in plan['profile']['pieces'] if p['lo']<=l<r<=p['hi']]
        if len(pieces)!=1:raise ValueError('Cell crosses CDF breakpoint')
        wrho=P.integral([F(0),*pieces[0]['W']],l,r)
        alpha=(r*w-wrho)/(r-l);beta=(wrho-l*w)/(r-l)
        assert alpha>=0 and beta>=0 and alpha+beta==w
        Rl=F(1) if l==0 else P.Q(catalog[str(l/plan['a'])]['upper'])
        Rr=P.Q(catalog[str(r/plan['a'])]['upper'])
        Rl=min(Rl,Rr)  # monotonicity permits this safe rounding tightening.
        upper+=alpha*Rl+beta*Rr
        records.append({'left':l,'right':r,'alpha':alpha,'beta':beta})
    upper=min(upper,base['J_upper']);assert base['J_lower']<=upper
    profile=plan['profile'];a=plan['a'];R=P.Q(catalog[str(profile['h']/a)]['upper'])
    factor=min(base['factor_upper'],R,1+profile['mu_infinity']/a+F(25,24)*upper/a**2)
    return {**base,'J_upper':upper,'factor_upper':factor,'secant_J_upper':upper,
            'curvature_ratio_upper':upper/base['old_J_upper'],
            'darboux_gap':upper-base['J_lower'],'positive_endpoint_weights':records}
