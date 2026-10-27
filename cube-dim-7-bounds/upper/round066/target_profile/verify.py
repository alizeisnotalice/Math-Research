"""One fixed R066 tuple: exact rational gates, no search or sampled integration."""
from fractions import Fraction as F
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
def log_bounds(v, terms=64):
    assert v >= 1
    y = (v-1)/(v+1)
    power, total = y, F(0)
    for j in range(terms):
        total += 2*power/(2*j+1)
        power *= y*y
    return total, total+2*power/((2*terms+1)*(1-y*y))

s, qmin, qmax = F(3,2), F(2,3), F(917,1000)
delta, e, eret, ea = F(1,2**30), F(1,2**20), F(1,4096), F(1,65536)
rlo, rhi = 1+delta, s*qmax-delta
sigma, chi = F(83,1000), F(9,200)
P0 = 1-e-F(1,2**60)
Ulo, Uhi = P0*(1-e)/(1+e), (1+e)/(P0*(1-e))
M = 2*(1+e)*Uhi
b = 2*(1-eret)*(1-e)*Ulo-2*e
mL, gamma = F(637,2000), F(317,1000)
ll, lu = log_bounds(rhi/rlo)
log6lo, log6hi = log_bounds(F(6,5))
mu0lo, mu0hi = F(1,2)-F(5,2)*log6hi, F(1,2)-F(5,2)*log6lo
profile_base = F(221,5000)+F(1,100000)
profile_all = F(2211,50000)
muSimple = F(223,5000)
Q5 = F(34187,187500)
budget = (1-e)*mL*b/M
eta = (1-qmax)/(2*qmax)
raw = gamma-sigma/(1-sigma)
gates = {
    'relative_domain': F(0)<qmin<qmax<1 and 1<rlo<rhi<F(3,2),
    'mass': (1-ea)*ll>mL,
    'P0_positive': 0<P0<1,
    'degrees': M>2 and b>0,
    'profile_mean': mu0hi<F(221,5000),
    'profile_bridge_budget': profile_all-profile_base==F(1,100000),
    'source_mean': profile_all/(1-e)<muSimple,
    'pin_log': log6hi<Q5<F(1,5),
    'source_variance': (2*s/3)*Q5**3/(1-e)<F(1,150),
    'retention_gap': chi-muSimple-F(1,5000)==F(1,5000),
    'deficit': budget>gamma,
    'cap_gap': eta-chi==F(47,183400)>0,
    'tau_gap': sigma/2-chi*(1-sigma)==F(47,200000)>0,
    'raw_cost': raw==F(207689,917000),
    'complete_cost_slack': raw-F(9,40)==F(341,229250)>0,
    'product_shift': 512*F(1,2**80)<F(1,2**60),
}
assert all(gates.values()), gates
out = {
    'status':'PASS', 'scope':'one fixed exact rational scalar/profile certificate; full original-operator proof is separate',
    'gates':gates, 'mass_lower':str(mL), 'gamma':str(gamma),
    'budget':str(budget), 'budget_minus_gamma':str(budget-gamma),
    'profile_mu0_interval':[str(mu0lo),str(mu0hi)],
    'profile_mu_upper':str(F(221,5000)),
    'source_mean_upper':str(muSimple), 'source_variance_upper':'1/150',
    'eta':str(eta), 'cap_gap':str(eta-chi),
    'tau_gap':str(sigma/2-chi*(1-sigma)),
    'raw_cost':str(raw), 'final_c':'9/40', 'complete_cost_slack':str(raw-F(9,40)),
    'explicit_dimension_threshold':None,
}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','gates':len(gates),'gamma':str(gamma),'final_c':'9/40'}))
