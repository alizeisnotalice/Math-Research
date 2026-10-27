"""New exact reconstruction of one fixed layer certificate; no search/import of producer code."""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse
import json

HERE = Path(__file__).resolve().parent
default_input = HERE.parents[1] / "round062" / "compute" / "candidate_certificate.json"
parser = argparse.ArgumentParser()
parser.add_argument("--certificate", type=Path, default=default_input)
args = parser.parse_args()
raw = args.certificate.read_bytes()
input_hash = sha256(raw).hexdigest()
assert input_hash == "b889e35ef909bd979937752c67410bf3913eb76cdf7dc6c720556e093d4c0711"
c = json.loads(raw)
den = 2**200
def down(v):
    return Q((v.numerator * den)//v.denominator, den)
def up(v):
    return -down(-v)
def log_interval(v):
    assert v >= 1
    y = (v-1)/(v+1)
    y2 = y*y
    power = y
    partial = Q(0)
    for j in range(48):
        partial += 2*power/Q(2*j+1)
        power *= y2
    tail = 2*power/(97*(1-y2))
    return down(partial), up(partial+tail)

d, sigma = Q(7,25), Q(3,50)
J, trim = 128, Q(1,2**30)
e, eret, ea = Q(1,2**20), Q(1,4096), Q(1,65536)
P0 = 1 - e - 1280*Q(1,2**48)
Ulo, Uhi = P0*(1-e)/(1+e), (1+e)/(P0*(1-e))
assert d <= Q(1,3) and 4*d >= 1 and P0 > 0
assert Q(c['constant']['d']) == d and c['constant']['J'] == J
assert Q(c['constant']['sigma']) == sigma
assert Q(c['constant']['Ulo']) == Ulo and Q(c['constant']['Uhi']) == Uhi
assert len(c['rows']) == J
total = Q(0)
eta_gaps, tau_gaps, mu_gaps, retention_x = [], [], [], []
last_hi = None
for j, row in enumerate(c['rows']):
    qlo = Q(18,25) + Q(6,25)*Q(j,J)
    qhi = Q(18,25) + Q(6,25)*Q(j+1,J)
    rlo, rhi = qlo+d+trim, qhi+d-trim
    assert row['j'] == j and 0 < qlo < qhi < 1 and rhi > rlo > 1
    if last_hi is not None:
        assert last_hi < rlo
    last_hi = rhi
    for key, value in (('qlo',qlo),('qhi',qhi),('rlo',rlo),('rhi',rhi)):
        assert Q(row[key]) == value
    ll, lu = log_interval(rhi/rlo)
    _, phi = log_interval(1/(1-d/(2*(qlo+d))))
    assert Q(row['log_lower']) == ll and Q(row['log_upper']) == lu
    assert Q(row['phi_upper']) == phi
    mu = up((qlo+d)/qlo*phi**2/(1-e))
    assert Q(row['mu_upper']) == mu
    microceil = (mu.numerator*10**6+mu.denominator-1)//mu.denominator
    chi = Q(microceil,10**6) + Q(1,10000)
    assert Q(row['chi']) == chi and chi > mu
    eta = (1-qhi)/(2*qhi)
    eta_gap = eta-chi
    tau_gap = sigma/2-chi*(1-sigma)
    assert Q(row['eta']) == eta and Q(row['eta_margin']) == eta_gap > 0
    assert Q(row['sigma_margin']) == tau_gap > 0
    rowb = (1-eret)*2*Ulo/qhi-2*e
    M = 2*Uhi/qlo
    assert rowb > 0 and M > 2 and Q(row['rowb']) == rowb
    gamma = down((1-e)*(1-ea)*ll*rowb/M)
    assert gamma > 0 and Q(row['gamma_lower']) == gamma
    total += gamma
    eta_gaps.append(eta_gap); tau_gaps.append(tau_gap); mu_gaps.append(chi-mu)
    V = Q(2,3)*(qlo+d)/qlo*phi**3/(1-e)
    retention_x.append(max(2*phi/(chi-mu),4*V/(eret*(chi-mu)**2)))

assert total == Q(c['gamma_total_lower'])
clean = Q(429,2000)
assert total > clean
raw_c = clean-Q(3,47)
slack = raw_c-Q(3,20)
assert raw_c == Q(14163,94000) and slack == Q(63,94000) > 0
assert Q(c['gamma_minus_3_47_minus_3_20_lower']) == total-Q(3,47)-Q(3,20)
out = {
    "status":"PASS",
    "scope":"new exact reconstruction of one fixed conditional scalar certificate; the accompanying theory supplies original-history/operator validity",
    "input_sha256": input_hash,
    "rows_checked": J,
    "gamma_total_lower": str(total),
    "clean_gamma": str(clean),
    "clean_gamma_margin": str(total-clean),
    "raw_cost_coefficient": str(raw_c),
    "final_target_coefficient": "3/20",
    "full_scalar_slack": str(slack),
    "min_eta_gap": str(min(eta_gaps)),
    "min_tau_gap": str(min(tau_gaps)),
    "min_chi_mu_gap": str(min(mu_gaps)),
    "max_retention_x_sufficient": str(max(retention_x)),
    "no_explicit_dimension_threshold_claimed": True,
}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({"status":"PASS","rows_checked":J,"clean_gamma":str(clean),"final_c":"3/20"}))
