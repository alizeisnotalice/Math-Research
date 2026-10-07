from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import hashlib
import json

BASE = Path(__file__).resolve().parent
PREFIX = 'packet_coverage_feasibility_20261007'
records = []
total = 0
previous = {}

def digest_fraction(v):
    raw = v.numerator.to_bytes((v.numerator.bit_length()+7)//8, 'big') + b'/' + v.denominator.to_bytes((v.denominator.bit_length()+7)//8, 'big')
    return {'sha256': hashlib.sha256(raw).hexdigest(), 'numerator_bits': v.numerator.bit_length(), 'denominator_bits': v.denominator.bit_length()}

for stage, (n, terms, radii) in enumerate([(64,24,(6,10)),(256,32,(7,12)),(1024,40,(8,14))],1):
    checks = []
    def check(name, predicate):
        global total
        assert predicate, (stage, name)
        checks.append(name)
        total += 1
    root = isqrt(n)
    theta, eta = F(1,root), F(1,4)
    eps, c, tau = F(1,n*n), 1-F(1,4*n*n), 1-F(1,n)
    D, a, b, h = 2*n**3, 1, 2, F(1,2)
    check('square dimension', root*root == n)
    check('positive curvature density and threshold', 0 < tau < c < 1)
    check('full source is high-density for stripe input', n*c > tau*root)
    check('source quadratic loss bound', eps*F(1,4) == 1-c)
    check('both selected candidate endpoints share stripe phase', F(a,2)/h == 1 and F(b,2)/h == 2)
    check('source box endpoints share stripe phase', F(D,2)/h == 2*n**3)
    cap = 32
    check('packing cap strict impossibility', (cap+1)**2*eta**2 > (cap+1)+(cap+1)*cap*eta**2/2)
    edge = F(D+b,D-b)**n-1
    check('boundary upper bound', 0 < edge < F(3,n*n))
    output = {'stage':stage,'n':n,'exp_terms':terms,'checks':checks,'edge_approx':float(edge),'inputs':[]}
    for label, density_cap, radius_upper in [('curved_uniform',1,radii[0]),('high_stripe',n,radii[1])]:
        arg = F(2*density_cap,1)/(eta**2*theta*c)
        exp_lower = sum((F(radius_upper,1)**k/F(factorial(k),1) for k in range(terms+1)), F(0))
        check(label+' strict log-radius envelope', exp_lower > arg)
        first = F(D,D-b)**n*cap*F((2*radius_upper)**n,1)/(theta*c*factorial(n))
        bound = first+edge
        check(label+' positive coverage bound', first > 0 and bound > 0)
        check(label+' total coverage upper envelope', bound < F(4,n*n))
        if label in previous:
            check(label+' bound strictly decreases across rounds', bound < previous[label])
        previous[label] = bound
        output['inputs'].append({'name':label,'density_cap':density_cap,'radius_upper':radius_upper,'coverage_bound_approx':float(bound),'exact_bound':digest_fraction(bound)})
    # Weak ratio is bracketed by core and support boxes, using c D^n <= W <= D^n.
    lower = tau*F(D-b,D)**n
    upper = tau/c*F(D+b,D)**n
    check('weak-ratio bracket lies near one', 1-F(2,n) < lower < 1 and upper < 1+F(2,n))
    output['weak_ratio_bracket_approx'] = [float(lower),float(upper)]
    output['check_count'] = len(checks)
    records.append(output)

result = {'status':'PASS','total_checks':total,'execution_scope':'Exact scalar parameter and packing/boundary envelopes only; not actual gates, spatial integral, packet optimizer, or general weak endpoint.','records':records}
result_path = BASE/(PREFIX+'_results.json')
result_path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
receipt = {'status':'PASS','executions':1,'checks':total,'registration_sha256':hashlib.sha256((BASE/(PREFIX+'_registration.md')).read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'results_sha256':hashlib.sha256(result_path.read_bytes()).hexdigest()}
(BASE/(PREFIX+'_receipt.json')).write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':total,'rounds':[(r['n'],r['check_count'],[(i['name'],i['coverage_bound_approx']) for i in r['inputs']]) for r in records]},ensure_ascii=False))
