"""Read-only postprocessing of frozen PB run plus independent tiny Fraction check."""
from pathlib import Path
from fractions import Fraction as F
import itertools, json, math, hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
rows=json.loads((HERE/'results.json').read_text())['rows']
checks=[]
for C in [0.5606,1.0]:
    max_ks_ratio=max_bin_ratio=0.
    for row in rows:
        eps=C*row['Lyapunov_ratio']
        assert row['Kolmogorov_jump_endpoint_error'] <= eps+1e-13
        max_ks_ratio=max(max_ks_ratio,row['Kolmogorov_jump_endpoint_error']/eps)
        for b in row['bins']:
            budget=b['finite_CDF_endpoints']*eps
            assert b['absolute_error'] <= budget+1e-13
            if budget: max_bin_ratio=max(max_bin_ratio,b['absolute_error']/budget)
    checks.append({'C':C,'all_57_KS_and_456_bin_checks_pass':True,
                   'max_observed_KS_over_bound':max_ks_ratio,
                   'max_observed_bin_error_over_bound':max_bin_ratio})

# Fraction verification of the macro p range, independent finite configurations.
macro_cases=0
for eta in [F(2,5),F(1,2),F(3,5)]:
    for q in [1,2,5]:
        for eplus in [0,eta/3,eta]:
            for eminus in [0,(1-eta)/3,1-eta]:
                p=(q*eta+eplus)/(q+eplus+eminus)
                assert eta/(2-eta) <= p <= 2*eta/(1+eta)
                macro_cases+=1

# Three independent hidden bits; four strongly correlated kernels per coordinate.
# Two tables have complementary permitted count segments; all arithmetic exact.
n=3
z=[[[[F(1)+F((i+b+1)*(1 if h else -1),20),
        F(3,4)+F((i+2*b+1)*(1 if h else -1),30)]
      for b in range(2)] for h in range(2)] for i in range(n)]
def table_poly(kernels,b):
    return sum((math.prod(kernels[i][b][bits[i]] for i in range(n))
                for bits in itertools.product([0,1],repeat=n)
                if (sum(bits)<=1 if b==0 else sum(bits)>=2)),F(0))
values=[]
for hidden in itertools.product([0,1],repeat=n):
    kernels=[z[i][hidden[i]] for i in range(n)]
    values.append(sum((table_poly(kernels,b) for b in range(2)),F(0)))
mean=sum(values)/len(values)
var=sum((v-mean)**2 for v in values)/len(values)
mu=[]; rel=[]
for i in range(n):
    coord=[]; rels=[]
    for b in range(2):
        coord.append([])
        for s in [0,1]:
            m=(z[i][0][b][s]+z[i][1][b][s])/2
            v=sum((z[i][h][b][s]-m)**2 for h in [0,1])/2
            coord[b].append(m);rels.append(v/m**2)
    mu.append(coord);rel.append(max(rels))
predicted=sum((table_poly(mu,b) for b in range(2)),F(0))
cap=math.prod(1+t for t in rel)-1
assert mean==predicted and var<=cap*mean**2

# Independently verify all saved rational invariants, without rerunning their queries.
other=HERE.parent/'round038_conditional'
saved=json.loads((other/'RESULTS.json').read_text())
layer=json.loads((other/'LAYERED_RESULTS.json').read_text())
assert hashlib.sha256((other/'RESULTS.json').read_bytes()).hexdigest()==layer['input_SHA256']
for case,lc in zip(saved['cases'],layer['cases']):
    assert (case['n'],case['R'],case['L'])==(lc['n'],lc['R'],lc['L'])
    m=F(case['conditional_mean']);v=F(case['conditional_variance'])
    relv=F(case['relative_variance']);capv=F(case['relative_variance_cap'])
    assert v==relv*m*m and relv<=capv
    assert sum(F(s['prior'])for s in case['segments'])==1
    table_m=list(map(F,lc['table_gated_conditional_means']))
    table_t=list(map(F,lc['table_relative_variance_caps']))
    assert sum(table_m)==m
    low=[];up=[];scale=10**40
    for t in table_t:
        a=F(math.isqrt(t.numerator*scale**2//t.denominator),scale)
        b=a if a*a==t else a+F(1,scale)
        assert a*a<=t<=b*b
        low.append(a);up.append(b)
    lower=(sum(w*a for w,a in zip(table_m,low))/m)**2
    upper=(sum(w*a for w,a in zip(table_m,up))/m)**2
    assert lower==F(lc['layered_cap_lower_rational'])
    assert upper==F(lc['layered_cap_upper_rational'])
    assert relv<=upper<=capv

rarest=min(((b['PB_probability'],row['n'],row['label'],b['b'],
             .5606*row['Lyapunov_ratio']*b['finite_CDF_endpoints'])
            for row in rows for b in row['bins'] if b['PB_probability']>0))
out={'scope':'zero new tail MC; frozen DP postprocessing and 8-point independent Fraction polynomial',
     'numpy':np.__version__,'BE_checks':checks,'macro_range_Fraction_cases':macro_cases,
     'tiny_polynomial':{'n':n,'hidden_queries':len(values),'mean':str(mean),
                       'variance':str(var),'relative_variance':str(var/mean**2),
                       'cap':str(cap),'mean_substitution_exact':True,'cap_valid_exact':True},
     'saved_conditional_cases_checked':len(saved['cases']),
     'rarest_positive_saved_bin':{'probability':rarest[0],'n':rarest[1],
                                'label':rarest[2],'b':rarest[3],'absolute_BE_budget':rarest[4]},
     'DP_results_SHA256':hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest(),
     'conditional_results_SHA256':layer['input_SHA256']}
path=HERE/'checks.json'
if path.exists():raise RuntimeError('refusing replacement')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
