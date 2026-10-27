#!/usr/bin/env python3
"""Independent finite-a scalar audit using the output of independent_n7.py."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json,hashlib,os,time,resource
start=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent
prior=here/'independent_n7.json'; d=json.loads(prior.read_text()); n=7;t=Q(d['t']);s=Q(d['sigma']);U=Q(d['U']);Bl,Bu=map(Q,d['B']);Pl=Q(d['P'][0]);D=Q(d['D'])
a=2**512; alpha=Q(a,a+1);u=Q(1,a);beta=1+1/t;R0=t*s;R=R0+Q(1,2**80)
# 0 <= ell <= 1 on [R0,R], so P(R) >= P(R0)-n log(R/R0) >= P(R0)-n (R-R0)/R0.
assert 1<R0<R<t
assert 0 < 1+1/(2*t)-1/(2*R0) and 1+1/(2*t)-1/(2*R)<=1
pl=Pl-n*(R-R0)/R0
# log(t) <= t-1 supplies a coarse fully rational error bound.
delta=((1+u)**n-1)*(1+n*(t-1))+u**n
E=Bu-D+delta;Aa=(beta+u)**n
S=10**40
assert E>0
x=Aa*E;m=isqrt(x.numerator*S*S//x.denominator);tau=Q(m+1,S)
guard=alpha*beta-(1+s/R+Q(2,a+1));budget=alpha**(-n)*(pl-D)-tau
cost=(1+Aa/tau)*(tau+E)
costup=Q(-(-cost.numerator*S//cost.denominator),S)
checks={'guard_positive':guard>0,'finite_budget_positive':budget>0,'finite_cost_upper_lt_U':costup<U}
assert all(checks.values())
sources=[Path(__file__),prior,here.parents[2]/'upper/round030/compute/n30b_finite.py']
out={'checks':checks,'a':'2^512','R':str(R),'P_guard_lower_method':'P(R0) - n*(R-R0)/R0','cost_upper':str(costup),'strict_U':str(U),'guard_margin_approx':float(guard),'budget_margin_approx':float(budget),'U_minus_cost_upper_approx':float(U-costup),'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'finite_n7.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
