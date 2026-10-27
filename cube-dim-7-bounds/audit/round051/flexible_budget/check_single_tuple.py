#!/usr/bin/env python3
"""Independent one-tuple Fraction check; no grid, LP, or production rerun."""
import hashlib
import json
import platform
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
GENERIC = Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md')
LEDGER = Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round052/compute/results.json')
EXPECTED = 'aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6'
assert hashlib.sha256(GENERIC.read_bytes()).hexdigest()==EXPECTED

a,s,d = Q(73,8),Q(293,16),Q(1023,2**20)
ej,ea,ed,es,et = Q(1,1024),Q(1,64),Q(1,1024),Q(1,1024),Q(1,1024)
p = Q(3579,4096)
ulo = p/(2*(1+ej))
uhi = 1/(p*(1-d/2)*(1-ej))
mass = (1-ea)*d/(2*(1+3*d/4))
degree = max(Q(4),2*uhi/(1-3*d/4))
bracket = (1-3/a)*ulo-es-et
Gamma = mass*bracket/degree
gamma = (1-ed)*Gamma
sigma = s*d*d
penalty = sigma/(1-sigma)
craw = gamma-penalty
domain = {'a_gt_3':a>3,'s_gt_2a':s>2*a,'d_positive':d>0,'d_le_1_over_16':d<=Q(1,16),
    'd_lt_1_over_32a':d<Q(1,32)/a,'sigma_lt_half':sigma<Q(1,2),
    'epsJ_in_unit':0<ej<1,'epsA_in_unit':0<ea<1,'epsD_in_unit':0<ed<1,
    'epsS_positive':es>0,'epsT_positive':et>0,'bracket_positive':bracket>0,
    'source_cap_exponent_positive':d/32-a*d*d>0,
    'tau_exponent_positive':(s/2-a)*d*d>0}
assert all(domain.values())
assert degree==4 and craw>Q(1,2**16)
assert craw==Q(6468305130094850733787980147,370780099533215320770594052505600)

record=json.loads(LEDGER.read_text())['top64'][0]
assert (Q(record['a']),Q(record['s']),Q(record['d']))==(a,s,d)
computed={'u_lo':ulo,'u_hi':uhi,'mass':mass,'M':degree,'bracket':bracket,
    'delta_lim':Gamma,'delta':gamma,'sigma':sigma,'penalty':penalty,'c_raw':craw}
assert all(Q(record[key])==value for key,value in computed.items())

# Worst local endpoints corroborate the analytic monotonicity bounds in the audit.
rmin,dmax=Q(1,2),Q(1,16)
local_checks={
    'fR_coefficient_below_8':2/(2*rmin-dmax)**2<8,
    'fg_below_4':2*rmin/(2*rmin-dmax)**2<4,
    'f_below_2':2*rmin/(2*rmin-dmax)<2,
    'Jprime_bulk_coefficient_below_32':Q(24)+Q(4,3)<32,
    'Jprime_tail_coefficient_below_32':1+4*dmax<32,
    'matching_log_error_exact':1280*Q(1,2**20)==Q(5,4096),
    'product_lower_exact':1-(Q(1,8)+Q(5,4096))==p,
}
assert all(local_checks.values())
files=[GENERIC,LEDGER,HERE/'proof_review.md',Path(__file__).resolve()]
receipt={'status':'PASS','python':platform.python_version(),'arithmetic':'Fraction',
    'tuple_count_checked':1,'full_grid_repeated':False,'new_LP_solved':False,
    'tuple':{'a':str(a),'s':str(s),'d':str(d),'epsJ':str(ej),'epsA':str(ea),
        'epsD':str(ed),'epsS':str(es),'epsT':str(et)},
    'domain':domain,'computed':{k:str(v) for k,v in computed.items()},
    'dyadic_c':str(Q(1,2**16)),'strict_c_margin':str(craw-Q(1,2**16)),
    'local_constant_checks':local_checks,
    'generic_source_sha256':EXPECTED,
    'sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
    'scope':'One exact fixed tuple and rational local constants. Generic operator validity is established by the independent analytic review, not by this computation. No finite n threshold, grid completeness, continuous optimum or production resource claim is checked.'}
(HERE/'single_tuple_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'c_raw':str(craw),'dyadic_c':receipt['dyadic_c'],
    'strict_c_margin':receipt['strict_c_margin'],'domain':domain},indent=2))
