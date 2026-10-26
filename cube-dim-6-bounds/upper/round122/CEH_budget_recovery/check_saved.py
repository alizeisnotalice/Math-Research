#!/usr/bin/env python3
"""Exact saved-only verification. No producer import or defining integral run."""
import sys
sys.set_int_max_str_digits(0)
import hashlib,json
from pathlib import Path
from fractions import Fraction as Q
r=Path(__file__).resolve().parent
inp=json.loads((r/'upstream.json').read_text());out=json.loads((r/'results.json').read_text())
assert out['status']=='COMPLETE_FOUR_FIXED_CONFIGURATIONS'
p,c,h,l=out['configurations']
for src,row in zip(inp['peaks033'],p['nine_cases']):
    rank=int(src['dimension']);D=Q(src['Gram_determinant']);R=Q(src['R_exact'])
    sq=(2+Q(2*rank,3))**rank/(16**(2*rank)*D*R**2)
    assert sq==Q(row['ratio_BL_over_exact_squared'])>1
assert p['Crem_times_pi6_exact']=='44743079/2970'
assert c['B_actual_le_p_before_enlargement'] and not c['B_bound_le_p_asserted']
s=list(map(Q,c['squares']));av=c['time_coefficients'];C0=Q(c['C0']);p0=Q(c['p0'])
def poly(u):
    a=Q(1)
    for si,ai in zip(s,av):a*=si+ai*u
    return C0*u**3-p0*a
ul,uu=map(Q,c['u_star_interval']);tl,tu=map(Q,c['T_interval'])
assert 0<ul<uu and poly(ul)<=0<=poly(uu)
assert tl>0 and tl**2<=ul and tu**2>=uu
assert c['root_bisections']==192 and c['partial_fraction_residual']=='0'
def sb(v):return [Q(int(v[z+'_mantissa']))*Q(2)**int(v[z+'_binary_exponent']) for z in ('lower','upper')]
piL,piU=sb(inp['high062']['constants']['pi'])
fl,fu=map(Q,c['finite_primitive_enclosure']);jl,ju=map(Q,c['B_new_time_integral'])
assert jl==C0/piU**6*fl+(p0/piU**6)/(5*tu**5)
assert ju==C0/piL**6*fu+(p0/piL**6)/(5*tl**5)
vp=sb(inp['high062']['constants']['V_power_5_over6']);ps=sb(inp['high062']['Pshort']);pl=sb(inp['high062']['Plong'])
bl,bu=map(Q,h['B_new_norm']);nl,nu=map(Q,h['new_H_norm']);il,iu=map(Q,h['new_H_I'])
assert bl==vp[0]*jl and bu==vp[1]*ju
assert nl==ps[0]+pl[0]+bl and nu==ps[1]+pl[1]+bu
assert il**5<=nl**6 and iu**5>=nu**6
for name in ('C_root','E_root','H_root'):
    z=inp['ledger'][name];assert Q(z['lower'])**int(z['degree'])<=Q(z['radicand'])<=Q(z['upper'])**int(z['degree'])
ce=Q(inp['ledger']['C_root']['upper'])+Q(inp['ledger']['E_root']['upper']);oldH=Q(inp['ledger']['H_root']['upper'])
old=Q(l['old_CEH_I_upper']);new=Q(l['new_CEH_I_upper']);delta=Q(l['I_recovery_exact'])
assert ce+oldH==old==Q(inp['paid_CEH_111'])
assert new==ce+min(oldH,iu) and delta==old-new>=0
rl,ru=map(Q,l['old_far_allocation']);newrl,newru=map(Q,l['new_far_sufficient_allocation'])
assert newrl==rl+delta and newru==ru+delta
assert not l['g0_globally_proved'] and not l['actual_far_paid'] and not l['new_C6_bound']
print(json.dumps({'status':'PASS_SAVED_ONLY','results_sha256':hashlib.sha256((r/'results.json').read_bytes()).hexdigest(),'new_science_configurations':0,'producer_or_old_program_runs':0}))
