#!/usr/bin/env python3
"""Final three actual scalar items; no kernel imports or passes."""
from pathlib import Path
from fractions import Fraction
from flint import arb,fmpq,ctx
import json,datetime,hashlib
ctx.prec=512
O=Path(__file__).resolve().parent; U=O.parents[1]
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round023/partial_fraction_kernel')
F=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/round023/partial_fraction_review')
def read(p):return json.loads(p.read_text())
def A(s):
    f=Fraction(str(s)); return arb(fmpq(f.numerator,f.denominator))
def ball(v):return v.str(120,more=True)
# Read all source fields before any new scalar item.
b=read(B/'results.json')['Gamma_G']; f=read(F/'results.json')['Gamma_G']
tail_text=read(U/'round005/primary/new_tail_budget.json')['conditional_targets'][0]['tail']
Q_text=read(U/'round006/integral/gauss3_relative_results.json')['quadrature_Q_interval']
items=[]
def save(name,status,**data):
    items.append({'name':name,'evidence_status':status,**{k:(ball(v) if isinstance(v,arb) else v) for k,v in data.items()}})
    out={'status':'RECOVERED_PARTIAL_OUTPUT' if len(items)<3 else 'COMPLETE_NO_MODEL_EXCLUSION',
         'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'precision_bits':512,
         'actual_scalar_items_total':12,'first_attempt_actual_items':4,'second_attempt_actual_items':5,'final_attempt_items':len(items),'items':items,
         'kernel_passes':0,'native_calls':0,'spatial_grids':0,'new_time_values':0,'matlab_calls':0,
         'failure_note':'First attempt stopped after4 actual items on missing upper/ path. Second stopped after5 actual items on decimal-to-fmpq parser before any item6 arithmetic. Their unsaved results are not claimed as retained evidence. Only these final3 items are preserved; total actual12. No further scalar run.',
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (O/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(items[-1],indent=2),flush=True)
Gamma=fmpq(f); assert Gamma==fmpq(b)
G=arb(Gamma)/(128*arb.pi()**5)
save('original_item1_Gamma_and_point_G','STRICT_POINT_G_AND_ESSENTIAL_G_UPPER',Gamma_matches=True,G=G,above_target=bool(G>A('1/60000')),below_1_over_55000=bool(G<A('1/55000')))
g0=A('1/60000'); normcap=A('1/11250'); q=A('6/5'); Icap=normcap**q
save('original_item3_target_competition_gate','STRICT_CONDITIONAL_GATE',g_target=g0,norm_threshold=normcap,I_threshold=Icap)
tail=arb(tail_text); Q=arb(Q_text); Hgoal=A('1/14000'); Ngoal=Hgoal+tail
cost=(6*Ngoal/g0)**q; W=Hgoal-Q/Hgoal**A('1/5')
assert cost<64 and W>0
save('original_item12_relaxed_H_goal','CONDITIONAL_NOT_CERTIFIED',H_norm_goal=Hgoal,tail_from_existing_certificate=tail,G_norm_goal=Ngoal,I_upper_goal=Ngoal**q,required_g_for64=3*Ngoal/16,cost=cost,Gauss_absolute_error_allowance=Hgoal**q-Q,mixed_W_allowance=W)
assert len(items)==3
