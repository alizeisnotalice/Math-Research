from pathlib import Path
from flint import arb,ctx
import json,hashlib,math
ctx.prec=384
here=Path(__file__).resolve().parent;src=here.parent/'local_box/arb384_results.json'
d=json.loads(src.read_text());rows=sorted(d['outputs'],key=lambda x:x['index']);assert [x['index'] for x in rows]==list(range(5))
h=arb(1)/1000;e=arb(1)/100;q=arb(6)/5;L=arb(3)/100000
loss=((arb(129)/128).sqrt()*h/e).cosh()**5
mean=arb(0)
for k,row in enumerate(rows):
 assert row['epsilon']=='1/100' and row['denominator']==1000
 assert sorted(row['point_numerators'])==sorted([1000]+[751]*k+[749]*(4-k)+[1])
 value=arb(row['H_epsilon_arb_384']);assert value>L*loss
 mean+=math.comb(4,k)*value**q/16
upper_expression=(2*h)**4*h*(q*5*h*h/(32*e*e)).exp()*mean
assert upper_expression<arb(7)/10**20
# Conditional global target, not attained by this box.
eI=arb(1)/1000;K=(arb(6)/7)*(arb(32)/243)**(arb(1)/6)*arb.pi()**(-arb(11)/6)
tail=K*eI**(arb(7)/6);target_g=arb(1)/60000;target_Hnorm=arb(1)/16000
cost=6**q*((target_Hnorm+tail)/target_g)**q
assert cost<64
out={'status':'PASS_local_continuous_box_and_conditional_global_targets','box':'x1=1; x2..x5 in [749/1000,751/1000]; x6 in [0,1/1000]','epsilon':'1/100','local_G_lower_strict':'3/100000','local_H_power_integral_upper_strict':'7/10^20','H_power_integral_bound_expression':str(upper_expression),'ODE_loss':str(loss),'conditional_global_targets':{'g_ge':'1/60000','H_1_over1000_L6_over5_norm_le':'1/16000','tail_norm_budget':str(tail),'conditional_cost':str(cost),'targets_achieved':False},'scope':'local box only; I bound is for H_epsilon, not G. Global targets remain unproved; no new C6 bound.','input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
