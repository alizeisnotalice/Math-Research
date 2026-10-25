#!/usr/bin/env python3
"""Arb384 arithmetic for a true whole-face H_epsilon norm upper budget."""
from pathlib import Path
from flint import arb,ctx
import json,hashlib
ctx.prec=384
HERE=Path(__file__).resolve().parent
A=arb(243)**(-arb(1)/6)*arb.pi()**(-arb(11)/6)
e=arb(1)/1000
# Integral_0^infty u^(7/6-1)/(1+u^2)du = pi/(2 sin(7pi/12)).
full=A*arb(32)**(-arb(7)/6)*arb.pi()/(2*(7*arb.pi()/12).sin())
# subtract a LOWER bound on the missing [0,e] integral.
subtract=A*(arb(6)/7)*e**(arb(7)/6)/(1+1024*e*e)
B=full-subtract
# Independent elementary majorant: split at t=1/32, use denom>=1 or1024t².
crude=A*arb(32)**(-arb(7)/6)*(arb(6)/7+arb(6)/5)
assert B>0 and B<full and full<crude and B>arb(1)/16000
# A slightly rounded rational upper endpoint, convenient for a ledger.
assert B<arb(7)/5000 and full<arb(7)/5000
out={'status':'PASS_TRUE_WHOLE_FACE_BUDGET_NOT_COMPETITIVE','epsilon':'1/1000','pointwise_time_norm_envelope':'A*t^(1/6)/(1+1024*t^2), A=243^(-1/6)*pi^(-11/6)','global_H_norm_upper_expression':str(B),'strict_rational_H_norm_upper':'7/5000','full_G_norm_budget_expression':str(full),'strict_rational_G_norm_upper':'7/5000','elementary_crosscheck_budget':str(crude),'short_interval_subtracted_lower':str(subtract),'target':'1/16000','target_not_proved':True,'scope':'all y in [0,1]^5, all t>=epsilon; S2 and overlaps included; integral upper budget, not true norm lower bound','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'global_budget_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
