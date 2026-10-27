#!/usr/bin/env python3
"""Read-only exact Radau algebra/conditional-cost interfaces, not mesh replay."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
import sympy as sp
P=Path(__file__).resolve().parents[1]
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
def main():
 c=json.loads((P/'upper/round005/radau_mellin/certificate.json').read_text());old=json.loads((P/'upper/round004/pair_mellin/certificate.json').read_text())
 assert c['reused_q_1_8']==old['q']
 lo,hi=dy(c['I_new_upper']['lower']),dy(c['I_new_upper']['upper']);cap=Q(7730422817512,10**18)
 assert 0<lo<=hi<cap<dy(old['I_pair']['lower'])
 assert dy(c['delta_lower']['lower'])>0 and dy(c['delta_lower']['upper'])<1
 mu,nu,G,H=sp.symbols('mu nu G H',positive=True);d=nu/mu
 alpha=(2*(G-1)-d*H)/d;beta=(1-G+d*H)/d**2
 assert sp.simplify(1+alpha*mu+beta*nu-(1+(G-1)*mu**2/nu))==0
 g=Q(1,219000);cost6=Q(7)**7*hi**6/g**7;assert cost6<128**6
 print(json.dumps(dict(status='PASS primary exact Hermite expectation algebra and endpoint comparison',new_C7_upper=False,global_g_proved=False,conditional_g_target=str(g),conditional_cost_sixth_power_exact=str(cost6),conditional_below128=True,I_finer_upper=str(hi),I_rational_upper=str(cap),certificate_sha256=hashlib.sha256((P/'upper/round005/radau_mellin/certificate.json').read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':main()
