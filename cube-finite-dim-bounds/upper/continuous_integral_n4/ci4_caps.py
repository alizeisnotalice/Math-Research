"""Arb enclosures for the analytic line cap and a crude whole-octant bound.
Only calculation evidence; analytic radial domination is in ci4_report.md.
"""
from flint import arb,ctx
from pathlib import Path
import hashlib,json,argparse
ctx.prec=256
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')
LOCK={'upper/u11_green.py':'fdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d','proofs/round011_face_caps.md':'b1564bf4c0cb4a6e2f72786ea5ffae30c922ad45c0cb9874aae6bfa55d4fa8de','audit/030_audit.md':'dcbc9f31561668e0b31922990b659786038f205cdb106230921ace53a72db640'}
def interval(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display':str(x)}
def run(save):
 for p,h in LOCK.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,(p,'source drift')
 pi=arb.pi();D=1/(48*pi**4);B=arb(5632)*32/pi**3
 def p43(x):return x.root(3)**4
 # avoid precedence ambiguities; r^(5/3)=(cuberoot(r))^5.
 def E(e):return arb(2).root(3)*4*pi*(p43(D)*e**3/3+arb(3)/5*p43(B)*e.root(3)**5)
 eps=arb(1)/10**10;cap=E(eps)
 out={'status':'analytic budget expression intervals, not actual integral intervals','precision_bits':ctx.prec,'source_hashes':LOCK,'epsilon':'1/10000000000','D':interval(D),'B':interval(B),'full_ball_cap_E4':interval(cap),'positive_octant_cap':interval(cap/8),'ordered_wedge_cap':interval(cap/48),'all_faces_caps_refined':interval(8*cap),'I4_crude_global_upper_expression':interval(E(arb(3).sqrt())/8),'target':'19/100000','strict_cap_slack_comparison':bool(cap/8<arb(19)/100000)}
 if save:Path(save).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--save');run(p.parse_args().save)
