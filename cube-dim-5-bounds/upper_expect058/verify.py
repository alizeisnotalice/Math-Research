#!/usr/bin/env python3
"""Independent rational checks of the emitted finite certificate."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, os, time

H=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    start=time.monotonic()
    r=json.loads((H/'results.json').read_text())
    reg=json.loads((H/'registration.json').read_text())
    assert sha(__file__)==reg['code_sha256']['verify.py']
    for p,h in r['sources'].items(): assert sha(p)==h
    b=json.loads(Path(reg['budget_source']).read_text())
    pl,pu=map(F,b['pi_machin80_interval'])
    bc=F(r['B0_fifth_power_coefficient_of_pi4'])
    gc=F(r['Jgate_fifth_power_coefficient_of_pi7'])
    assert bc==2*F(2027025,256)**5/(5040**4*F(654729075,1024))
    assert gc==F(60**4*32**9,3125*10000**5)
    expected=[bc*pl**4,bc*pu**4,gc*pl**7,gc*pu**7,
              243*bc**2*pl**8,243*bc**2*pu**8]
    assert len(r['power_bounds'])==6
    for i,row in enumerate(r['power_bounds']):
        x=F(row['input']); lo=F(row['lower']); hi=F(row['upper']); n=row['degree']
        assert x==expected[i] and n==([5]*4+[10]*2)[i]
        assert lo**n<=x<hi**n
    blo,bhi=map(F,r['B0_interval']); glo,ghi=map(F,r['Jgate_interval'])
    assert blo==F(r['power_bounds'][0]['lower']) and bhi==F(r['power_bounds'][1]['upper'])
    assert glo==F(r['power_bounds'][2]['lower']) and ghi==F(r['power_bounds'][3]['upper'])
    assert blo-ghi==F(r['B0_lower_minus_Jgate_upper'])>0
    assert bc-gc*pu**3==F(r['fifth_power_gap_coefficient_lower'])>0
    assert F(r['locked_Jflag_upper']) < F(r['global_Jensen_cap_expression_interval'][0])
    result={'status':'PASS_exact_rational_replay','pid':os.getpid(),
            'seconds':time.monotonic()-start,'results_sha256':sha(H/'results.json'),
            'checks':['source hashes','six outward integer-root brackets',
                      'independent Gamma rational formula','strict fifth-power gate separation',
                      'strict endpoint gap','new crude cap weaker than locked Jflag'],
            'scope':'arithmetic certificate; analytic proof is in report.md'}
    out=H/'verification.json'
    if out.exists():
        print(json.dumps({**result,'prior_verification_left_unchanged':True})); return
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
