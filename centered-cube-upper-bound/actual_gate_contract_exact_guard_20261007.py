#!/usr/bin/env python3
"""Exact three-round old actual-gate exclusion; no phi/FIRST/history rerun."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent

def rec(x):
    x=F(x)
    if x.numerator.bit_length()+x.denominator.bit_length()<500:
        return f'{x.numerator}/{x.denominator}'
    num=x.numerator.to_bytes((x.numerator.bit_length()+7)//8,'big')
    den=x.denominator.to_bytes((x.denominator.bit_length()+7)//8,'big')
    raw=len(num).to_bytes(8,'big')+num+len(den).to_bytes(8,'big')+den
    return {'sha256_framed_binary':hashlib.sha256(raw).hexdigest(),
            'numerator_bits':x.numerator.bit_length(),'denominator_bits':x.denominator.bit_length()}

def run(n):
    k=(n+1).bit_length() # ceil log2(n+2)
    assert 2**(k-1)<n+2<=2**k
    floor_log=(n+2).bit_length()-1
    # log(n+2)>=floor_log*log2>=floor_log/2, proven log2>=1/2.
    # Hmark/a > (1664/3) log(n+2) >= (832/3)*floor_log.
    H_lower=F(832*floor_log,3)
    assert F(k)<H_lower
    R0=1+F(16,n)
    rows=[]
    for delta in (F(0),F(1,64*n),F(1,32*n),F(1,16*n)):
        R=R0+delta
        qRn=(R/R0)**n/3
        # Bound valid for every h(P)<=W=1 and every CP v_j<=1.
        assert qRn>=F(1,3)
        assert n*qRn*qRn>1
        massrows=[]
        for h in (F(0),F(1,2**n),F(1,3),F(1)):
            assert h*h<n*qRn*qRn
            massrows.append({'h':rec(h),'strict_CP_possible':False,
                             'GP_strict_for_this_h':h>qRn,'exact_PASS':True})
        rows.append({'R_over_R0':rec(R/R0),'qR_to_n':rec(qRn),
                     'squared_CP_threshold_lower':rec(n*qRn*qRn),
                     'all_h_le_W_all_v_le_1_CP_paid':True,'mass_components':massrows})
    return {'n':n,'k':k,'source_l1_diameter_over_a':k,
            'Hmark_over_a_certified_lower':rec(H_lower),
            'l1_far_gate_impossible':True,'CP_global_square_lower':rec(F(n,9)),
            'CP_all_parents_impossible':True,'receiver_box_components':rows,
            'exact_PASS':True}

def main():
    data={'seed':None,'scope':'Only exact l1-diameter/old MM-near and global mass/strict-CP exclusion. All original cube phi/FIRST records frozen; no new actual-history qualification or general endpoint claim.',
          'analytic_prerequisites':['MM Hmark positive constant and log coefficient 1664/3','ln2>=1/2','probe Ls=Rh=R>=R0 and q=R0^-n/3','every eligible h(P)<=W=1 and v_j<=1'],
          'rounds':[run(n) for n in (512,1024,4096)],'all_exact_PASS':True}
    out=HERE/'actual_gate_contract_exact_guard_20261007_results.json'
    out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'rounds':3,'receiver_scale_rows':12,'parent_mass_components':48,
                     'l1_certificates':3,'global_CP_certificates':3,'all_exact_PASS':True,
                     'results':str(out)},ensure_ascii=False))

if __name__=='__main__':main()
