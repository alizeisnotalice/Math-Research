#!/usr/bin/env python3
"""Exact certificate for the BG11 scaled compact-bump family.

With y=sqrt(epsilon)=2^-k, take
  w_epsilon=1+|x|^(-1/2)+(64/y)*1_(y^2,2y^2), E=(0,y^2), alpha=1/4.
Its bump L1 mass is 64y and its witness ratio is
  (3y+66+2sqrt(2))/(y+2).
"""
from fractions import Fraction as F
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent

def main():
    rows=[]
    # Difference from the k=0 lower witness is
    # (60+2sqrt(2))*(1-y)/(3*(y+2)), which is nonnegative for 0<y<=1.
    # The strict base comparison is 45-22sqrt(2)>0, certified by squares.
    assert 45**2-2*22**2==1057
    for k in range(61):
        y=F(1,2**k); epsilon=y*y
        bump_height=F(64,1)/y; bump_length=epsilon
        bump_l1=bump_height*bump_length
        assert bump_l1==64*y and 0<y<=1
        assert (1-y)>=0 and y+2>0
        # Compute the two rational coefficients of A+B*sqrt(2) directly
        # from the actual support length, background height and pure masses.
        cost=epsilon+2*y
        mass_a=3*epsilon+64*y+2*y
        mass_b=2*y
        ratio_a,ratio_b=mass_a/cost,mass_b/cost
        assert ratio_a==(3*y+66)/(y+2) and ratio_b==2/(y+2)
        assert ratio_a-23==20*(1-y)/(y+2)>=0
        assert ratio_b-F(2,3)==F(2,3)*(1-y)/(y+2)>=0
        rows.append({'k':k,'y':str(y),'epsilon':str(epsilon),
                     'bump_height':str(bump_height),'bump_support':'(y^2,2y^2)',
                     'bump_L1_mass':str(bump_l1),
                     'ratio_formula':'(3y+66+2sqrt(2))/(y+2)',
                     'uniform_lower':'(69+2sqrt(2))/3',
                     'difference_numerator':'(60+2sqrt(2))*(1-y)'})
    out={'family':'w_epsilon=1+|x|^(-1/2)+(64/y)*1_(y^2,2y^2), E=(0,y^2), y=2^-k',
         'scope':'exact algebraic certificate for k=0,...,60; no optimality claim',
         'halo_note':'H_(1/4)(0,y^2)=(-y^2,2y^2) modulo endpoints by scaling.',
         'proof_note':'The stated difference is nonnegative; base strictness uses 45^2-2*22^2=1057.',
         'count':len(rows),'rows':rows}
    (ROOT/'scaled_bump_certificate.json').write_text(json.dumps(out,indent=2))
    print('scaled_BG11_cases=%d L1=64y ratio>=uniform_lower PASS' % len(rows))

if __name__=='__main__': main()
