#!/usr/bin/env python3
"""Exact constant-chain audit for the analytic orthogonal-cloud obstruction.
No sampling of high-dimensional volumes or enumeration of exponential subsets.
"""
from fractions import Fraction as F
import json

def main():
    results=[]
    # Three size batches; all quantities here are exact rational/integer.
    for batch,ms in enumerate(((512,1024,2048),(4096,8192,16384),(32768,65536,131072))):
        for m in ms:
            n=512*m;j=(18*m-1).bit_length();p=F(6,2**j)
            near=F(3,2)*m*p;far=1-near;av=F(8,2**j) # alpha*omega_n for R_j=1
            assert near==F(9,32) and far==F(23,32) and p/av==F(3,4)
            good=1-F(128*m,n)-F(64,n)
            assert good>=F(1,2)
            # Bernoulli lower bounds certify all radial power inequalities.
            assert 1+F(n,16)>=18*m
            assert 1+F(n,32)>16*m # R_1^2<17/16 also controls all prefix capture margins
            assert 1+F(n,32)>=F(9*m,8)*F(1000,999)
            assert 1+F(n,14)>F(16,3) # (8/7)^(n/2)>16/3
            assert 1+F(n,7)>F(8,3) # singleton r_min>=7/8
            assert 1+F(n,31)>=2 # 2^(-1/n)>=31/32
            assert (1+F(31,32))**2>F(15,4)
            assert F(9,8)**8>2
            exponent=2*m-n//16
            assert exponent<=-1 # epsilon<=2^exponent<=1/2, no giant powers needed
            # Distances: other source >=sqrt(9/8); good pair in [sqrt(3/2),sqrt(15/4)].
            assert F(7,8)+F(1,2)-F(1,4)==F(9,8)
            assert F(1,2)+F(7,4)-F(1,4)-F(1,2)==F(3,2)
            assert 1+2+F(1,4)+F(1,2)==F(15,4)
            # Fixed record overlap and nonzero theta despite zero common source.
            W_over_alpha=F(1,8);lost=p;gain=p/2
            theta=W_over_alpha*av/lost
            assert theta==F(1,6) and W_over_alpha*av/(lost-gain)==F(1,3)
            # Explicit compact smoothing h=1/(1000n).
            h=F(1,1000*n);power_up=F(875,874)
            assert F(17,16)+F(3,1000)<F(9,8)
            assert power_up*F(3,2)<2 # total near superlevel volume <=2*near/alpha
            assert F(8,7)/power_up>1 and F(8,5)/F(874,875)<2
            assert F(1,16)/F(999,1000)<F(1,8)
            assert F(21,32)*power_up<1 and F(21,64)*power_up<F(1,2)
            assert F(875,874)*F(1000,999)*F(3,2)<2
            assert far-F(1,16)>F(1,2)
            assert F(1,2)-4*F(1,1000)>F(1,4)
            assert F(8,5)/F(874,875)+F(1,1000)<2
            assert near/4+F(1,1000)<F(1,8)
            assert F(998,1000)>far*F(1000,999)/4 # far singleton halo remains inside J=1 exclusion
            atomic_numerator=F(1,2)*good*m*m*(p/2)*(p/4)/av
            smooth_numerator=atomic_numerator/4
            denominator=2*near
            assert atomic_numerator/denominator>=F(m,128)
            assert smooth_numerator/denominator>=F(m,512)
            results.append(dict(batch=batch,m=m,n=n,j=j,D=j+1,near_mass=str(near),far_mass=str(far),good_pair_fraction_lower=str(good),subset_error_log2_upper=exponent,smoothing_radius=str(h),compact_theta='1/6',compact_common_source_mass='0',atomic_lag_ratio_lower=str(F(m,128)),smooth_lag_ratio_lower=str(F(m,512)),smooth_ratio_over_dimension=str(F(1,262144)),positive_density_support_ratio_lower=str(F(m,1024))))
    return dict(status='passed',batch_counts=[3,3,3],cases=results,method='Fraction and integer checks of proved inequalities; no numerical high-dimensional integration',scope='Counters source-mass-deleted lag bounds only; not a maximal weak-type lower bound',main_theorem_proved=False)
if __name__=='__main__':print(json.dumps(main(),indent=2))
