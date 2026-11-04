from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
import json,os
from pathlib import Path

def ceilf(x):return -(-x.numerator//x.denominator)
def bits(x):
 x=F(x);return {'numerator':abs(x.numerator).bit_length(),'denominator':x.denominator.bit_length()}

def check(job):
 n,D=job;b=ceilf(F(4**n+4)*(D+1))+1;J=b+5;v=F(1,2**J);N=2**(4*J);theta=F(1);a=v**n;g=v*(1-b*v-v**3)/(2*(1+v));h=(2*g)**n;rb=2*g/v;rup=1-(b-2*D)*v
 gap=1-2*g-(rup+1)/2;gaplower=(F(b,2)-D-1)*v
 weights=[theta*v*F(1,2**j) for j in range(1,J+1)];w0=sum(weights,F(0));w1=sum((j*w for j,w in enumerate(weights,1)),F(0))
 m0=theta*v*(1-F(1,2**J));m1=theta*v*(2-F(J+2,2**J));qstar=a/(1+v)**n
 checks={
 'weight_sum_identity':w0==m0,'weighted_index_identity':w1==m1,'baseline_weight_positive':0<w0<v and 1-w0>0,
 'parameters_positive':0<v<1 and a>0 and 0<g and 0<h<a and 0<D*v<1 and 0<rup<1,
 'neighbor_capacity':(a-h)/h<F(1,4**n),
 'witness_radius_exceeds_hole_width':rb>2*g,
 'r_upper_by_nth_power': (1-b*v-v**3)**n < rup**n*(1+v)**n*(1-D*v),
 'gap_lower_bound': gap>=gaplower,
 'gap_lower_strictly_positive':gaplower>0,
 'column_e_bound_strict':gap/(2*4**(n-1)*rup)>D*v,
 'column_conservative_gap_strict':gaplower/(2*4**(n-1)*rup)>D*v,
 'protection_layer':1-F(N,N+2)**n<=F(2*n,N),
 'rho_divided_difference_at_most_four':0<=h<=qstar<F(1,2),
 }
 return {'n':n,'D':str(D),'b':b,'J':J,'theta':str(theta),'N_is_2_power':4*J,'v_denominator_bitlength':v.denominator.bit_length(),'N_bitlength':N.bit_length(),'h_size_bits':bits(h),'gap_size_bits':bits(gap),'checks':checks,'all_strict_checks':all(checks.values())}
if __name__=='__main__':
 jobs=[(n,D) for n in (2,3,4) for D in (F(1,4),F(1,2),F(1),F(2),F(4))]
 with ProcessPoolExecutor(max_workers=15) as pool:results=list(pool.map(check,jobs))
 out={'formula':'a=v^n; v=2^-(b+5); b=ceil((4^n+4)(D+1))+1; J=b+5; N=2^(4J); theta=1; g=v(1-bv-v^3)/(2(1+v)); h=(2g)^n; Rb=2g/v; Rup=1-(b-2D)v; gap=1-2g-(Rup+1)/2', 'rho_bound':'for 0<=x,y<=1/2, |x/(1-x)-y/(1-y)|=|x-y|/((1-x)(1-y))<=4|x-y|', 'cases':results,'all_pass':all(r['all_strict_checks'] for r in results)}
 Path(__file__).with_name('tail-parameters-result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'case_count':len(results),'all_pass':out['all_pass'],'cases':results},indent=2))
