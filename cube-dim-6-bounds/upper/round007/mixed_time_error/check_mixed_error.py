#!/usr/bin/env python3
"""Bounded exact m<=12 constants and ONE bivariate polynomial tensor self-check."""
from fractions import Fraction as F
from math import comb,factorial
from pathlib import Path
from time import perf_counter
import hashlib,json,os
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
ctx.prec=384

def coefficients(M=12):
    q=F(6,5);a=[F(0)]+[F(comb(63+k,k),32**k) for k in range(1,M+1)]
    power=[F(1)]+[F(0)]*M;out=[F(0)]*(M+1);binq=F(1)
    for j in range(1,M+1):
        power=[sum(power[i]*a[k-i] for i in range(k+1)) for k in range(M+1)]
        binq*=F(q-j+1,j)
        for k in range(M+1):out[k]+=abs(binq)*power[k]
    return [out[k]*factorial(k) for k in range(M+1)]

def gauss3_moment(k):
    # Symmetric nodes 1/2 +-sqrt(15)/10; odd powers of sqrt(15) cancel.
    pair=sum(F(comb(k,j)*15**(j//2),2**(k-j)*10**j) for j in range(0,k+1,2))
    return F(5,9)*pair+F(4,9)*F(1,2)**k

def exact_2d_check():
    # F(x,y)=(1+x+y)^12; this is the single analytic self-check.
    I=[F(1,k+1)for k in range(13)];Q=[gauss3_moment(k)for k in range(13)];E=[Q[k]-I[k]for k in range(13)]
    assert all(E[k]==0 for k in range(6))
    II=F(0);QQ=F(0);EI=F(0);IE=F(0);EE=F(0)
    for a in range(13):
      for b in range(13-a):
        c=F(factorial(12),factorial(a)*factorial(b)*factorial(12-a-b))
        II+=c*I[a]*I[b];QQ+=c*Q[a]*Q[b]
        EI+=c*E[a]*I[b];IE+=c*I[a]*E[b];EE+=c*E[a]*E[b]
    beta=F(factorial(3)**4,7*factorial(6)**3)
    assert E[6]==-factorial(6)*beta
    assert EE==factorial(12)*beta**2 and EE>0
    assert QQ-II==EI+IE+EE
    assert QQ-II!=EI+IE
    # Exact full-Lebesgue mixed derivative L1 is 12! for this polynomial.
    c=F(1,2**6*factorial(6));mixed_envelope=c*c*factorial(12)
    assert EE<=mixed_envelope
    return {'function':'(1+x+y)^12','II':str(II),'QQ':str(QQ),'QQ_minus_II':str(QQ-II),'E1I2':str(EI),'I1E2':str(IE),'E1E2':str(EE),'Peano_integral_beta':str(beta),'mixed_derivative_L1':str(factorial(12)),'mixed_error_sup_kernel_envelope':str(mixed_envelope),'full_identity_pass':True,'omitting_mixed_term_is_false':True}

def av(f):return arb(f.numerator)/f.denominator

def main():
    start=perf_counter();D=coefficients();c=F(1,2**6*factorial(6));e=arb(1)/1000;B=arb(1)/16000
    assert D[6]==F(73354040202609,16384000000)
    assert D[12]>D[6]**2
    A=arb(243)**(-arb(1)/6)*arb.pi()**(-arb(11)/6)
    rows=[]
    for s,m in [(1,6),(2,12)]:
      alpha=arb(m)-arb(7)/6
      upper=A*e**(-alpha)/alpha
      lower=upper*(1-2**(-alpha))/(1+4096*e*e)
      Rlo=e**m*lower/B;Rhi=e**m*upper/B
      rho=c**s*D[m]
      # Relative ALLOCATED error at conditional target B and common h=2e.
      factor=av(rho)*2**m
      rows.append({'axes':s,'total_order':m,'D_m_exact':str(D[m]),'Peano_product_sup_exact':str(c**s),'rho_exact':str(rho),'R_lower':str(Rlo),'R_upper':str(Rhi),'per_subset_relative_budget_at_h_2epsilon_lower':str(factor*Rlo),'per_subset_relative_budget_at_h_2epsilon_upper':str(factor*Rhi),'five_axis_subset_multiplicity':comb(5,s)})
    out={'status':'PASS_EXACT_MIXED_ORDER12_AND_ONE_2D_CHECK','pid':os.getpid(),'max_total_order':12,'mixed_D12_over_D6_squared':str(D[12]/D[6]**2),'rows':rows,'self_check':exact_2d_check(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'new_body32_kernel_nodes':0,'matlab_called':False,'elapsed_seconds':perf_counter()-start}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
