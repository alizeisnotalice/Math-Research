#!/usr/bin/env python3
"""Read-only extraction arithmetic; does not test the accepted lower theorem."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
from math import isqrt,log,pi,sqrt
from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
PAPER=HERE.parent.parent/'hl_cube_paper_20261013/paper.tex'

def ceilfrac(x):return -(-x.numerator//x.denominator)
def decimal(x):return Decimal(x.numerator)/Decimal(x.denominator)
def number(x):
    return float(x) if abs(x)<F(10**300) else str(x)
def row(n,p,w):
    with localcontext()as ctx:
        ctx.prec=80;lnn=Decimal(n).ln();threshold=16*lnn
        cp2=F(9,29)*p*p
        d=1
        while decimal(cp2*n/F((2**(d+1)-1)**2))>=threshold:d+=1
        K=2**d-1;v=cp2*n/F(K*K)
        assert decimal(v)>=threshold and decimal(cp2*n/F((2*K+1)**2))<threshold
        B=1
        while B<v:B*=2
        M=isqrt(v.numerator//v.denominator)
        if F(M*M)<v:M+=1
        delta=1/(w*K*v);H=sum((F(1,j)for j in range(1,M+1)),F(0))
        N=B-1;T=N*(N+1)//2;RB=N*(N-1);HB=N+T*(T-1)+2*RB*RB+4*RB
        positivity_error=F(20,3)*p*H/K
        etaD=F(1,3*B*B)+F(2,K)
        etaN=F(2,3*B*B)+F(2,K)
        cx=40/9*((1-1/B)*sqrt(K+1)/K*sum(1/j**2 for j in range(1,M+1))
                +pi*(1-1/B)**2/(2*K)*float(H))
        e=cx+float(F(2,3*w)/v)
        tm=F(48,29*M**3)
        vmin=F(100,261)-tm/3-etaD-F.from_float(e)
        Q=ceilfrac(6*w*K*v)
        bits=0
        while F(1,2**bits)>delta/64:bits+=1
        return dict(n=str(n),p=str(p),w=w,d=d,K=str(K),v=number(v),B=B,b=B.bit_length()-1,M=M,
            v_exact=str(v),delta_exact=str(delta),h_exact=str(v/2-2),
            delta=number(delta),h=number(v/2-2),lambda_log_h=number(v/2-2),
            sixteen_log_n=float(threshold),maximal_K_rule=True,
            gates=dict(B_ge_v=B>=v,Kplus1_gt_HB=K+1>HB,Kplus1_gt_3B4=K+1>3*B**4,
                M_lt_K=M<K,positive_lower_at_least_half_background=positivity_error<=(1-p)/2,
                box_inside_cell=delta<F(1,2*B*K),Vmin_ge_quarter_diagnostic=vmin>=F(1,4)),
            density_triangle_lower=number(1-p-positivity_error),
            amplitude_bound=number(F(77,15)*p/K),amplitude_times_v=number(F(77,15)*p*v/K),
            v_squared_over_n=number(v*v/n),v_etaD=number(v*etaD),v_etaN=number(v*etaN),
            Cx_diagnostic=cx,e_diagnostic=e,Vmin_diagnostic=number(vmin),
            field_threshold_HB=str(HB),K_over_3B4=number(F(K+1,3*B**4)),
            one_K_label_uint16_bytes=str(2*K),ratio_set_uint64_bytes=8*RB,
            field_element_storage_bytes=8*((d+63)//64),
            ratio_set_field_payload_bytes=8*((d+63)//64)*RB,
            max_query_bin_width_delta_over6_bins=str(Q),one_bin_counts_uint64_bytes=str(8*Q),
            query_phase_bits_for_mesh_le_delta_over64=bits,
            pair_label_states=B*B,
            basic_full_count_DP_states_log10=n*log(2) if n<1000 else 'exponential in n; not allocated',
            full_local_grid_nodes_per_unit_t_Kv=number(K*v),
            validity='Algebraic construction gates only. No qualified basis generated, no safe nodes built, no explicit asymptotic n0 or C_n certificate.')

def sampler_probe():
    path=HERE.parent.parent/'hl_cube_paper_20261013/numerics/inputs/source_d10.npz'
    labels=np.load(path)['labels'];K=len(labels);B=8;v=F(1024);w=64;p=F(1,4);M=32
    atoms=(np.arange(K)+(labels+.5)/B)/K;delta=float(1/(w*K*v))
    m=np.arange(1,M+1,dtype=float);alpha=729/29
    bm=4*alpha*((np.sin(np.pi*m*10/27)/(np.pi*m))**4+3*(np.sin(np.pi*m/9)/(np.pi*m))**4)
    cm=np.sqrt(2*bm/3);coeff=np.pi*m*cm
    HM=sum((F(1,j)for j in range(1,M+1)),F(0));C=F(20,3)*p*HM/K
    hmin=1-p-C;hmax=1-p+C;assert hmin>0
    rng=np.random.default_rng(450045);N=100000
    peak=rng.random(N)<float(p);draw=np.empty(N)
    draw[peak]=atoms[rng.integers(K,size=int(peak.sum()))]+(rng.random(int(peak.sum()))-.5)*delta
    needed=int((~peak).sum());accepted=[];proposals=0;accepted_total=0
    while accepted_total<needed:
        batch=max(1024,int((needed-accepted_total)*1.03)+100)
        x=rng.random(batch);h=1-float(p)+float(p)/K*(np.cos(2*np.pi*x[:,None]*m)@coeff)
        ok=rng.random(batch)<h/float(hmax);accepted.append(x[ok]);accepted_total+=int(ok.sum());proposals+=batch
    draw[~peak]=np.concatenate(accepted)[:needed]
    grid=np.linspace(0,1,257)
    peakcdf=np.clip((grid[:,None]-(atoms-delta/2))/delta,0,1).sum(axis=1)/K
    cdf=(1-float(p))*grid+float(p)*peakcdf+float(p)/(2*K)*(np.sin(2*np.pi*grid[:,None]*m)@cm)
    empirical=np.searchsorted(np.sort(draw),grid,side='right')/N
    return dict(K=K,B=B,v=str(v),p=str(p),w=w,M=M,delta=delta,
        label_input_SHA256=hashlib.sha256(path.read_bytes()).hexdigest(),
        mu_mass_exact='1',DM_zero_mode='absent; frequencies1..32',DM_integral_exact='0',rho_mass_exact='1',
        h_triangle_lower_exact=str(hmin),h_triangle_upper_exact=str(hmax),
        h_triangle_lower=float(hmin),h_triangle_upper=float(hmax),
        ideal_background_acceptance=float((1-p)/hmax),
        N=N,seed=450045,peak_draws=int(peak.sum()),background_needed=needed,background_proposals=proposals,
        background_accepted_before_truncation=accepted_total,CDF_grid_points=len(grid),
        CDF_endpoints=[float(cdf[0]),float(cdf[-1])],max_empirical_CDF_grid_error=float(abs(cdf-empirical).max()),
        qualification='Finite legitimate source sampler/CDF diagnostic only. B8<v1024: not matched growth family. Float trig and PRNG are not interval/randomness certificates. These are rho draws, NOT Haar query draws.')

def main():
    rows=[row(10**power,p,64)for p in (F(1,2),F(99,100))for power in (6,12,18,24,30,32,34,35,36,42,60)]
    pilots=[]
    for n,K,B,p,v in [(6704800,255,4,F(1,2),F(8)),(55249289216,1023,8,F(1,4),F(1024)),(885285273600,4095,8,F(1,4),F(1024))]:
        exactv=F(9,29)*p*p*n/(K*K)
        pilots.append(dict(n=n,K=K,B=B,p=str(p),v_exact=str(exactv),
            B_ge_v=B>=exactv,v_ge_16logn=float(exactv)>=16*log(n),
            next_K_fails=float(F(9,29)*p*p*n/((2*K+1)**2))<16*log(n),
            B_required_min_power_two=2**ceilfrac(F.from_float(log(float(exactv),2)))) )
    C=87/(220*sqrt(58*pi))
    out=dict(source_paper_SHA256=hashlib.sha256(PAPER.read_bytes()).hexdigest(),rows=rows,old_pilots=pilots,source_sampler_probe=sampler_probe(),
        constant=dict(Cstar_float=C,one_over35=1/35,p_min_for_pCstar_above_one_over35=(1/35)/C,
            exact_comparison_margin=str(F(87,2992)-F(1,35))),
        computation_scope='Actual local Python Decimal80/Fraction parameter arithmetic and cost estimates, not source instantiation or theorem re-audit')
    (HERE/'PARAMETER_PROBE.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
