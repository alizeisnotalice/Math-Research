#!/usr/bin/env python3
"""Only postprocess registered saved cells; no geometry oracle or sampling."""
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import gzip, hashlib, json, time, sys
sys.set_int_max_str_digits(100000)

BASE=Path(__file__).resolve().parent
PREFIX='nearflat_cube_geometry_probe_20261007'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(q): return Decimal(q.numerator)/Decimal(q.denominator)
def newjson(p,obj):
    with p.open('x') as h: json.dump(obj,h,ensure_ascii=False,indent=2);h.write('\n')

def atanh_log_interval(r,N=80):
    z=(r-1)/(r+1)
    assert 0<=z<=F(1,3)
    partial=2*sum((z**(2*j+1)/F(2*j+1) for j in range(N)),F(0))
    tail=2*z**(2*N+1)/(F(2*N+1)*(1-z*z))
    return partial,partial+tail

@lru_cache(None)
def log_interval(r):
    if r<=1:return F(0),F(0)
    k=0;s=r
    while s>=2:s/=2;k+=1
    lo,hi=atanh_log_interval(s)
    a,b=atanh_log_interval(F(2))
    return lo+k*a,hi+k*b

def phi(hist,tau):
    low=high=F(0);numeric=Decimal(0)
    for M,vol in hist.items():
        if M<=tau:continue
        l,h=log_interval(M/tau)
        low+=tau*vol*l;high+=tau*vol*h
        numeric+=dec(tau*vol)*dec(M/tau).ln()
    return low,high,numeric

def main():
    started=time.perf_counter()
    rp=BASE/f'{PREFIX}_registration.json';reg=json.loads(rp.read_text())
    allrounds=[]
    with localcontext() as ctx:
        ctx.prec=160
        for entry in reg['source_profiles']:
            p=Path(entry['path']);assert sha(p)==entry['sha256']
            saved=json.load(gzip.open(p,'rt'));inp=saved['input'];n=inp['n']
            tau0=F(inp['tau']);tau=tau0/2;J=list(map(F,inp['J']))
            weights=list(map(F,inp['weights']));v=list(map(F,inp['v']))
            atoms=[tuple(map(F,a)) for a in inp['atoms']]
            phys=sorted(set(atoms));idx={a:i for i,a in enumerate(phys)}
            w=[F(0)]*len(phys);signed=[F(0)]*len(phys)
            for y,a,vi in zip(atoms,weights,v):w[idx[y]]+=a;signed[idx[y]]+=a*vi
            vp=[s/a for s,a in zip(signed,w)]
            assert all(abs(vi)<=1 for vi in vp) and sum(signed)==0 and sum(w)==1
            ts=list(map(F,reg['t_positive']));hists={F(0):defaultdict(F)}
            for t in ts:hists[t]=defaultdict(F);hists[-t]=defaultdict(F)
            Evolume=F(0);tievolume=F(0);A=F(0);records=[]
            checks=0;ecount=tiecount=0
            for ci,c in enumerate(saved['cells']):
                vol=F(c['volume']);M=F(c['M']);rr=list(map(F,c['full_response_by_scale']))
                assert M==max(rr)
                hists[F(0)][M]+=vol
                for pert in c['perturbations']:
                    t=F(pert['t']);hists[t][F(pert['M'])]+=vol
                if M<=tau:continue
                Evolume+=vol;ecount+=1
                ties=[j for j,r in enumerate(rr) if r==M];j0,j1=ties[0],ties[-1]
                cap0={idx[atoms[i]] for i in c['full_capture_labels_by_scale'][j0]}
                cap1={idx[atoms[i]] for i in c['full_capture_labels_by_scale'][j1]}
                assert cap0<=cap1
                m0=sum((w[i] for i in cap0),F(0));m1=sum((w[i] for i in cap1),F(0))
                assert m0==F(c['full_mass_by_scale'][j0]) and m1==F(c['full_mass_by_scale'][j1])
                pi0=[w[i]/m0 if i in cap0 else F(0) for i in range(len(w))]
                pi1=[w[i]/m1 if i in cap1 else F(0) for i in range(len(w))]
                assert sum(pi0)==sum(pi1)==1
                beta=m0/m1;assert beta==(J[j0]/J[j1])**n
                dif=[a-b for a,b in zip(pi0,pi1)]
                assert all(dif[i]==(1-beta)*pi0[i] for i in cap0)
                p0=sum((a*vi for a,vi in zip(pi0,vp)),F(0))
                p1=sum((a*vi for a,vi in zip(pi1,vp)),F(0))
                collision=sum((a*a for a in pi0),F(0));normsq=sum((a*a for a in dif),F(0))
                assert normsq>=(1-beta)**2*collision
                A+=tau*vol*abs(p0-p1);checks+=6
                if j0!=j1:tiecount+=1;tievolume+=vol
                records.append({'saved_cell_index':ci,'volume':str(vol),'Rmin':str(J[j0]),'Rmax':str(J[j1]),'physical_atom_capture_min':sorted(cap0),'physical_atom_capture_max':sorted(cap1),'beta':str(beta),'posterior_min':[str(x) for x in pi0],'posterior_max':[str(x) for x in pi1],'delta_p':str(p0-p1),'collision_min':str(collision),'posterior_delta_l2_squared':str(normsq)})
            I=tau*Evolume
            phis={t:phi(h,tau) for t,h in hists.items()};lo0,hi0,num0=phis[F(0)]
            rows=[]
            for t in ts:
                lp,hp,np=phis[t];lm,hm,nm=phis[-t]
                rhs=t*A-t*t*I/(1-t)**2
                slo=lp+lm-2*hi0;shi=hp+hm-2*lo0;nd=np+nm-2*num0
                assert slo>=rhs
                rows.append({'t':str(t),'rational_switch_RHS':str(rhs),'switch_RHS_decimal':str(dec(rhs)),'pair_Phi_delta_decimal':str(nd),'pair_delta_lower_decimal':str(dec(slo)),'pair_delta_upper_decimal':str(dec(shi)),'interval_width_upper_decimal':str(dec(shi-slo)),'certified_lower_bound_residual_decimal':str(dec(slo-rhs)),'switch_inequality_rational_interval_certified':True,'strict_positive_RHS':rhs>0,'strict_positive_pair_delta_interval_certified':slo>0,'pair_delta_lower_rational':str(slo),'pair_delta_upper_rational':str(shi)})
            rr={'n':n,'source_profile_sha256':sha(p),'original_tau':str(tau0),'new_tau':str(tau),'physical_atom_count':len(phys),'physical_v':[str(x) for x in vp],'saved_cell_count':len(saved['cells']),'new_E_cell_count':ecount,'new_E_exact_tie_cell_count':tiecount,'new_E_volume':str(Evolume),'new_E_tie_volume':str(tievolume),'I':str(I),'tau_integral_abs_delta_p':str(A),'exact_posterior_checks':checks,'Phi_base_decimal':str(num0),'variation_pairs':rows}
            allrounds.append(rr)
            cp=BASE/f'{PREFIX}_n{n}_posterior_cells.json.gz'
            payload={'n':n,'tau':str(tau),'source_profile_sha256':sha(p),'cells':records}
            if cp.exists():
                with gzip.open(cp,'rt',encoding='utf8') as h:assert json.load(h)==payload
            else:
                with gzip.open(cp,'xt',encoding='utf8') as h:json.dump(payload,h,separators=(',',':'))
            rr['posterior_cells_path']=str(cp);rr['posterior_cells_sha256']=sha(cp)
            print(json.dumps({'n':n,'E_ties':tiecount,'A':str(A),'positive_RHS_pairs':sum(r['strict_positive_RHS'] for r in rows),'certified_positive_pair_deltas':sum(r['strict_positive_pair_delta_interval_certified'] for r in rows)}),flush=True)
    out={'status':'passed','registration_sha256':sha(rp),'script_sha256':sha(Path(__file__)),'rounds':allrounds,'elapsed_seconds':time.perf_counter()-started,'saved_cells_only':True,'old_oracle_or_MC_rerun':False,'log_interval_method':'Exact Fraction 80-term atanh with rigorous geometric tail; dyadic reduction and sign-aware subtraction. Decimal160 is a separate diagnostic.','nearoptimal_input_or_known_K_claim':False,'dimension_order_claim':False}
    newjson(BASE/f'{PREFIX}_results.json',out)
    print(json.dumps({'status':'passed','seconds':out['elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
