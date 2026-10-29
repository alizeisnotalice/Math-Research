#!/usr/bin/env python3
"""Deterministic one-dimensional periodic class-channel experiment.
Only auxiliary channel Grams; never computes or estimates C_n.
"""
from pathlib import Path
import json,sys,time,platform,hashlib
from fractions import Fraction
import numpy as np
from numpy.polynomial.legendre import leggauss
OUT=Path(__file__).resolve().parent

def source(K,kind,shift=0.,eps=.15,width=.22):
    # Each (class,left,right,density) component belongs to the original cell.
    c=[]
    for i in range(K):
        l=i/K;r=(i+1)/K
        if kind=='uniform': c.append((i,l,r,1.));continue
        c.append((i,l,r,eps))
        phase=(((7*i)%K)/K+shift)%1
        left=(i+(1-width)*phase)/K
        c.append((i,left,left+width/K,(1-eps)/width))
    return c

def integral(c,K,y):
    y=np.asarray(y);base=np.floor(y);t=y-base
    ans=np.broadcast_to(base[:,None]/K,(len(y),K)).copy()
    for i,l,r,d in c: ans[:,i]+=d*np.clip(t-l,0,r-l)
    return ans

def nodes(c,ell,order):
    knots=sorted(set([0.,1.]+[(z+s*ell/2)%1 for _,l,r,_ in c for z in (l,r) for s in (-1,1)]))
    knots=np.asarray(knots);g,w=leggauss(order)
    widths=np.diff(knots);mid=(knots[1:]+knots[:-1])/2
    return (mid[:,None]+widths[:,None]*g/2).ravel(),(widths[:,None]*np.broadcast_to(w,(len(widths),order))/2).ravel()

def channel(K,q,ell,kind,shift=0.,order=16,midpoint=None):
    c=source(K,kind,shift)
    if midpoint:
        x=(np.arange(midpoint)+.5)/midpoint;w=np.full(midpoint,1/midpoint)
    else:x,w=nodes(c,ell,order)
    r=integral(c,K,x+ell/2)-integral(c,K,x-ell/2)
    rt=r.sum(axis=1);L=q+ell;a=(q+rt)/L
    z=(q/K+r)/L;p=z/a[:,None]
    H=(p.T*(w*a))@p*K
    u=np.ones(K)/np.sqrt(K);A=(H-np.outer(u,u));A=(A+A.T)/2
    # Canonical basis of u-perp. Common class labels are the original cells.
    B=np.linalg.qr(np.column_stack((u,np.eye(K)[:,:K-1])))[0][:,1:]
    Ac=B.T@A@B;ev=np.linalg.eigvalsh(Ac)
    full=float(np.sum(w*rt*(1-rt)/(L*(q+rt))))
    return {'K':K,'q':q,'ell':ell,'kind':kind,'shift':shift,'quadrature_nodes':len(x),
      'mass_residual':float(abs(np.sum(w*a)-1)),
      'prior_residual':float(np.max(abs(np.sum((w*a)[:,None]*p,axis=0)-1/K))),
      'trivial_eigen_residual':float(np.max(abs(H@u-u))),
      'trace_nontrivial':float(np.trace(Ac)),'eigenvalues_nontrivial':ev.tolist(),
      'lambda_max':float(ev[-1]),'lambda_min':float(ev[0]),'full_source_chi2':full,
      'jensen_bound':ell*(1-ell)/L**2,'q_uniform_bound':1/(4*q*(q+1)),
      'histogram_trace_infinite_product':float(np.prod(1/(1-ev))),
      'histogram_trace_from_trace':float(1/(1-np.trace(Ac)))},Ac,H

def common_domination(mats):
    mean=np.mean(mats,axis=0);ev,V=np.linalg.eigh(mean)
    floor=max(ev[-1]*1e-12,1e-16);reg=V@np.diag(np.maximum(ev,floor))@V.T
    R=V@np.diag(1/np.sqrt(np.maximum(ev,floor)))@V.T
    gamma=max(float(np.linalg.eigvalsh(R@a@R)[-1]) for a in mats)*(1+1e-8)
    D=gamma*reg+1e-10*np.eye(len(ev))
    return {'count':len(mats),'gamma_over_mean':gamma,'trace_D':float(np.trace(D)),
      'minimum_psd_margin':min(float(np.linalg.eigvalsh(D-a)[0]) for a in mats),
      'isotropic_trace':len(ev)*max(float(np.linalg.eigvalsh(a)[-1]) for a in mats),
      'sum_trace':sum(float(np.trace(a)) for a in mats),
      'histogram_trace_bound_if_trace_D_lt_1':float(1/(1-np.trace(D))) if np.trace(D)<1 else None,
      'D_eigenvalues':np.linalg.eigvalsh(D).tolist(),'D_matrix_in_canonical_u_perp_basis':D.tolist(),
      'max_commutator_frobenius':max(float(np.linalg.norm(a@b-b@a)) for a in mats[:6] for b in mats[:6])}

def hist_checks(H):
    K=len(H);cols=[]
    for i in range(K):
        for j in range(i,K):
            col=np.zeros(K*K);col[i*K+j]=1 if i==j else 1/np.sqrt(2)
            if i!=j:col[j*K+i]=1/np.sqrt(2)
            cols.append(col)
    S=np.stack(cols,axis=1);direct=float(np.trace(S.T@np.kron(H,H)@S))
    lam=np.linalg.eigvalsh(H);formula=float((lam.sum()**2+(lam**2).sum())/2)
    return {'n':2,'direct_symmetric_tensor_trace':direct,'h2_eigenvalue_formula':formula,'residual':abs(direct-formula)}

def main():
    start=time.time();rows=[];common=[];allq={};refH=None
    ells=np.unique(np.r_[np.linspace(.01,.99,41),1/6,1/3,.5,2/3])
    for kind in ['uniform','positive_box']:
      for K in [4,8,16,32]:
       for q in [1,2,3]:
        across=[]
        for ell in ells:
          mats=[]
          for shift in [0.,1/32,1/16,1/8,1/4,.5]:
            row,A,H=channel(K,q,float(ell),kind,shift)
            rows.append(row);mats.append(A);across.append(A)
            if refH is None:refH=H
          d=common_domination(mats);d.update(K=K,q=q,ell=float(ell),kind=kind);common.append(d)
        d=common_domination(across);d.update(K=K,q=q,kind=kind,ell_range='all scanned');allq[f'{kind}_K{K}_q{q}']=d
    # Independent uniform analytic n=1 case: K=2,q=1,ell=1/2 => nontrivial eigenvalue 1/27.
    exact,A,H=channel(2,1,.5,'uniform',order=32)
    exact['exact_trace_fraction']=str(Fraction(1,27));exact['exact_residual']=abs(exact['trace_nontrivial']-1/27)
    # Piecewise quadrature refinement and independent midpoint integration.
    refine=[]
    for K in [4,32]:
      for kind in ['uniform','positive_box']:
       for ell in [1/6,1/3,.509,.99]:
        triples=[channel(K,1,ell,kind,1/16,order=o) for o in [8,16,32]]
        vals=[v[0] for v in triples]
        midtuple=channel(K,1,ell,kind,1/16,midpoint=131072);mid=midtuple[0]
        refine.append({'K':K,'kind':kind,'ell':ell,'trace_orders_8_16_32':[r['trace_nontrivial'] for r in vals],
          'max_order_spread':max(r['trace_nontrivial'] for r in vals)-min(r['trace_nontrivial'] for r in vals),
          'midpoint_N':131072,'midpoint_trace':mid['trace_nontrivial'],
          'midpoint_minus_order32':mid['trace_nontrivial']-vals[-1]['trace_nontrivial'],
          'gram_order8_minus_order32_operator_norm':float(np.linalg.norm(triples[0][1]-triples[-1][1],2)),
          'gram_midpoint_minus_order32_operator_norm':float(np.linalg.norm(midtuple[1]-triples[-1][1],2))})
    summary={'scan_rows':len(rows),'max_mass_residual':max(r['mass_residual'] for r in rows),
      'max_prior_residual':max(r['prior_residual'] for r in rows),'max_trivial_eigen_residual':max(r['trivial_eigen_residual'] for r in rows),
      'minimum_nontrivial_eigenvalue':min(r['lambda_min'] for r in rows),
      'max_trace':max(r['trace_nontrivial'] for r in rows),'max_full_minus_jensen':max(r['full_source_chi2']-r['jensen_bound'] for r in rows),
      'max_class_minus_full':max(r['trace_nontrivial']-r['full_source_chi2'] for r in rows),
      'max_common_fixed_ell_trace_D':max(r['trace_D'] for r in common),
      'max_common_all_ell_trace_D':max(r['trace_D'] for r in allq.values()),
      'max_refinement_spread':max(r['max_order_spread'] for r in refine),
      'max_midpoint_discrepancy':max(abs(r['midpoint_minus_order32']) for r in refine)}
    worst=max(rows,key=lambda r:r['trace_nontrivial']);coef=np.zeros(101);coef[0]=1
    for eigen in worst['eigenvalues_nontrivial']:
      for degree in range(1,101):coef[degree]+=eigen*coef[degree-1]
    finite_hist={str(n):float(coef[:n+1].sum()) for n in [1,2,10,100]}
    payload={'finite_histogram_trace_worst_row':finite_hist,'scope':'Auxiliary class channel; no C_n inference','runtime':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'seconds':time.time()-start},
      'source_spec':{'period':1,'class_cells':'[i/K,(i+1)/K)','class_mass':'1/K','box_baseline':.15,'box_width_fraction':.22,
        'box_left_cell_phase':'(1-.22)*(((7*i mod K)/K+shift) mod 1)','shifts':[0,1/32,1/16,1/8,1/4,.5]},
      'summary':summary,'n1_exact_check':exact,'n2_histogram_check':hist_checks(refH),
      'refinement':refine,'common_psd_fixed_ell':common,'common_psd_all_scanned_ell':allq,'scan':rows}
    (OUT/'RESULTS.json').write_text(json.dumps(payload,indent=2))
    write_summary(payload)
    print(json.dumps(summary,indent=2))
def write_summary(payload):
    rows=payload['scan']
    def compact_common(d):
        return {k:v for k,v in d.items() if k not in ['D_matrix_in_canonical_u_perp_basis','D_eigenvalues']}
    reps=[r for r in rows if r['q']==1 and r['ell']==1/3 and r['shift']==0]
    short={'scope':payload['scope'],
      'status':'Deterministic exploratory integration, not preregistered or blinded; no Monte Carlo confidence claims.',
      'chronology':'Definitions and candidate Jensen algebra were considered before execution; full scans and refinements ran before the formal proof was written after inspecting results.',
      'probe_SHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'full_local_results_filename':'RESULTS.json (not required for compact archive; regenerate with probe.py)',
      'runtime':payload['runtime'],'source_spec':payload['source_spec'],
      'configuration':{'K':[4,8,16,32],'q':[1,2,3],
        'ell':sorted(set(r['ell'] for r in rows)),'kinds':['uniform','positive_box'],
        'quadrature_order':16,'refinement_orders':[8,16,32],
        'independent_midpoint_N':131072,'row_count':len(rows)},
      'summary':payload['summary'],'n1_exact_check':payload['n1_exact_check'],
      'n2_histogram_check':payload['n2_histogram_check'],
      'finite_histogram_trace_worst_row':payload['finite_histogram_trace_worst_row'],
      'refinement':payload['refinement'],'fixed_representatives_q1_ell_one_third_shift0':reps,
      'worst_trace_row':max(rows,key=lambda r:r['trace_nontrivial']),
      'common_psd_scope':'Each D covers only the explicitly scanned finite family in the common canonical u-perp class basis; float matrices are not interval certificates. No common channel eigenbasis is assumed. Different tables may have different query marginals.',
      'common_psd_fixed_ell_worst_trace':compact_common(max(payload['common_psd_fixed_ell'],key=lambda r:r['trace_D'])),
      'common_psd_fixed_ell_minimum_margin':min(d['minimum_psd_margin'] for d in payload['common_psd_fixed_ell']),
      'common_psd_max_commutator_frobenius':max(d['max_commutator_frobenius'] for d in payload['common_psd_fixed_ell']),
      'common_psd_all_scanned_ell':{k:compact_common(v) for k,v in payload['common_psd_all_scanned_ell'].items()}}
    (OUT/'SUMMARY.json').write_text(json.dumps(short,indent=2))

if __name__=='__main__':main()
