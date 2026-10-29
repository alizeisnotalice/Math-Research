"""Independent Decimal80 check of saved confidence formula and full finite profiles."""
from pathlib import Path
from decimal import Decimal as D, localcontext
import json,math,hashlib
H=Path(__file__).resolve().parent
S=H.parent/'round045_prior/SUMMARY.json'
def main():
 outpath=H/'CI_PROFILE.json';assert not outpath.exists()
 s=json.loads(S.read_text());N=s['config']['test_N'];hs=s['config']['log_thresholds'];records=[];maxerr=0.
 with localcontext() as c:
  c.prec=80;a=D('0.05')/D(336);ell=(D(2)/a).ln()
  assert s['directions']==336 and math.isclose(s['alpha_per_direction'],float(a),abs_tol=1e-20)
  for r in s['results']:
   for rr in r['rows']:
    for mu,var,ci in zip(rr['IS_weak'],rr['sample_variance'],rr['EB_CI']):
     rad=(D(2)*D(str(var))*ell/D(N)).sqrt()+D(7)*D(str(rr['integrand_cap']))*ell/(D(3)*D(N-1))
     got=[max(0.,mu-float(rad)),min(rr['integrand_cap'],mu+float(rad))]
     maxerr=max(maxerr,*[abs(x-y) for x,y in zip(ci,got)])
   mix=r['rows'][0];fixed=r['rows'][1:];assert all(mix['IS_weak'][j] <= sum(rr['IS_weak'][j] for rr in fixed)+1e-12 for j in range(len(hs)));j=max(range(len(hs)),key=lambda j:mix['IS_weak'][j]);b,k=max(((b,k)for b in range(4)for k in range(len(hs))),key=lambda bk:fixed[bk[0]]['IS_weak'][bk[1]])
   records.append(dict(n=r['n'],mixture_max=mix['IS_weak'][j],mixture_max_h=hs[j],bestfixed_max=fixed[b]['IS_weak'][k],bestfixed_b=b,bestfixed_h=hs[k],empirical_gap=mix['IS_weak'][j]-fixed[b]['IS_weak'][k],finite_R_pointwise_union_numerical_check=True,relative_to_bestfixed_profile=mix['IS_weak'][j]/fixed[b]['IS_weak'][k],simultaneous_profile_lower_separation=max(v[0]for v in mix['EB_CI'])-max(v[1]for rr in fixed for v in rr['EB_CI']),high_h32_mixture=dict(mean=mix['IS_weak'][-1],hits=mix['hits'][-1],EB_CI=mix['EB_CI'][-1])))
 assert maxerr<1e-12
 out=dict(status='PASS',models_thresholds_checked=165,directions=336,alpha_per_direction=str(a),Decimal_precision=80,max_EB_endpoint_difference=maxerr,profiles=records,summary_SHA256=hashlib.sha256(S.read_bytes()).hexdigest(),script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Audit of formula arithmetic and finite saved grid. EB interpretation assumes exact iid and exact weights/classifier; not an interval certificate or a continuum weak norm.')
 outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
