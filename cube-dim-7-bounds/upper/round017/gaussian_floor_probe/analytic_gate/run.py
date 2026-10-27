#!/usr/bin/env python3
from flint import arb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,json,time,datetime,resource,sys,flint
ctx.prec=512;D=Path(__file__).resolve().parent;now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def A(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Fraction) else arb(x)
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def main():
 t=time.perf_counter();cpu=time.process_time();start=now();mf=json.loads((D/'r016-manifest.json').read_text());assert hashlib.sha256((D/'r016-result.json').read_bytes()).hexdigest()==mf['files']['result.json'];r=json.loads((D/'r016-result.json').read_text());J=Fraction(r['B_RB_rational_upper']);assert J==Fraction(4018421141119,10**18)
 pi=arb.pi();C=1/(2*pi**3*pi.sqrt());Dconst=2*arb(fmpq(65,2)).gamma()/arb(32).gamma();F=1/(105*(arb(2).log()*arb(fmpq(25,6))).exp()*(pi.log()*arb(fmpq(13,12))).exp());q=arb(fmpq(7,6));M0=((A(J)/F).log()*arb(fmpq(6,7))).exp();proxy=C*M0/343;required=7*(A(J/128).log()*arb(fmpq(6,7))).exp();barrier=(arb(7).log()*arb(fmpq(14,3))).exp()*pi**3/840;other=(arb(7).log()*q).exp()*F*((343/C).log()*q).exp();shape=C*Dconst/(343*arb(2048).sqrt());rhs=(A(Fraction(128))/barrier).log()*arb(fmpq(6,7));assert 0<proxy<required and barrier>128 and shape>0 and barrier.overlaps(other) and (proxy/required).overlaps(rhs.exp())
 enclosures={'K_barrier':ep(barrier),'U_proxy':ep(proxy),'g_required':ep(required),'U_shape':ep(shape)};out=dict(status='completed',J0=str(J),scalar_enclosures=4,precision_bits=512,enclosures=enclosures,proxy_strictly_below_required=True,barrier_strictly_above_128=True,barrier_forms_overlap=True,ratio_identity_overlap=True,scope='only g_cert<=E L_trace and F M_U^(7/6) numerator; not actual g upper or Cn lower')
 (D/'result.json').write_text(json.dumps(out,indent=2)+'\n');rc=dict(status='completed',start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,precision_bits=512,named_scalar_enclosures=4,new_shapes=0,old_cache_evaluations=0,python=sys.version,flint=flint.__version__);(D/'receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2));print(json.dumps({k:v['display'] for k,v in enclosures.items()},indent=2))
if __name__=='__main__':main()
