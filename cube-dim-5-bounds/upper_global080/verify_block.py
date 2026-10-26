#!/usr/bin/env python3
"""Read a finite <=4096 block: exact identities, source reuse, independent chord."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod,factorial,ceil
from collections import Counter
import argparse,gzip,json,hashlib,time,os,resource,signal,sys
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_global080');N=1<<160

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def alarm(sig,frame):raise TimeoutError('590second verify alarm')
def main():
 start=time.monotonic();signal.signal(signal.SIGALRM,alarm);signal.alarm(590);ctx.prec=384;ap=argparse.ArgumentParser();ap.add_argument('--block',type=int,required=True);ap.add_argument('--full-bands',action='store_true');args=ap.parse_args();bid=f'{args.block:04d}';reg=json.loads((H/'registration.json').read_text());breg=json.loads((H/('block'+bid+'_registration.json')).read_text());receipt=json.loads((H/('block'+bid+'_receipt.json')).read_text());ledgerpath=S/'blocks'/bid/'ledger.json.gz';assert sha(ledgerpath)==receipt['ledger_sha256'];rows=json.loads(gzip.decompress(ledgerpath.read_bytes()));roots=json.loads(gzip.decompress((H/'roots.json.gz').read_bytes()));assert [r['rank'] for r in rows]==list(range(breg['start'],breg['end']));assert len(rows)<=4096
 sums=[0,0,0];capcount=0;maxratio=F(0);stages={};minmass=None
 for row in rows:
  rank=row['rank'];root=roots[rank];assert row['rootid']==root['rootid'] and row['weight']==root['weight'];assert row['weight']==factorial(4)//prod(factorial(m) for m in Counter(root['indices']).values());b=row['bounds'];assert b['positive_bands']==73;capcount+=b['global_cap_used'];maxratio=max(maxratio,F(b['relative_baseline_upper_i'],N));l,u,m=[F(b[k],N) for k in ('L_i','U_i','M_i')];assert 0<l<=m<=u<=3;q=arb(5)/4;v=arb(1)/1048576;lp=arb(rq(l))**q;up=arb(rq(u))**q;chord=v*(lp*arb(rq((u-m)/(u-l)))+up*arb(rq((m-l)/(u-l))))
  # Two equivalent Arb expressions can differ below2^-160; saved ceil has ample gap.
  assert F(b['integral_lower_i'],N)<=ep(v*lp);assert ep(chord,True)<=F(b['integral_upper_i'],N)
  sums[0]+=row['weight']*b['integral_lower_i'];sums[1]+=row['weight']*b['integral_upper_i'];sums[2]+=row['weight'];mm=F(b['minimum_band_mass_lower']);assert mm>0;minmass=mm if minmass is None else min(mm,minmass)
  if args.full_bands:
   path=S/'evidence'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz';assert sha(path)==row['evidence_sha256'];r=json.loads(gzip.decompress(path.read_bytes()));assert r['point']==root['point'] and r['indices']==root['indices'];ends=r['endpoints'];assert [x['t'] for x in ends]==reg['time_endpoints'];saved=ends[24];assert saved['t']=='1/64' and saved['origin']=='cached066_H1over64' and saved['H']==root['cached_H1over64'];assert sum(x['origin']=='cached066_H1over64' for x in ends)==1
   masses=[]
   for i,x in enumerate(ends):
    nxt=ends[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};ml=F(x['H']['lower'])-F(nxt['upper']);mu=F(x['H']['upper'])-F(nxt['lower']);assert 0<ml<=mu;masses.append(ml)
   assert min(masses)==mm
   # Independent direct second moments and all exact corner maxima.
   c=list(map(F,r['point']))[1:];lo=[F(i,32)-x for i,x in zip(root['indices'],c)];hi=[F(i+1,32)-x for i,x in zip(root['indices'],c)];means=[(l+h)/2 for l,h in zip(lo,hi)];second=[(l*l+l*h+h*h)/3 for l,h in zip(lo,hi)];Q=[[F(i==j,4)-F(1,36) for j in range(4)] for i in range(4)];eta=sum(Q[i][j]*(second[i] if i==j else means[i]*means[j]) for i in range(4) for j in range(4));B2=max(sum(Q[i][j]*w[i]*w[j] for i in range(4) for j in range(4)) for w in product(*zip(lo,hi)));assert eta==F(r['geometry']['eta']) and B2==F(r['geometry']['B2'])
  if time.monotonic()-start>=540:raise TimeoutError('finite readback540s soft stop')
 result={'status':'pass','block':args.block,'count':len(rows),'full_bands_checked':args.full_bands,'positive_band_count':len(rows)*73,'cached_endpoint_records_checked':len(rows) if args.full_bands else 0,'weighted_lower_units':sums[0],'weighted_upper_units':sums[1],'sum_weights':sums[2],'global_cap_used_count':capcount,'max_relative_baseline_display':float(maxratio),'minimum_band_mass_lower':str(minmass),'new_H_log_or_jet_evaluations':0,'seconds':time.monotonic()-start,'pid':os.getpid(),'source_sha256':sha(__file__),'ledger_sha256':sha(ledgerpath),'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(H/('block'+bid+'_verification.json')).write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
