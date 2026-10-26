#!/usr/bin/env python3
"""Independent recursive DD, alternate envelopes, iterative hierarchy readback."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
from math import factorial,prod
import argparse,gzip,json,hashlib,multiprocessing as mp,os,time,resource,sys
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;OLD=H.parent/'upper_global080';S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adapt082');N=1<<160;TABLE=None;CERTS=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def key(L,i):return str(L)+':'+','.join(map(str,i))
def w(i):return factorial(4)//prod(factorial(v) for v in Counter(i).values())
def init(table,certs):
 global TABLE,CERTS
 TABLE=table;CERTS=certs;ctx.prec=384
def recursive(nodes):
 out={}
 for cone in TABLE['simplices']:
  ks=sorted(nodes[i] for i in cone['vertices']);tab=[{k*k:fmpq(1,2)} for k in ks]
  for order in range(1,5):
   nt=[]
   for i in range(5-order):nt.append({k:(tab[i+1].get(k,fmpq(0))-tab[i].get(k,fmpq(0)))/(ks[i+order]-ks[i]) for k in set(tab[i])|set(tab[i+1])})
   tab=nt
  for k,a in tab[0].items():out[k]=out.get(k,fmpq(0))-cone['det_numerator']*a
 return {k:a for k,a in out.items() if a}
def one(j):
 t0=time.monotonic();k=j['key'];path=S/'evidence'/(k.replace(':','_').replace(',','-')+'.json.gz');cert=CERTS[k];assert sha(path)==cert['evidence_sha256'];ev=json.loads(gzip.decompress(path.read_bytes()));assert ev['point']==j['point'] and ev['level']==j['level'] and ev['indices']==j['indices'];D=ev['D'];point=list(map(F,j['point']));nodes=[sum(int(x*D)*a for x,a in zip(point,z)) for z in TABLE['vertices']];assert len(set(nodes))==280;weights={int(k):rq(a) for k,a in ev['weights'].items()};assert recursive(nodes)==weights;assert sum(weights.values(),fmpq(0))==0 and sum((a*k for k,a in weights.items()),fmpq(0))==0;assert sum((a*k*k for k,a in weights.items()),fmpq(0))==-fmpq(sum(c['det_numerator'] for c in TABLE['simplices']),2)
 L=j['level'];lo=[F(i,2**L)-x for i,x in zip(j['indices'],point[1:])];hi=[F(i+1,2**L)-x for i,x in zip(j['indices'],point[1:])];mu=[(a+b)/2 for a,b in zip(lo,hi)];sec=[(a*a+a*b+b*b)/3 for a,b in zip(lo,hi)];Q=[[F(i==l,4)-F(1,36) for l in range(4)] for i in range(4)];eta=sum(Q[i][l]*(sec[i] if i==l else mu[i]*mu[l]) for i in range(4) for l in range(4));B2=max(sum(Q[i][l]*v[i]*v[l] for i in range(4) for l in range(4)) for v in product(*zip(lo,hi)));h=max(map(abs,lo+hi));lam=max(map(abs,mu));assert eta==F(ev['geometry']['eta']) and B2==F(ev['geometry']['B2']);ends=ev['endpoints'];assert len(ends)==73
 if j['cached_H1over64'] is not None:assert ends[24]['H']==j['cached_H1over64'] and ends[24]['origin']=='cached066_H1over64'
 else:assert all(x['origin']=='new082_scalar' for x in ends)
 low=arb(0);up=arb(0);mean=arb(0)
 for i,x in enumerate(ends):
  a=F(x['t']);nxt=ends[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};ml=F(x['H']['lower'])-F(nxt['upper']);mh=F(x['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh;z=arb(rq(1+B2/(10*a*a)));P=16*z**5-20*z**3+5*z;old=P+(P*P-1).sqrt();ex=arb(rq(h/(2*a))).exp()**2;R=min(ep(old,True),ep(ex,True));C=min(R,1+lam/a+F(21,8)*eta*R/a**2);low+=arb(rq(ml/R));up+=arb(rq(mh*R));mean+=arb(rq(mh*C))
 b=cert['bounds'];assert b==ev['bounds'];l=F(b['L_i'],N);u=F(b['U_i'],N);m=F(b['M_i'],N);assert l<=ep(low) and min(F(3),ep(up,True))<=u and min(u,ep(mean,True))<=m;assert 0<l<=m<=u<=3;v=arb(rq(F(1,2**(4*L))));q=arb(5)/4;lp=arb(rq(l))**q;upow=arb(rq(u))**q;chord=v*(lp*arb(rq((u-m)/(u-l)))+upow*arb(rq((m-l)/(u-l))));assert F(b['integral_lower_i'],N)<=ep(v*lp) and ep(chord,True)<=F(b['integral_upper_i'],N)
 return {'key':k,'pid':os.getpid(),'seconds':time.monotonic()-t0,'positive_bands':73}
def rebuild(last):
 old=json.loads(gzip.decompress((OLD/'initial_cells.json.gz').read_bytes()));roots=json.loads(gzip.decompress((OLD/'roots.json.gz').read_bytes()));nodes={}
 for xs in old['rows']:
  r=dict(zip(old['columns'],xs));i=roots[r['rank']]['indices'];k=key(5,i);nodes[k]={'level':5,'indices':i,'w':r['weight'],'lo':r['weight']*r['integral_lower_i'],'up':r['weight']*r['integral_upper_i'],'children':None}
 rootkeys=set(nodes);front=set(nodes)
 for bid in range(last+1):
  tag=f'{bid:04d}';br=json.loads((H/('batch'+tag+'_registration.json')).read_text());rc=json.loads((H/('batch'+tag+'_receipt.json')).read_text());certs=json.loads(gzip.decompress((H/('batch'+tag+'_cells.json.gz')).read_bytes()));lookup={c['key']:c for c in certs};expected=sorted(front,key=lambda k:(-(nodes[k]['up']-nodes[k]['lo']),nodes[k]['level'],nodes[k]['indices']))[:br['parent_count']];assert expected==[x['key'] for x in br['selected_parents']]
  for p,cs in br['children_by_parent'].items():
   assert p in front and nodes[p]['children'] is None;par=nodes[p];expectedcs={key(par['level']+1,sorted((2*a+b for a,b in zip(par['indices'],bits)),reverse=True)) for bits in product((0,1),repeat=4)};assert set(cs)==expectedcs;assert sum(lookup[c]['weight'] for c in cs)==16*par['w']
   for c in cs:
    assert c not in nodes;r=lookup[c];assert r['weight']==w(r['indices']);assert [i//2 for i in r['indices']]==par['indices'];b=r['bounds'];nodes[c]={'level':r['level'],'indices':r['indices'],'w':r['weight'],'lo':r['weight']*b['integral_lower_i'],'up':r['weight']*b['integral_upper_i'],'children':None}
   par['children']=cs;front.remove(p);front.update(cs)
  for k in sorted(nodes,key=lambda k:nodes[k]['level'],reverse=True):
   n=nodes[k];n['elo']=max(n['lo'],sum(nodes[c]['elo'] for c in n['children'])) if n['children'] else n['lo'];n['eup']=min(n['up'],sum(nodes[c]['eup'] for c in n['children'])) if n['children'] else n['up'];assert n['elo']<=n['eup']
  maxlevel=max(nodes[k]['level'] for k in front);assert sum(nodes[k]['w']*2**(4*(maxlevel-nodes[k]['level'])) for k in front)==2**(4*maxlevel);lo=F(sum(nodes[k]['elo'] for k in rootkeys),N);up=F(sum(nodes[k]['eup'] for k in rootkeys),N);assert lo==F(rc['after']['LH_truncated_integral']) and up==F(rc['after']['UH_truncated_integral']);assert sorted(front)==json.loads(gzip.decompress((H/('batch'+tag+'_frontier.json.gz')).read_bytes()));E=F(rc['after']['fixed_E076_upper']);B=ep(arb(rq(up))**(arb(4)/5)+arb(rq(E)),True);assert B<=F(rc['after']['B_upper']);g=F(1,10000);assert (5**5*F(rc['after']['B_upper'])**5<32**4*g**5)==rc['after']['conditional_cost_strictly_below32_exact']
 return {'LH':str(lo),'UH':str(up),'frontier_count':len(front),'weighted_volume':'1','node_count':len(nodes)}
def main():
 t0=time.monotonic();ctx.prec=384;ap=argparse.ArgumentParser();ap.add_argument('--batch',type=int,required=True);args=ap.parse_args();tag=f'{args.batch:04d}';reg=json.loads((H/'registration.json').read_text());assert sha(OLD/'manifest.json')==reg['source080_manifest_sha256'];jobs=json.loads(gzip.decompress((H/('batch'+tag+'_jobs.json.gz')).read_bytes()));certs={c['key']:c for c in json.loads(gzip.decompress((H/('batch'+tag+'_cells.json.gz')).read_bytes()))};table=json.loads(gzip.decompress((OLD/'geometry017.json.gz').read_bytes()));pool=mp.get_context('spawn').Pool(12,initializer=init,initargs=(table,certs))
 try:proofs=list(pool.imap_unordered(one,jobs));pool.close();pool.join()
 except BaseException:pool.terminate();pool.join();raise
 hierarchy=rebuild(args.batch);result={'status':'pass','batch':args.batch,'children_checked':len(proofs),'independent_recursive_DD_count':len(proofs),'positive_bands':sum(x['positive_bands'] for x in proofs),'alternate_Chebyshev_R_and_exp_checked':True,'independent_moments_chords_checked':True,'hierarchy':hierarchy,'all12_workers_joined':True,'worker_pids':sorted({x['pid'] for x in proofs}),'pid':os.getpid(),'seconds':time.monotonic()-t0,'source_sha256':sha(__file__),'new_H_evaluations':0};(H/('batch'+tag+'_verification.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
