"""Independent recursive DD, rescaled roots, largest-coordinate moment and full graph."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
from math import factorial,prod,lcm,ceil,floor
from functools import lru_cache
import gzip,json,hashlib,importlib.util,multiprocessing as mp,os,time,sys
from flint import fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive095');N=1<<160;TABLE=None;CERT=None;ALT=None;TS=None

def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
V=module('verify092',H.parent/'upper_range092/verify.py');V.H=H
M=module('verify088',H.parent/'upper_reenvelope088/verify.py')
D=module('verify082',H.parent/'upper_adapt082/verify.py')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def key(L,ii):return str(L)+':'+','.join(map(str,ii))
def weight(ii):return factorial(4)//prod(factorial(v) for v in Counter(ii).values())
@lru_cache(maxsize=32)
def moment(spans):return M.conditional_maximum_moment(spans)
def init(table,cert,alt,ts):
 global TABLE,CERT,ALT,TS
 TABLE=table;D.TABLE=table;CERT=cert;ALT=alt;TS=ts

def one(j):
 k=j['key'];p=S/'evidence'/(k.replace(':','_').replace(',','-')+'.json.gz');assert sha(p)==CERT[k]['evidence_sha256'];e=json.loads(gzip.decompress(p.read_bytes()));assert e['point']==j['point'] and e['bounds']==CERT[k]['bounds'] and e['indices']==j['indices'] and e['parent']==j['parent'];point=list(map(F,j['point']));den=lcm(*(x.denominator for x in point));assert den==e['D'];nodes=[sum(int(x*den)*z for x,z in zip(point,v)) for v in TABLE['vertices']];assert len(set(nodes))==280;ww={int(x):fmpq(v) for x,v in e['weights'].items()};assert D.recursive(nodes)==ww;assert sum(ww.values(),fmpq(0))==sum((a*x for x,a in ww.items()),fmpq(0))==0;assert sum((a*x*x for x,a in ww.items()),fmpq(0))==-fmpq(sum(c['det_numerator'] for c in TABLE['simplices']),2)
 L=j['level'];spans=tuple((F(i,2**L)-x,F(i+1,2**L)-x) for i,x in zip(j['indices'],point[1:]));assert all(a<0<b for a,b in spans);h=max(abs(x) for span in spans for x in span);mu=max(abs((a+b)/2) for a,b in spans);m2=moment(spans);g=e['geometry'];assert g==j['geometry'] and h==F(g['h_infinity']) and mu==F(g['mu_infinity']) and m2==F(g['M2']);vol=F(1,2**(4*L));assert vol==F(g['volume']);ends=e['endpoints'];assert [x['t'] for x in ends]==TS;low=up=mean=0;Q=1<<256;reused=0
 for i,x in enumerate(ends):
  a=F(x['t']);assert F(x['H']['lower'])>0 and F(x['H']['lower'])<=F(x['H']['upper']);nxt=ends[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};ml=F(x['H']['lower'])-F(nxt['upper']);mh=F(x['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh;R=ALT[str(h/a)];c=min(R,1+mu/a+F(25,48)*m2*R/a**2);low+=floor(ml/R*Q);up+=ceil(mh*R*Q);mean+=ceil(mh*c*Q)
  if x['t'] in j['reuse_H']:assert x['H']==j['reuse_H'][x['t']] and x['origin']=='exact_registered_cache';reused+=1
  else:assert x['origin']=='genuine_new_H'
 b=e['bounds'];l,u,m=F(b['L_i'],N),F(b['U_i'],N),F(b['M_i'],N);assert l<=F(low,Q) and min(F(3),F(up,Q))<=u and min(u,F(mean,Q))<=m;assert 0<l<=m<=u<=3;il=F(b['integral_lower_i'],N);assert (il/vol)**4<=l**5;ll,lu=V.fourth_root_interval(l**5);ul,uu=V.fourth_root_interval(u**5);chord=vol*lu if l==u else vol*((u-m)*lu+(m-l)*uu)/(u-l);assert chord<=F(b['integral_upper_i'],N);assert e['runtime']['new_H_calls']==73-reused and e['runtime']['reused_H_values']==reused
 return {'key':k,'pid':os.getpid(),'reused':reused,'mathhash':hashlib.sha256(serial({'point':e['point'],'weights':e['weights'],'endpoints':ends,'bounds':b}).encode()).hexdigest()}

def hierarchy(reg):
 assert sha(reg['initial_state_path'])==reg['initial_state_sha256'];old=json.loads(gzip.decompress(Path(reg['initial_state_path']).read_bytes()));nodes=old['nodes'];front=set(old['frontier']);pre=json.loads((H/'preflight_registration.json').read_text());parents=pre['parents'];assert [x['key'] for x in parents]==sorted(front,key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:1024];events=[json.loads(x) for x in (S/'parent_events.jsonl').read_text().splitlines()];assert len(events)==1024;expectededges=pre['children_by_parent'];prev='0'*64;lo=sum(nodes[k]['effL'] for k in old['roots']);up=sum(nodes[k]['effU'] for k in old['roots']);initialL,initialU=lo,up
 for seq,ev in enumerate(events,1):
  hashval=ev.pop('event_hash');assert hashlib.sha256(serial(ev).encode()).hexdigest()==hashval and ev['previous_hash']==prev and ev['seq']==seq;prev=hashval;p=ev['parent'];cs=ev['children'];par=nodes[p];assert p in front and par['children'] is None;ii=par['indices'];assert set(cs)=={key(par['level']+1,sorted((2*a+b for a,b in zip(ii,bits)),reverse=True)) for bits in product((0,1),repeat=4)} and cs==expectededges[p];root=p
  while nodes[root]['parent'] is not None:root=nodes[root]['parent']
  oL,oU=nodes[root]['effL'],nodes[root]['effU']
  for c in cs:
   j=JOBS[c];b=CERT[c]['bounds'];assert c not in nodes and j['weight']==weight(j['indices']) and [x//2 for x in j['indices']]==ii;cl=j['weight']*b['integral_lower_i'];cu=j['weight']*b['integral_upper_i'];nodes[c]={'key':c,'level':j['level'],'indices':j['indices'],'weight':j['weight'],'selfL':cl,'selfU':cu,'effL':cl,'effU':cu,'children':None,'parent':p,'source':'095'}
  assert sum(nodes[c]['weight'] for c in cs)==16*par['weight'];assert sum(nodes[c]['selfL'] for c in cs)==ev['child_raw_L_units'] and sum(nodes[c]['selfU'] for c in cs)==ev['child_raw_U_units'];par['children']=cs;front.remove(p);front.update(cs);k=p
  while k is not None:
   n=nodes[k];n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in n['children']));n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in n['children']));assert n['effL']<=n['effU'];k=n['parent']
  dL,dU=nodes[root]['effL']-oL,nodes[root]['effU']-oU;assert dL==ev['root_delta_L_units']>=0 and dU==ev['root_delta_U_units']<=0;lo+=dL;up+=dU;assert lo==ev['global_L_units'] and up==ev['global_U_units'];ml=max(nodes[k]['level'] for k in front);assert sum(nodes[c]['weight']*2**(4*(ml-nodes[c]['level'])) for c in front)==2**(4*ml)
 actual=json.loads(gzip.decompress((S/'state.json.gz').read_bytes()));assert actual['nodes']==nodes and actual['frontier']==sorted(front);rc=json.loads((H/'run_receipt.json').read_text());assert rc['status']=='complete' and rc['workers_joined'] and rc['event_hash']==prev;assert sha(S/'state.json.gz')==rc['state_sha256'];b=rc['after'];lh,uh=F(lo,N),F(up,N);assert lh==F(b['LH_truncated_integral']) and uh==F(b['UH_truncated_integral']);B,I,C,E,g=map(F,[b['B_upper'],b['I_upper'],b['cost_upper'],reg['E076_upper'],reg['valid_g']]);assert B>E and (B-E)**5>=uh**4 and I**4>=B**5 and C**4>=(5*B/g)**5;assert (32**4*g**5-5**5*B**5)==F(b['exact_32_gate_margin']);assert (32**4*g**5>5**5*B**5)==b['strict32'];return {'LH':str(lh),'UH':str(uh),'physical_frontier':len(front),'nodes':len(nodes),'all1024_parent_prefixes_volume':'1','I_upper':str(I),'cost_upper':str(C),'strict32':b['strict32'],'initial_L_units':initialL,'initial_U_units':initialU,'event_hash':prev}

def main():
 global JOBS,CERT
 t=time.monotonic();reg=json.loads((H/'registration.json').read_text());assert sha(H/'run.py')==reg['run_sha256'] and sha(H/'kernel.py')==reg['kernel_sha256'] and sha(H/'jobs.json.gz')==reg['jobs_sha256'] and sha(H/'ratio_catalog.json')==reg['catalog_sha256'];jobs=json.loads(gzip.decompress((H/'jobs.json.gz').read_bytes()));JOBS={j['key']:j for j in jobs};CERT={x['key']:x for x in json.loads(gzip.decompress((H/'children_certificates.json.gz').read_bytes()))};alt=V.qualify_catalog();table=json.loads(gzip.decompress(Path(reg['table_path']).read_bytes()));assert sha(reg['table_path'])==reg['table_sha256'];pool=mp.get_context('spawn').Pool(12,initializer=init,initargs=(table,CERT,alt,reg['time_endpoints']));proof=[]
 try:
  for x in pool.imap_unordered(one,jobs):
   proof.append(x)
   if len(proof)%2048==0:print('VERIFIED',len(proof),round(time.monotonic()-t,3),flush=True)
  pool.close();pool.join()
 except BaseException:pool.terminate();pool.join();raise
 graph=hierarchy(reg);out={'status':'pass','children':len(proof),'new_H_evaluations':0,'independent_recursive_DD':len(proof),'all73_bands':73*len(proof),'reused_H64':sum(x['reused'] for x in proof),'rescaled_ratio_roots224':True,'integer_quarter_roots224':True,'conditional_maximum_moment':True,'all1024_prefixes_exact_hierarchy_and_volume':True,'hierarchy':graph,'pid':os.getpid(),'worker_pids':sorted({x['pid'] for x in proof}),'all_workers_joined':True,'seconds':time.monotonic()-t,'source_sha256':sha(__file__)};(H/'verification.json').write_text(json.dumps(out,indent=2)+'\n');(H/'mathematical_hashes.json.gz').write_bytes(gzip.compress(serial(sorted(proof,key=lambda x:x['key'])).encode(),mtime=0));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
