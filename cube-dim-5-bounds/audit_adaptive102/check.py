#!/usr/bin/env python3
"""102 independent finite audit. Integer/Fraction only; never evaluates H/log."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,lcm,factorial,prod
from itertools import product,permutations,combinations
from collections import Counter
import hashlib,json,gzip,time,os,multiprocessing as mp,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent;T=Path('/Users/zhengzhihao/Desktop/T');B=T/'Math-Research-cube5-alt-cont/cube-dim-5-bounds';H=B/'upper_adaptive095';R=B/'upper_range092';A=B/'upper_adapt082';S=T/'cube-finite-dim-scratch/upper_adaptive095';N=1<<160;Q=1<<256;SHIFT=1<<96
TABLE=None;CAT=None;TS=None;MOM=None

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while z:=f.read(1<<20):h.update(z)
 return h.hexdigest()
def read(p):
 p=Path(p);return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def save(p,x):
 p=Path(p);t=Path(str(p)+'.tmp');t.write_bytes(gzip.compress(serial(x).encode(),mtime=0)) if p.suffix=='.gz' else t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def weight(ii):return factorial(4)//prod(factorial(v) for v in Counter(ii).values())
def key(L,ii):return str(L)+':'+','.join(map(str,ii))
def rawpath(k):return S/'evidence'/(k.replace(':','_').replace(',','-')+'.json.gz')
def offsets(j):return tuple((F(i,1<<j['level'])-F(x),F(i+1,1<<j['level'])-F(x)) for i,x in zip(j['indices'],j['point'][1:]))
def ospans(sp):return '|'.join(str(a)+','+str(b) for a,b in sp)
def sqrtpair(q,bits=224):
 D=1<<bits;n=isqrt(q.numerator*D*D//q.denominator);assert n*n*q.denominator<=q.numerator*D*D<(n+1)**2*q.denominator
 return F(n,D),F(n+1,D)
def manifest(d,expected):
 assert sha(d/'manifest.json')==expected
 fs=read(d/'manifest.json')['files']
 if isinstance(fs,dict):fs=[{'path':k,'sha256':v.get('sha256') if isinstance(v,dict) else v} for k,v in fs.items()]
 for r in fs:assert sha(d/r.get('path',r.get('name')))==r['sha256']
 return len(fs)
def fourth(c):
 b=c['bounds'];li,ui,mi=b['L_i'],b['U_i'],b['M_i'];assert 0<li<=mi<=ui<=3*N
 a=isqrt(isqrt(li**5));z=isqrt(isqrt(ui**5));assert a**4<=li**5<(a+1)**4 and z**4<=ui**5<(z+1)**4
 au=a+(a**4<li**5);zu=z+(z**4<ui**5);D=1<<(40+4*c['level'])
 assert b['integral_lower_i']*D<=a
 if li==ui:assert b['integral_upper_i']*D>=au
 else:assert b['integral_upper_i']*D*(ui-li)>=(ui-mi)*au+(mi-li)*zu

def moment(sp):
 # For strict interior anchors, Fj(r)=[min(r,-lo)+min(r,hi)]/width.
 assert all(a<0<b for a,b in sp);cuts=sorted({F(0)}|{abs(x) for ab in sp for x in ab});integ=F(0)
 for a,b in zip(cuts,cuts[1:]):
  mid=(a+b)/2;poly=[F(1)]
  for lo,hi in sp:
   thresholds=[-lo,hi];slope=sum(mid<t for t in thresholds);constant=sum(t for t in thresholds if t<mid);width=hi-lo
   out=[F(0)]*(len(poly)+1)
   for i,c in enumerate(poly):out[i]+=c*constant/width;out[i+1]+=c*slope/width
   poly=out
  integ+=sum(2*c*(b**(i+2)-a**(i+2))/F(i+2) for i,c in enumerate(poly))
 return cuts[-1]**2-integ

def budget_check(b,lo,up,reg):
 assert F(b['LH_truncated_integral'])==F(lo,N) and F(b['UH_truncated_integral'])==F(up,N)
 bv,iv,cv,ev,g=map(F,[b['B_upper'],b['I_upper'],b['cost_upper'],reg['E076_upper'],reg['valid_g']]);assert bv>ev and iv>0 and cv>0
 assert (bv-ev)**5>=F(up,N)**4 and iv**4>=bv**5 and cv**4>=(5*bv/g)**5
 assert F(b['source_E076_upper'])==ev and F(b['valid_g'])==g==F(1,10000)
 margin=32**4*g**5-5**5*bv**5;assert margin==F(b['exact_32_gate_margin']) and (margin>0)==b['strict32']
 assert b['no_tube_added'] and (iv<F(reg['old_complete_I_cap']))==b['I_strictly_below_old_cap']
 return {'LH':str(F(lo,N)),'UH':str(F(up,N)),'I_upper':str(iv),'cost_upper':str(cv),'strict32':margin>0}

def stage_ledger(reg,jobs,cert,pre):
 began=time.monotonic();nodes={}
 # Rebuild accepted old self ledger without reevaluating any old chord/H.
 for i in range(14):
  for c in read(R/f'block{i:04d}_cells.json.gz'):
   nodes[c['key']]={'key':c['key'],'indices':c['indices'],'level':c['level'],'weight':c['weight'],'selfL':c['new_self_L'],'selfU':c['new_self_U'],'children':None,'parent':None}
 roots=[key(5,r['indices']) for r in read(B/'upper_global080/roots.json.gz')];assert len(nodes)==56139 and len(roots)==52360
 for i in range(2):
  for par,cs in read(A/f'batch{i:04d}_registration.json')['children_by_parent'].items():
   nodes[par]['children']=cs
   for c in cs:assert nodes[c]['parent'] is None;nodes[c]['parent']=par
 def dfs(k):
  n=nodes[k]
  if n['children']:
   v=[dfs(c) for c in n['children']];lo=max(n['selfL'],sum(a for a,b in v));up=min(n['selfU'],sum(b for a,b in v))
  else:lo,up=n['selfL'],n['selfU']
  assert lo<=up;n['effL'],n['effU']=lo,up;return lo,up
 for k in roots:dfs(k)
 initial=read(reg['initial_state_path']);assert set(nodes)==set(initial['nodes']) and set(roots)==set(initial['roots'])
 for k,n in nodes.items():
  for f,v in n.items():assert initial['nodes'][k][f]==v,(k,f)
 front={k for k,n in nodes.items() if not n['children']};assert front==set(initial['frontier']) and len(front)==55819
 selected=sorted(front,key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:1024]
 assert selected==[r['key'] for r in pre['parents']]
 for r in pre['parents']:
  n=nodes[r['key']]
  for f in ['level','indices','weight','selfL','selfU']:assert r[f]==n[f]
  assert r['gap_units']==n['selfU']-n['selfL']
 seen=set();ordered_cells=0
 for pk in selected:
  par=nodes[pk];cs=pre['children_by_parent'][pk];expected=set()
  # Actual ordered child boxes, not only the weight identity.
  for ii in set(permutations(par['indices'])):
   for bits in product((0,1),repeat=4):expected.add(tuple(2*i+b for i,b in zip(ii,bits)))
  actual=set()
  for c in cs:
   assert c not in seen and c not in nodes;seen.add(c);j=jobs[c];v=cert[c];assert j['parent']==v['parent']==pk
   for f in ['level','indices','weight','key']:assert j[f]==v[f]
   ii=j['indices'];assert ii==sorted(ii,reverse=True) and [i//2 for i in ii]==par['indices'] and j['level']==par['level']+1
   assert c==key(j['level'],ii) and j['weight']==weight(ii)
   orbit=set(permutations(ii));assert actual.isdisjoint(orbit);actual.update(orbit);fourth(v)
  assert actual==expected and len(actual)==16*par['weight']==sum(jobs[c]['weight'] for c in cs);ordered_cells+=len(actual)
 assert seen==set(jobs)==set(cert) and len(seen)==13406
 rc=read(H/'run_receipt.json');oldL=sum(nodes[k]['effL'] for k in roots);oldU=sum(nodes[k]['effU'] for k in roots);budget_check(rc['before'],oldL,oldU,reg)
 volume=sum(F(nodes[k]['weight'],1<<(4*nodes[k]['level'])) for k in front);assert volume==1
 events=[json.loads(x) for x in gzip.decompress((H/'parent_events.jsonl.gz').read_bytes()).decode().splitlines()];assert len(events)==1024
 prev='0'*64;done=set();trace=[];unchanged={k:(nodes[k]['selfL'],nodes[k]['selfU']) for k in selected}
 for seq,event in enumerate(events,1):
  ev=dict(event);eh=ev.pop('event_hash');assert hashlib.sha256(serial(ev).encode()).hexdigest()==eh and ev['previous_hash']==prev and ev['seq']==seq;prev=eh
  pk=ev['parent'];assert pk in set(selected)-done and pk in front;par=nodes[pk];assert par['children'] is None;cs=pre['children_by_parent'][pk];assert ev['children']==cs
  for c in cs:
   j=jobs[c];b=cert[c]['bounds'];nodes[c]={'key':c,'indices':j['indices'],'level':j['level'],'weight':j['weight'],'selfL':j['weight']*b['integral_lower_i'],'selfU':j['weight']*b['integral_upper_i'],'children':None,'parent':pk}
  assert sum(nodes[c]['selfL'] for c in cs)==ev['child_raw_L_units'] and sum(nodes[c]['selfU'] for c in cs)==ev['child_raw_U_units']
  par['children']=cs;front.remove(pk);front.update(cs);done.add(pk)
  volume-=F(par['weight'],1<<(4*par['level']));volume+=sum(F(nodes[c]['weight'],1<<(4*nodes[c]['level'])) for c in cs);assert volume==1
  root=pk
  while nodes[root]['parent'] is not None:root=nodes[root]['parent']
  dfs(root) # recompute the entire affected root subtree, not the ancestor path
  lo=sum(nodes[k]['effL'] for k in roots);up=sum(nodes[k]['effU'] for k in roots)
  assert lo==ev['global_L_units'] and up==ev['global_U_units'] and lo-oldL==ev['root_delta_L_units']>=0 and up-oldU==ev['root_delta_U_units']<=0
  for k in set(selected)-done:assert nodes[k]['children'] is None and (nodes[k]['selfL'],nodes[k]['selfU'])==unchanged[k]
  trace.append({'seq':seq,'parent':pk,'global_L_units':lo,'global_U_units':up,'uncompleted_parents_paid_old':1024-seq,'volume':'1'});oldL,oldU=lo,up
 final=read(S/'state.json.gz');assert len(nodes)==69545 and len(front)==68201 and sorted(front)==read(H/'physical_frontier.json.gz')==final['frontier']
 assert set(final['nodes'])==set(nodes)
 for k,n in nodes.items():
  for f,v in n.items():assert final['nodes'][k][f]==v,(k,f)
 assert rc['event_hash']==prev and sha(S/'state.json.gz')==rc['state_sha256'] and rc['status']=='complete' and rc['workers_joined']
 bound=budget_check(rc['after'],oldL,oldU,reg);cap=F(54702517,10**12);assert F(bound['I_upper'])<=cap<F(reg['old_complete_I_cap']) and cap**4>=F(rc['after']['B_upper'])**5
 result={'status':'pass','initial_nodes':56139,'new_children':13406,'final_nodes':69545,'physical_frontier':68201,'ordered_new_child_boxes':ordered_cells,'all_parent_prefixes':1024,'volume':'1','event_hash':prev,'budget':bound,'I_simple':str(cap),'seconds':time.monotonic()-began,'trace':trace}
 save(P/'ledger_results.json',result);print('LEDGER_PASS',json.dumps({k:v for k,v in result.items() if k!='trace'}),flush=True)

def init(table,cat,ts,mom):
 global TABLE,CAT,TS,MOM
 TABLE=table;CAT={F(k):F(v['upper']) for k,v in cat.items()};TS=list(map(F,ts));MOM={k:F(v) for k,v in mom.items()}

def one(task):
 j,c=task;k=j['key']
 try:
  p=rawpath(k);assert sha(p)==c['evidence_sha256'];e=read(p)
  for f in ['key','level','indices','weight','parent','point','geometry']:assert e[f]==j[f]
  assert e['bounds']==c['bounds'];pt=list(map(F,j['point']));D=lcm(*(v.denominator for v in pt));assert D==e['D'] and pt[0]==1
  iv=[int(v*D) for v in pt];nn=sorted(sum(a*b for a,b in zip(iv,v)) for v in TABLE['vertices']);assert len(nn)==len(set(nn))==280 and 0 not in nn
  gap=F(min(b-a for a,b in zip(nn,nn[1:])),D);sp=offsets(j);assert all(a<0<b for a,b in sp);h=max(abs(v) for ab in sp for v in ab);mu=max(abs((a+b)/2) for a,b in sp);m2=MOM[ospans(sp)];g=j['geometry']
  assert g['offsets']==[[str(a),str(b)] for a,b in sp] and h==F(g['h_infinity']) and mu==F(g['mu_infinity']) and m2==F(g['M2']) and F(g['volume'])==F(1,1<<(4*j['level']))
  ends=e['endpoints'];assert [F(x['t']) for x in ends]==TS and len(ends)==73;pairs=[[0,0],[0,0],[0,0]];masses=[(F(x['H']['lower']),F(x['H']['upper'])) for x in ends];assert all(0<a<=b for a,b in masses);reuse=0
  for i,(a,(lo,up)) in enumerate(zip(TS,masses)):
   nl,nu=masses[i+1] if i<72 else (F(0),F(0));ml,mh=lo-nu,up-nl;assert 0<ml<=mh;rr=CAT[h/a];cf=min(rr,1+mu/a+F(25,48)*m2*rr/a**2)
   for pair,q in zip(pairs,[ml/rr,mh*rr,mh*cf]):
    num=q.numerator*Q;den=q.denominator;pair[0]+=num//den;pair[1]+=-(-num//den)
   if ends[i]['t'] in j['reuse_H']:assert a==F(1,64) and ends[i]['H']==j['reuse_H'][ends[i]['t']] and ends[i]['origin']=='exact_registered_cache';reuse+=1
   else:assert ends[i]['origin']=='genuine_new_H'
  bins=[(pairs[0][0]//SHIFT,pairs[0][1]//SHIFT)]+[(-(-a//SHIFT),-(-b//SHIFT)) for a,b in pairs[1:]];fallback=False
  if any(a!=b for a,b in bins):
   fallback=True;su=[F(0)]*3
   for i,(a,(lo,up)) in enumerate(zip(TS,masses)):
    nl,nu=masses[i+1] if i<72 else (F(0),F(0));rr=CAT[h/a];cf=min(rr,1+mu/a+F(25,48)*m2*rr/a**2)
    for z,q in enumerate([(lo-nu)/rr,(up-nl)*rr,(up-nl)*cf]):su[z]+=q
   vv=[(su[0]*N).__floor__(),(su[1]*N).__ceil__(),(su[2]*N).__ceil__()]
  else:vv=[a for a,b in bins]
  li,ui,mi=c['bounds']['L_i'],c['bounds']['U_i'],c['bounds']['M_i'];assert li<=vv[0] and min(3*N,vv[1])<=ui and min(ui,vv[2])<=mi and 0<li<=mi<=ui<=3*N
  assert e['runtime']['new_H_calls']==73-reuse and e['runtime']['reused_H_values']==reuse
  ww={int(a):F(b) for a,b in e['weights'].items()};tv=D**4*sum(map(abs,ww.values()))
  body=[(1,)+x for x in product((-1,1),repeat=4)];dots=[sum(a*b for a,b in zip(iv,v)) for v in body];norm=sum(v*v for v in iv);ds=[]
  for a,b in combinations(range(16),2):
   dot=sum(x*y for x,y in zip(body[a],body[b]));den=25-dot*dot;ds.append(F(norm*den-(5*dots[a]**2-2*dot*dots[a]*dots[b]+5*dots[b]**2),den*D*D))
  return {'key':k,'status':'pass','projection_gap':str(gap),'scaled_weight_variation':str(tv),'S2_distance_squared':str(min(ds)),'H512_lower':str(masses[0][0]),'H512_upper':str(masses[0][1]),'reused':reuse,'exact_bins':(li,ui,mi)==(vv[0],min(3*N,vv[1]),min(ui,vv[2])),'fallback':fallback,'cap_used':c['bounds']['global_cap_used']}
 except BaseException as exc:
  import traceback
  return {'key':k,'status':'failed','error':repr(exc),'traceback':traceback.format_exc()}

# Streaming just one top-level object prevents loading the 698 MiB cache as a Python tree.
def cache_entries(path):
 import mmap
 with Path(path).open('rb') as f:
  with mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as m:offset=m.find(b'"cache":{');assert offset>=0
  f.seek(offset+len(b'"cache":{'));buf='';pos=0;dec=json.JSONDecoder();eof=False
  while True:
   if pos>1<<19:buf=buf[pos:];pos=0
   while pos>=len(buf):
    x=f.read(1<<20);assert x;buf+=x.decode()
   while buf[pos] in ' \n\r\t,':pos+=1
   if buf[pos]=='}':return
   start=pos
   try:
    k,q=dec.raw_decode(buf,pos)
    while buf[q] in ' \n\r\t':q+=1
    assert buf[q]==':';q+=1
    while buf[q] in ' \n\r\t':q+=1
    v,end=dec.raw_decode(buf,q)
   except (json.JSONDecodeError,IndexError):
    buf=buf[start:];pos=0;x=f.read(1<<20);assert x,'truncated cache';buf+=x.decode();continue
   pos=end;yield k,v

def registry(jobs,reg):
 began=time.monotonic();G=T/'cube-finite-dim-scratch/upper_g105089';gr=read(G/'registration.json');ps=read(H/'preflight_summary.json');assert gr['epsilon']=='1/64' and gr['source_sha256']==ps['source089_base_sha256'];assert sha(gr['source'])==gr['source_sha256'];bm=read(T/'Math-Research-cube5-upper/cube-dim-5-bounds/upper_g105089/batch0008_manifest.json');overlayp=G/'batch0008/state.json';assert sha(overlayp)==bm['files'][str(overlayp)]==ps['source089_batch8_overlay_sha256'];assert bm['epsilon']=='1/64';over=read(overlayp)['new_cache'];assert len(over)==4744
 wanted={','.join(j['point']):j for j in jobs.values()};assert len(wanted)==13406;hits={};basecount=0
 for pk,v in cache_entries(gr['source']):
  basecount+=1
  if pk in wanted:hits[pk]=v
 assert basecount==180040
 for pk,v in over.items():
  if pk in wanted:hits[pk]=v
 saved=read(H/'H64_reuse_values.json.gz');assert set(hits)==set(saved) and len(hits)==32
 for pk,v in hits.items():
  assert v['point']==saved[pk]['point']==wanted[pk]['point'] and saved[pk]['t']=='1/64'
  assert saved[pk]['H']=={'lower':v['values'][0][0],'upper':v['values'][0][1]}==wanted[pk]['reuse_H']['1/64']
 for pk,j in wanted.items():assert bool(j['reuse_H'])==(pk in hits) and not j['reuse_H73']
 oldkeys={','.join(x['point']) for x in read(B/'upper_global080/roots.json.gz')}
 for bid in range(2):oldkeys.update(','.join(x['point']) for x in read(A/f'batch{bid:04d}_jobs.json.gz'))
 assert len(oldkeys)==56139 and not set(wanted)&oldkeys
 pre=read(H/'preflight_registration.json');assert sha(pre['pending9_manifest_path'])==pre['pending9_manifest_sha256'];pending={','.join(x['point']) for x in read(pre['pending9_manifest_path'])['points']};assert len(pending)==12 and not set(wanted)&pending
 overlap=read(B/'upper_adaptive103/overlap095_098.json');assert overlap['time']=='1/128' and overlap['source095_jobs_sha256']==sha(H/'jobs.json.gz')
 assert sha(overlap['source098_H128_point_manifest_path'])==overlap['source098_H128_point_manifest_sha256'];m098=read(overlap['source098_H128_point_manifest_path']);keys098={','.join(x['point']) for x in m098['records']};assert len(keys098)==194 and all(x['epsilon']=='1/128' for x in m098['records']) and not set(wanted)&keys098
 result={'status':'pass','H128_098_points':194,'H128_095_098_intersection':0,'base_count':basecount,'overlay_count':len(over),'H73_count':len(oldkeys),'new_point_count':len(wanted),'actual_H64_matches':32,'pending9_count':12,'pending_H73_intersections':0,'epsilon':'1/64','slot':0,'seconds':time.monotonic()-began};save(P/'registry_results.json',result);print('REGISTRY_PASS',json.dumps(result),flush=True)

def sample_weights(selected,jobs,table):
 results=read(P/'sample_checkpoint.json') if (P/'sample_checkpoint.json').exists() else []
 assert [r['key'] for r in results]==selected[:len(results)]
 for k in selected[len(results):]:
  e=read(rawpath(k));pt=list(map(F,e['point']));D=lcm(*(x.denominator for x in pt));nn=[sum(int(x*D)*v for x,v in zip(pt,r)) for r in table['vertices']];agg={}
  for cone in table['simplices']:
   ks=[nn[i] for i in cone['vertices']];poly=[1]
   for x in ks:
    nxt=[0]*(len(poly)+1)
    for i,c in enumerate(poly):nxt[i]-=x*c;nxt[i+1]+=c
    poly=nxt
   for x in ks:
    deriv=0
    for i in range(5,0,-1):deriv=deriv*x+i*poly[i]
    assert deriv!=0;value=F(-cone['det_numerator'],2*deriv);agg[x*x]=agg.get(x*x,F(0))+value
  agg={x:v for x,v in agg.items() if v};saved={int(x):F(v) for x,v in e['weights'].items()};assert agg==saved
  assert sum(agg.values())==0 and sum(x*v for x,v in agg.items())==0 and sum(x*x*v for x,v in agg.items())==F(-sum(c['det_numerator'] for c in table['simplices']),2)
  results.append({'key':k,'cones':3504,'folded_weights':len(agg),'all_weights_exact_match':True,'moments_0_2_4':True,'new_H_calls':0});save(P/'sample_checkpoint.json',results)
 return results

def main():
 started=time.monotonic();scope=read(P/'registration.json');assert sha(__file__)==scope['code_sha256'];assert not (P/'results.json').exists(),'completed finite audit refuses replay'
 for r in scope['sources']:assert sha(r['path'])==r['sha256'],r['path']
 mf={d:manifest(Path(d),h) for d,h in scope['manifests'].items()}
 reg=read(H/'registration.json');jobs={j['key']:j for j in read(H/'jobs.json.gz')};cert={c['key']:c for c in read(H/'children_certificates.json.gz')};pre=read(H/'preflight_registration.json');table=read(B/'upper_global080/geometry017.json.gz');assert len(table['vertices'])==280 and len(table['simplices'])==3504 and sha(B/'upper_global080/geometry017.json.gz')==reg['table_sha256'];assert len(jobs)==len(cert)==13406
 assert read(P.parent/'audit_range094/whole_results.json')['status']=='pass' and reg['source092_manifest_sha256']==sha(R/'manifest.json')
 accept=read(H/'source092_primary_acceptance.json');assert accept['I_upper']==reg['old_complete_I_cap']=='6010267/100000000000' and accept['qualified_g']==reg['valid_g']=='1/10000'
 for row in accept['frozen_inputs'].values():assert sha(row['path'])==row['sha256']
 Epath=T/'Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_shorttime076/results.json';assert reg['E076_upper']==next(c for c in read(Epath)['cases'] if c['epsilon']=='1/512')['E_norm']['upper'];assert reg['epsilon']=='1/512' and reg['q']=='5/4' and reg['time_endpoints']==read(R/'registration.json')['time_endpoints']
 cat=read(H/'ratio_catalog.json');oldcat=read(R/'ratio_catalog.json');assert len(cat)==97 and len(oldcat)==89 and all(cat[k]==v for k,v in oldcat.items());new=sorted(set(cat)-set(oldcat),key=F);assert len(new)==8 and new==read(H/'ratio_registration.json')['new_x'];checks=[]
 for s in ([] if (P/'ratio_results.json').exists() else new):
  x=F(s);a,b=sqrtpair(16+x*x);c,d=sqrtpair(256+x*x);lo=(a+x)**2*(c+x)**8/16**9;up=(b+x)**2*(d+x)**8/16**9;assert F(cat[s]['lower'])<=lo<up<=F(cat[s]['upper']);checks.append({'x':s,'lower224':str(lo),'upper224':str(up)})
 required={F(j['geometry']['h_infinity'])/F(a) for j in jobs.values() for a in reg['time_endpoints']};assert required==set(map(F,read(H/'ratio_registration.json')['used_x'])) and required<=set(map(F,cat));
 if not (P/'ratio_results.json').exists():save(P/'ratio_results.json',{'status':'pass','old89_reused':89,'new8':checks,'used_count':len(required)})
 else:assert read(P/'ratio_results.json')['status']=='pass' and len(read(P/'ratio_results.json')['new8'])==8
 # Consume exact old offset moments where available; calculate only genuinely new geometries.
 oldrows=read(R/'nodes.json.gz');mom={};moment_sources=[]
 for L in sorted({j['level'] for j in jobs.values()}):
  old=next((r for r in oldrows if r['level']==L),None)
  if old:
   assert sha(old['original_raw_path'])==old['original_raw_sha256'];e=read(old['original_raw_path']);sp=offsets({'point':e['point'],'indices':old['indices'],'level':L});mom[ospans(sp)]=old['M2'];moment_sources.append({'kind':'accepted094_same_offsets','level':L,'key':old['key'],'M2':old['M2']})
 newmom=read(P/'moment_results.json')['new_geometries'] if (P/'moment_results.json').exists() else []
 for row in newmom:mom[ospans(tuple((F(a),F(b)) for a,b in row['offsets']))]=row['M2']
 for j in jobs.values():
  sp=offsets(j);s=ospans(sp)
  if s not in mom:mom[s]=str(moment(sp));newmom.append({'offsets':[[str(a),str(b)] for a,b in sp],'M2':mom[s]})
  assert F(mom[s])==F(j['geometry']['M2'])
 save(P/'moment_results.json',{'status':'pass','old_reused':moment_sources,'new_geometries':newmom,'qualified_children':13406})
 if not (P/'ledger_results.json').exists():stage_ledger(reg,jobs,cert,pre)
 if not (P/'registry_results.json').exists():registry(jobs,reg)
 ordered=sorted(jobs);allout=[]
 for start in range(0,13406,2048):
  outp=P/f'raw_{start:05d}.json.gz'
  if outp.exists():out=read(outp)
  else:
   with mp.get_context('spawn').Pool(8,initializer=init,initargs=(table,cat,reg['time_endpoints'],mom)) as pool:out=list(pool.imap_unordered(one,[(jobs[k],cert[k]) for k in ordered[start:start+2048]],chunksize=8))
   failed=[r for r in out if r['status']!='pass']
   if failed:save(P/'failure.json',failed);raise AssertionError(failed[:1])
   save(outp,out)
  assert len(out)==min(2048,13406-start) and all(r['status']=='pass' for r in out);allout.extend(out);print('RAW_PASS',len(allout),flush=True)
 bykey={r['key']:r for r in allout};selectors={};selectors['first']=ordered[0];selectors['last']=ordered[-1]
 for field,reverse,name in [('projection_gap',False,'mingap'),('scaled_weight_variation',True,'maxvariation'),('S2_distance_squared',False,'nearS2')]:
  sel=sorted(allout,key=lambda r:((-F(r[field]) if reverse else F(r[field])),r['key']))[:2]
  for i,v in enumerate(sel):selectors[name+str(i)]=v['key']
 selectors['lowH']=min(allout,key=lambda r:(F(r['H512_lower']),r['key']))['key'];selectors['highH']=min(allout,key=lambda r:(-F(r['H512_upper']),r['key']))['key'];rk=sorted(r['key'] for r in allout if r['reused']);selectors['firstreuse']=rk[0];selectors['lastreuse']=rk[-1];selected=sorted(set(selectors.values()));assert len(selected)<=12
 save(P/'sample_selection.json',{'selectors':selectors,'selected':selected,'records':[bykey[k] for k in selected]});proof=sample_weights(selected,jobs,table);save(P/'sample_results.json',{'status':'pass','selection':selectors,'samples':proof,'method':'expanded integer cone polynomial derivative residues, no log/H'})
 assert sum(r['reused'] for r in allout)==32 and len(allout)*73-32==978606
 result={'status':'pass','pid':os.getpid(),'seconds':time.monotonic()-started,'source_registration_sha256':sha(P/'registration.json'),'code_sha256':sha(__file__),'manifest_counts':mf,'new_ratio_count':8,'old_ratios_consumed':89,'new_children':13406,'new_endpoint_inputs':978638,'new_H_evaluated_by_auditor':0,'shared_new_H_production_count':978606,'exact_H64_reused':32,'all280_guards':13406,'all_new_range_mean_bins_exact':all(r['exact_bins'] for r in allout),'rational_sum_fallbacks':sum(r['fallback'] for r in allout),'global_cap_uses':sum(r['cap_used'] for r in allout),'old_moment_geometry_reused':len(moment_sources),'new_moment_geometries':len(newmom),'PF_stencil_samples':len(proof),'PF_cones_checked':3504*len(proof),'full_ledger':{k:v for k,v in read(P/'ledger_results.json').items() if k!='trace'},'shared_H_Arb_production':True,'new_H_log_exp_Arb_calls':0,'new_C5_Cn_claim':False}
 save(P/'results.json',result);print('ALL_PASS',json.dumps(result),flush=True)
if __name__=='__main__':main()
