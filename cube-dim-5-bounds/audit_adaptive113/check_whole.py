#!/usr/bin/env python3
"""113 finite independent103 audit; shared old095 qualification, zeroH/log."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,lcm,factorial,prod
from itertools import product,permutations,combinations
from collections import Counter
import hashlib,json,gzip,time,os,multiprocessing as mp,sys,copy
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent;T=Path('/Users/zhengzhihao/Desktop/T');B=T/'Math-Research-cube5-alt-cont/cube-dim-5-bounds';H=B/'upper_adaptive103';OLD=B/'upper_adaptive095';R=B/'upper_range092';A=B/'upper_adapt082';S=T/'cube-finite-dim-scratch/upper_adaptive103';N=1<<160;Q=1<<256;SHIFT=1<<96
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

def cache_entries(path,ts):
 import mmap
 with Path(path).open('rb') as f:
  with mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as m:start=m.find(b'"cache_by_epsilon":{');assert start>=0;marker=('\"'+ts+'\":{').encode();offset=m.find(marker,start);assert offset>=0
  f.seek(offset+len(marker));buf='';pos=0;dec=json.JSONDecoder();eof=False
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

def g105_check(b,up,event):
 g=F(21,200000);reg=read(H/'registration.json');E=F(reg['E076_upper']);Bv,Iv,Cv=map(F,[b['B_upper'],b['I_upper'],b['cost_upper']]);assert F(b['UH'])==F(up,N) and F(b['valid_g'])==g
 assert Bv>E and Iv>0 and Cv>0 and (Bv-E)**5>=F(up,N)**4 and Iv**4>=Bv**5 and Cv**4>=(5*Bv/g)**5
 margin=32**4*g**5-5**5*Bv**5;assert margin==F(b['exact_32_gate_margin']) and (margin>0)==b['strict32'];assert b['no_tube_added']
 assert b['parent_count']==event['seq'] and b['parent']==event['parent'] and b['event_hash']==event['event_hash']
 return {'I_upper':str(Iv),'cost_upper':str(Cv),'strict32':margin>0}

def ledger(reg,jobs,cert,pre):
 t=time.monotonic();initial=read(reg['initial_state_path']);assert sha(reg['initial_state_path'])==reg['initial_state_sha256'];oldproof=read(P.parent/'audit_adaptive102/results.json');assert oldproof['status']=='pass';nodes=copy.deepcopy(initial['nodes']);roots=initial['roots'];front=set(initial['frontier']);assert len(nodes)==69545 and len(front)==68201 and len(roots)==52360
 selected=sorted(front,key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:2048];assert selected==[r['key'] for r in pre['parents']]
 for r in pre['parents']:
  n=nodes[r['key']]
  for f in ['level','indices','weight','selfL','selfU']:assert r[f]==n[f]
  assert r['gap_units']==n['selfU']-n['selfL']
 seen=set();ordered_cells=0
 for pk in selected:
  par=nodes[pk];cs=pre['children_by_parent'][pk];expected=set()
  for ii in set(permutations(par['indices'])):
   for bits in product((0,1),repeat=4):expected.add(tuple(2*i+b for i,b in zip(ii,bits)))
  actual=set()
  for c in cs:
   assert c not in seen and c not in nodes;seen.add(c);j=jobs[c];v=cert[c];assert j['parent']==v['parent']==pk
   for f in ['level','indices','weight','key']:assert j[f]==v[f]
   ii=j['indices'];assert ii==sorted(ii,reverse=True) and [i//2 for i in ii]==par['indices'] and j['level']==par['level']+1
   assert c==key(j['level'],ii) and j['weight']==weight(ii);orb=set(permutations(ii));assert actual.isdisjoint(orb);actual.update(orb);fourth(v)
  assert actual==expected and len(actual)==16*par['weight']==sum(jobs[c]['weight'] for c in cs);ordered_cells+=len(actual)
 assert seen==set(jobs)==set(cert) and len(seen)==28427
 rc=read(H/'run_receipt.json');oldL=sum(nodes[k]['effL'] for k in roots);oldU=sum(nodes[k]['effU'] for k in roots);budget_check(rc['before'],oldL,oldU,reg)
 assert F(oldL,N)==F(oldproof['full_ledger']['budget']['LH']) and F(oldU,N)==F(oldproof['full_ledger']['budget']['UH'])
 volume=sum(F(nodes[k]['weight'],1<<(4*nodes[k]['level'])) for k in front);assert volume==1
 events=[json.loads(x) for x in gzip.decompress((H/'parent_events.jsonl.gz').read_bytes()).decode().splitlines()];gb=read(H/'g105_prefix_budgets.json.gz');assert len(events)==len(gb)==2048
 def dfs(k):
  n=nodes[k]
  if n['children']:
   z=[dfs(c) for c in n['children']];a=max(n['selfL'],sum(v for v,w in z));b=min(n['selfU'],sum(w for v,w in z))
  else:a,b=n['selfL'],n['selfU']
  assert a<=b;n['effL'],n['effU']=a,b;return a,b
 prev='0'*64;done=set();trace=[];unchanged={k:(nodes[k]['selfL'],nodes[k]['selfU']) for k in selected}
 for seq,event in enumerate(events,1):
  ev=dict(event);eh=ev.pop('event_hash');assert hashlib.sha256(serial(ev).encode()).hexdigest()==eh and ev['previous_hash']==prev and ev['seq']==seq;prev=eh
  pk=ev['parent'];assert pk in set(selected)-done and pk in front;par=nodes[pk];assert par['children'] is None;cs=pre['children_by_parent'][pk];assert ev['children']==cs
  for c in cs:
   j=jobs[c];b=cert[c]['bounds'];nodes[c]={'key':c,'indices':j['indices'],'level':j['level'],'weight':j['weight'],'selfL':j['weight']*b['integral_lower_i'],'selfU':j['weight']*b['integral_upper_i'],'children':None,'parent':pk,'source':'103'}
  assert sum(nodes[c]['selfL'] for c in cs)==ev['child_raw_L_units'] and sum(nodes[c]['selfU'] for c in cs)==ev['child_raw_U_units'];par['children']=cs;front.remove(pk);front.update(cs);done.add(pk)
  volume-=F(par['weight'],1<<(4*par['level']));volume+=sum(F(nodes[c]['weight'],1<<(4*nodes[c]['level'])) for c in cs);assert volume==1
  root=pk
  while nodes[root]['parent'] is not None:root=nodes[root]['parent']
  dfs(root);lo=sum(nodes[k]['effL'] for k in roots);up=sum(nodes[k]['effU'] for k in roots)
  assert lo==ev['global_L_units'] and up==ev['global_U_units'] and lo-oldL==ev['root_delta_L_units']>=0 and up-oldU==ev['root_delta_U_units']<=0
  for k in set(selected)-done:assert nodes[k]['children'] is None and (nodes[k]['selfL'],nodes[k]['selfU'])==unchanged[k]
  gproof=g105_check(gb[seq-1],up,event);trace.append({'seq':seq,'parent':pk,'global_L_units':lo,'global_U_units':up,'uncompleted_parents_paid_old':2048-seq,'volume':'1','g105_strict32':gproof['strict32']});oldL,oldU=lo,up
  if seq%256==0:save(P/'ledger_checkpoint.json',{'complete_prefixes':seq,'total':2048,'global_L_units':lo,'global_U_units':up})
 final=read(S/'state.json.gz');assert len(nodes)==97972 and len(front)==94580 and sorted(front)==read(H/'physical_frontier.json.gz')==final['frontier'];assert set(final['nodes'])==set(nodes)
 for k,n in nodes.items():
  for f,v in n.items():assert final['nodes'][k][f]==v,(k,f)
 assert rc['event_hash']==prev and sha(S/'state.json.gz')==rc['state_sha256'] and rc['status']=='complete' and rc['workers_joined']
 bound=budget_check(rc['after'],oldL,oldU,reg);gproof=g105_check(gb[-1],oldU,events[-1]);cap=F(51807888,10**12);assert F(bound['I_upper'])<=cap<F(reg['old_complete_I_cap']) and F(gproof['I_upper'])<=cap
 flo=read(P/'floor_results.json');TL,TU=F(flo['T_lower']),F(flo['T_upper']);IU=F(bound['I_upper']);IL=F(oldL,N);g=F(21,200000)
 assert 32**4*g**5-5**5*IU**4<0 and TL-IU<=TU-IU<0 and TL-IU<0<TU-IL
 result={'status':'pass','initial095_nodes_consumed':69545,'new_children':28427,'final_nodes':97972,'physical_frontier':94580,'ordered_new_child_boxes':ordered_cells,'all_parent_prefixes':2048,'volume':'1','event_hash':prev,'budget_g100':bound,'budget_g105':gproof,'I_simple':str(cap),'S_cert_lower':str(TL-IU),'S_cert_upper':str(TU-IU),'S_trueI_fixed_g105_lower':str(TL-IU),'S_trueI_fixed_g105_upper':str(TU-IL),'cap_improvement_over095':str(F(reg['old_complete_I_cap'])-cap),'seconds':time.monotonic()-t,'trace':trace}
 save(P/'ledger_results.json',result);print('LEDGER_PASS',json.dumps({k:v for k,v in result.items() if k!='trace'}),flush=True)

def registry(jobs,reg):
 t=time.monotonic();gp=T/'cube-finite-dim-scratch/upper_mixed098/final_state.json';assert sha(gp)==reg['source098_final_state_sha256'];wanted={','.join(j['point']):j for j in jobs.values()};assert len(wanted)==28427;hits={};counts={}
 for ts in ['1/64','1/128']:
  count=0;hit={}
  for pk,v in cache_entries(gp,ts):
   count+=1;assert ','.join(v['point'])==pk
   if pk in wanted:hit[pk]=v
  counts[ts]=count;hits[ts]=hit
 assert counts=={'1/64':184796,'1/128':194} and len(hits['1/64'])==230 and not hits['1/128']
 saved=read(H/'reuse_values.json.gz');assert set(saved)==set(hits['1/64'])
 for pk,v in hits['1/64'].items():
  h={'lower':v['values'][0][0],'upper':v['values'][0][1]};assert wanted[pk]['point']==v['point'] and wanted[pk]['reuse_H']=={'1/64':h};assert saved[pk]['1/64']['H']==h and saved[pk]['1/64']['t']=='1/64'
 for pk,j in wanted.items():assert bool(j['reuse_H'])==(pk in hits['1/64']) and not j['reuse_H73']
 old={','.join(x['point']) for x in read(B/'upper_global080/roots.json.gz')}
 for bid in range(2):old.update(','.join(x['point']) for x in read(A/f'batch{bid:04d}_jobs.json.gz'))
 old.update(','.join(x['point']) for x in read(OLD/'jobs.json.gz'));assert len(old)==69545 and not old&set(wanted)
 result={'status':'pass','H64_cache_points':184796,'H128_cache_points':194,'new_point_count':28427,'H64_actual_matches':230,'H128_actual_matches':0,'old_H73_points':69545,'old_H73_intersection':0,'actual_H64_same_point_time_slot0':True,'new_H_calls':0,'seconds':time.monotonic()-t};save(P/'registry_results.json',result);print('REGISTRY_PASS',json.dumps(result),flush=True)

def main():
 started=time.monotonic();scope=read(P/'whole_registration.json');assert sha(__file__)==scope['code_sha256'];assert not (P/'results.json').exists(),'completed audit refuses replay'
 for r in scope['sources']:assert sha(r['path'])==r['sha256'],r['path']
 mfs={d:manifest(Path(d),h) for d,h in scope['manifests'].items()};assert read(P/'floor_results.json')['status']=='pass' and read(P.parent/'audit_adaptive102/results.json')['status']=='pass'
 reg=read(H/'registration.json');jobs={j['key']:j for j in read(H/'jobs.json.gz')};cert={c['key']:c for c in read(H/'children_certificates.json.gz')};pre=read(H/'preflight_registration.json');table=read(B/'upper_global080/geometry017.json.gz');assert len(jobs)==len(cert)==28427 and len(table['vertices'])==280 and len(table['simplices'])==3504 and sha(B/'upper_global080/geometry017.json.gz')==reg['table_sha256']
 accept=read(H/'source095_primary_acceptance.json');assert accept['I_upper']==reg['old_complete_I_cap']=='54702517/1000000000000' and accept['qualified_g']=='21/200000';gp=H/'sourceg105_primary_acceptance.json';assert sha(gp)=='b249cc86042f4f5a8040cfdf21dcb7eca132813d0e0d531d2a40db1627045744' and read(gp)['qualified_g']=='21/200000'
 Epath=T/'Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_shorttime076/results.json';assert reg['E076_upper']==next(c for c in read(Epath)['cases'] if c['epsilon']=='1/512')['E_norm']['upper'];assert reg['epsilon']=='1/512' and reg['q']=='5/4' and reg['time_endpoints']==read(OLD/'registration.json')['time_endpoints']
 cat=read(H/'ratio_catalog.json');old=read(OLD/'ratio_catalog.json');assert len(cat)==170 and len(old)==97 and all(cat[k]==v for k,v in old.items());new=sorted(set(cat)-set(old),key=F);assert len(new)==73 and new==read(H/'ratio_registration.json')['new_x'];old091={r['x']:r for r in read(T/'Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_shift091/results.json')['ratios']};checks=[];overlap=[]
 if not (P/'ratio_results.json').exists():
  for s in new:
   if s in old091:
    assert [cat[s]['lower'],cat[s]['upper']]==old091[s]['integer_api_interval'];overlap.append(s);continue
   x=F(s);a,b=sqrtpair(16+x*x);c,d=sqrtpair(256+x*x);lo=(a+x)**2*(c+x)**8/16**9;up=(b+x)**2*(d+x)**8/16**9;assert F(cat[s]['lower'])<=lo<up<=F(cat[s]['upper']);checks.append({'x':s,'lower224':str(lo),'upper224':str(up)})
  save(P/'ratio_results.json',{'status':'pass','old97_reused':97,'table_new73':73,'old091_cases_reused':overlap,'actually_new_root_objects':len(checks),'checks':checks})
 required={F(j['geometry']['h_infinity'])/F(a) for j in jobs.values() for a in reg['time_endpoints']};assert required==set(map(F,read(H/'ratio_registration.json')['used_x'])) and required<=set(map(F,cat))
 oldjobs=read(OLD/'jobs.json.gz');mom={};sources=[]
 for j in oldjobs:
  sp=offsets(j);sig=ospans(sp)
  if sig not in mom:mom[sig]=j['geometry']['M2'];sources.append({'key':j['key'],'offsets':[[str(a),str(b)] for a,b in sp],'M2':mom[sig]})
 assert len(mom)==3;newmom=[];shapes=Counter()
 for j in jobs.values():
  sp=offsets(j);sig=ospans(sp)
  if sig not in mom:
   assert len(set(sp))==1 and all(a==-b for a,b in sp);h=sp[0][1];mom[sig]=str(F(2,3)*h*h);newmom.append({'offsets':[[str(a),str(b)] for a,b in sp],'M2':mom[sig],'method':'proved symmetric four-tail max moment 2h^2/3'})
  assert F(mom[sig])==F(j['geometry']['M2']);shapes[sig]+=1
 assert len(shapes)==4 and len(newmom)==1;save(P/'moment_results.json',{'status':'pass','old095_geometries_reused':sources,'new_geometries':newmom,'shape_counts':dict(shapes),'qualified_children':28427})
 if not (P/'ledger_results.json').exists():ledger(reg,jobs,cert,pre)
 if not (P/'registry_results.json').exists():registry(jobs,reg)
 ordered=sorted(jobs);allout=[]
 for start in range(0,28427,2048):
  outp=P/f'raw_{start:05d}.json.gz'
  if outp.exists():out=read(outp)
  else:
   with mp.get_context('spawn').Pool(8,initializer=init,initargs=(table,cat,reg['time_endpoints'],mom)) as pool:out=list(pool.imap_unordered(one,[(jobs[k],cert[k]) for k in ordered[start:start+2048]],chunksize=8))
   fail=[r for r in out if r['status']!='pass']
   if fail:save(P/'failure.json',fail);raise AssertionError(fail[:1])
   save(outp,out)
  assert len(out)==min(2048,28427-start) and all(r['status']=='pass' for r in out);allout.extend(out);print('RAW_PASS',len(allout),flush=True)
 bykey={r['key']:r for r in allout};selectors={'first':ordered[0],'last':ordered[-1]}
 for field,reverse,name in [('projection_gap',False,'mingap'),('scaled_weight_variation',True,'maxvariation'),('S2_distance_squared',False,'nearS2')]:
  sel=sorted(allout,key=lambda r:((-F(r[field]) if reverse else F(r[field])),r['key']))[:2]
  for i,v in enumerate(sel):selectors[name+str(i)]=v['key']
 selectors['lowH']=min(allout,key=lambda r:(F(r['H512_lower']),r['key']))['key'];selectors['highH']=min(allout,key=lambda r:(-F(r['H512_upper']),r['key']))['key'];rk=sorted(r['key'] for r in allout if r['reused']);selectors['firstreuse']=rk[0];selectors['lastreuse']=rk[-1];mid=sorted(k for k,j in jobs.items() if j['candidate']==0);assert len(mid)==2
 for i,k in enumerate(mid):selectors['midpoint'+str(i)]=k
 selected=sorted(set(selectors.values()));assert len(selected)<=14;save(P/'sample_selection.json',{'selectors':selectors,'selected':selected,'records':[bykey[k] for k in selected]});proof=sample_weights(selected,jobs,table);save(P/'sample_results.json',{'status':'pass','selection':selectors,'samples':proof,'method':'expanded integer polynomial derivative residues, zero H/log'})
 assert sum(r['reused'] for r in allout)==230 and len(allout)*73-230==2074941
 ratio=read(P/'ratio_results.json');led=read(P/'ledger_results.json');result={'status':'pass','pid':os.getpid(),'seconds':time.monotonic()-started,'registration_sha256':sha(P/'whole_registration.json'),'code_sha256':sha(__file__),'manifest_counts':mfs,'new_children':28427,'new_endpoint_inputs':2075171,'shared_new_H_production_count':2074941,'actual_H64_reused':230,'all280_guards':28427,'all_new_range_mean_bins_exact':all(r['exact_bins'] for r in allout),'rational_sum_fallbacks':sum(r['fallback'] for r in allout),'global_cap_uses':sum(r['cap_used'] for r in allout),'old97_ratios_consumed':97,'table_new73':73,'old091_ratio_cases_consumed':ratio['old091_cases_reused'],'actually_new_ratio_roots':ratio['actually_new_root_objects'],'old095_moment_geometries_consumed':3,'new_symmetric_moment_geometry':1,'PF_samples':len(proof),'PF_cones':3504*len(proof),'full_ledger':{k:v for k,v in led.items() if k!='trace'},'floor_first_stage_pass':True,'shared_H_Arb_production':True,'new_H_log_exp_Arb_calls':0,'new_C5_Cn_claim':False};save(P/'results.json',result);print('ALL_PASS',json.dumps(result),flush=True)
if __name__=='__main__':main()
