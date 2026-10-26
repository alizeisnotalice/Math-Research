"""Read-only all-ledger verification plus deterministic actual-band spot checks.
No H, R or geometry preflight rerun. Integer power checks certify saved roots.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,os,time,sys
import production_v2 as C
H=Path(__file__).resolve().parent;S=C.S

def main():
 t=time.monotonic();reg=C.read(H/'scalar_registration.json');old=C.read(reg['source103_state_path']);st=C.read(S/'state.json.gz');ns=st['nodes'];assert set(ns)==set(old['nodes']) and st['frontier']==old['frontier'] and st['roots']==old['roots'];rows={r['key']:r for r in C.read(H/'nodes.json.gz')};done=st['path112']['processed'];assert len(done)==len(set(done));certs={};blocks=[]
 for p in [S/'pilot_results.jsonl',*sorted(S.glob('block*_results.jsonl'))]:
  rr=[json.loads(z) for z in p.read_text().splitlines()];blocks.append({'path':str(p),'sha256':C.sha(p),'count':len(rr)})
  for r in rr:assert r['key'] not in certs;certs[r['key']]=r
 assert set(certs)==set(done)
 for k,n in ns.items():
  o=old['nodes'][k];r=rows[k];assert all(n[x]==o[x] for x in ['parent','children','weight','level','indices','selfL','effL']);assert (r['old_selfL'],r['old_selfU'])==(o['selfL'],o['selfU'])
  if k in certs:
   c=certs[k];assert c['packet_sha256']==r['packet_sha256'];assert r['old_L_i']<=c['new_M_i']<=min(r['old_M_i'],r['old_U_i']);assert n['selfU']==c['new_selfU']<=o['selfU'];assert c['old_selfU']==o['selfU']
  else:assert n['selfU']==o['selfU']
 # Independent bottom-up ledger, with identical fixed topology and lower side.
 L={};U={}
 for k in sorted(ns,key=lambda k:-ns[k]['level']):
  n=ns[k];cs=n['children'];L[k]=max(n['selfL'],sum(L[c] for c in cs)) if cs else n['selfL'];U[k]=min(n['selfU'],sum(U[c] for c in cs)) if cs else n['selfU'];assert L[k]==n['effL']<=n['effU']==U[k]
 ml=max(ns[k]['level'] for k in st['frontier']);assert sum(ns[k]['weight']*2**(4*(ml-ns[k]['level'])) for k in st['frontier'])==2**(4*ml)
 up=sum(U[k] for k in st['roots']);low=sum(L[k] for k in st['roots']);b=C.budget(up);E=F(b['E076']);B=F(b['B_upper']);I=F(b['I_upper']);g=F(b['g105']);T=F(b['S_lower'])+I;assert (B-E)**5>=F(up,C.N)**4 and I**4>=B**5 and T**4<=32**4*(g/5)**5;assert b['strict32']==(32**4*(g/5)**5>I**4)
 # Finite deterministic samples: first/last16 and32 smallest SHA ranks.
 selected=sorted(set(done[:16]+done[-16:]+sorted(done,key=lambda k:hashlib.sha256(k.encode()).hexdigest())[:32]));fac=C.read(H/'factor_catalog.json');samples=[]
 for k in selected:
  r=rows[k];e=C.read(r['packet_path']);assert C.sha(r['packet_path'])==r['packet_sha256'] and e['point']==r['point'];ends=e['endpoints'];assert [x['t'] for x in ends]==reg['time_endpoints'];total=F(0)
  for j,z in enumerate(ends):
   nxt=ends[j+1]['H'] if j+1<len(ends) else {'lower':'0','upper':'0'};lo=F(z['H']['lower'])-F(nxt['upper']);hi=F(z['H']['upper'])-F(nxt['lower']);assert lo>0 and hi>=lo;total+=hi*F(fac[r['geometry']][z['t']]['factor_upper'])
  mi=min(r['old_M_i'],r['old_U_i'],-(-total.numerator*C.N//total.denominator));assert mi==certs[k]['new_M_i'];l,u=F(r['old_L_i'],C.N),F(r['old_U_i'],C.N);theta=F(mi-r['old_L_i'],r['old_U_i']-r['old_L_i']) if u>l else F(0)
  # Independent high-precision inequality comparison to saved certificate;
  # binary224 rounding production must lie above exact384 lower chord.
  lp=C.root(l**5,4,384,False);upow=C.root(u**5,4,384,False);chord=((1-theta)*lp+theta*upow)/2**(4*r['level'])*r['weight'];assert F(certs[k]['new_selfU'],C.N)>=min(F(r['old_selfU'],C.N),chord)
  samples.append(k)
 out={'status':'pass','pid':os.getpid(),'seconds':time.monotonic()-t,'nodes':len(ns),'processed':len(done),'frontier':len(st['frontier']),'weighted_volume':'1','unchanged_range_L_U':True,'all_old_lower_and_topology_unchanged':True,'all_ancestor_min_self_childsum_pass':True,'all_processed_result_rows_checked':len(certs),'deterministic_actual_band_samples':samples,'sample_count':len(samples),'LH':str(F(low,C.N)),'budget':b,'results':blocks,'state_sha256':C.sha(S/'state.json.gz'),'verifier_sha256':C.sha(__file__),'new_H_R_jet_log_exp_asinh':0}
 C.write(H/'author_verification.json',out);print(C.ser({k:out[k] for k in ['status','nodes','processed','sample_count','seconds']}),flush=True)
if __name__=='__main__':main()
