"""All actual geometry/time preflight and finite scalar registration;zeroR,H."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import gzip,json,hashlib,time
import pathmean as P
H=Path(__file__).resolve().parent;B=H.parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_path112');A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds');N=1<<160

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ser(x):return json.dumps(x,sort_keys=True,separators=(',',':'),default=str)
def read(p):return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def gz(p,x):p.write_bytes(gzip.compress(ser(x).encode(),mtime=0))
def key(L,ii):return str(L)+':'+','.join(map(str,ii))
def main():
 t=time.monotonic();assert not (H/'scalar_registration.json').exists();S.mkdir(exist_ok=True);assert sha(B/'upper_adaptive103/manifest.json')=='8abb73101ed1b328c0ee658a835dc22792379bbfa4297e8a842fda6276d36cb6';sp=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive103/state.json.gz');rc=read(B/'upper_adaptive103/run_receipt.json');assert sha(sp)==rc['state_sha256'];st=read(sp);nodes=st['nodes'];times=read(B/'upper_adaptive103/registration.json')['time_endpoints'];points={};packets={};bounds={};oldgeo={}
 for r in read(B/'upper_global080/roots.json.gz'):points[key(5,r['indices'])]=r['point']
 for p in sorted((B/'upper_adapt082').glob('batch*_jobs.json.gz')):
  for r in read(p):points[r['key']]=r['point']
 for r in read(B/'upper_range092/nodes.json.gz'):packets[r['key']]=(r['packet_path'],r['packet_sha256']);oldgeo[r['key']]=(F(r['M2']),F(r['h_infinity']),F(r['mu_infinity']))
 for p in sorted((B/'upper_range092').glob('block*_cells.json.gz')):
  for r in read(p):bounds[r['key']]={'L_i':r['rangeL_i'],'U_i':r['rangeU_i'],'M_i':r['mean_i']}
 for name in ('upper_adaptive095','upper_adaptive103'):
  jj=read(B/name/'jobs.json.gz');cc={r['key']:r for r in read(B/name/'children_certificates.json.gz')}
  for r in jj:
   k=r['key'];points[k]=r['point'];raw=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch')/name/'evidence'/(k.replace(':','_').replace(',','-')+'.json.gz');packets[k]=(str(raw),cc[k]['evidence_sha256']);bounds[k]=cc[k]['bounds'];g=r['geometry'];oldgeo[k]=(F(g['M2']),F(g['h_infinity']),F(g['mu_infinity']))
 assert set(points)==set(packets)==set(bounds)==set(nodes) and len(nodes)==97972;geomnodes={};offsets={};nodegeom={}
 for k,n in nodes.items():
  L=n['level'];point=list(map(F,points[k]));off=tuple((F(i,2**L)-x,F(i+1,2**L)-x) for i,x in zip(n['indices'],point[1:]));assert all(a<0<b for a,b in off);gkey=ser(off);offsets[gkey]=off;nodegeom[k]=gkey;geomnodes.setdefault(gkey,[]).append(k)
 gnames={s:f'g{i:03d}' for i,s in enumerate(sorted(offsets))};profiles={};plans=[];xs=set();rows=[]
 for gkey,off in offsets.items():
  gid=gnames[gkey];pr=P.weight_profile(off);assert pr['J0']==pr['M2']/2;profiles[gid]={'offsets':off,'profile':pr,'nodes':len(geomnodes[gkey])}
  for k in geomnodes[gkey]:assert (pr['M2'],pr['h'],pr['mu_infinity'])==oldgeo[k]
  for a0 in times:
   plan=P.darboux_plan(off,F(a0),8);ab=[]
   for c in plan['cells']:
    l,r,w=c['left'],c['right'],c['weight'];piece=next(p for p in pr['pieces'] if p['lo']<=l<r<=p['hi']);rho=P.integral([F(0),*piece['W']],l,r);alpha=(r*w-rho)/(r-l);beta=(rho-l*w)/(r-l);assert alpha>=0 and beta>=0 and alpha+beta==w;ab.append({'left':l,'right':r,'alpha':alpha,'beta':beta})
   assert sum((c['weight'] for c in plan['cells']),F(0))==pr['J0'];xs.update(plan['ratio_nodes']);plans.append({'geometry':gid,'time':a0,'plan_sha256':hashlib.sha256(ser(plan).encode()).hexdigest(),'alpha_beta_sha256':hashlib.sha256(ser(ab).encode()).hexdigest(),'panels':len(ab),'J0':str(pr['J0'])})
 for k in sorted(nodes):
  n=nodes[k];rows.append({'key':k,'level':n['level'],'indices':n['indices'],'weight':n['weight'],'geometry':gnames[nodegeom[k]],'point':points[k],'packet_path':packets[k][0],'packet_sha256':packets[k][1],'old_L_i':bounds[k]['L_i'],'old_U_i':bounds[k]['U_i'],'old_M_i':bounds[k]['M_i'],'old_selfL':n['selfL'],'old_selfU':n['selfU']})
 catalogs={};sources=[]
 def add(x,lo,up,origin):
  x=str(F(x));r={'lower':str(F(lo)),'upper':str(F(up)),'origin':origin};assert F(r['lower'])<=F(r['upper'])
  if x not in catalogs or F(r['upper'])-F(r['lower'])<F(catalogs[x]['upper'])-F(catalogs[x]['lower']):catalogs[x]=r
 paths=[B/'upper_adaptive103/ratio_catalog.json',A/'upper_pathmean111/ratio_catalog.json']
 for p in paths:
  sources.append({'path':str(p),'sha256':sha(p)})
  for x,r in read(p).items():add(x,r['lower'],r['upper'],{'path':str(p),'sha256':sha(p),'x':x})
 p=A/'upper_shift091/results.json';sources.append({'path':str(p),'sha256':sha(p)})
 for r in read(p)['ratios']:add(r['x'],r['R_new']['lower'],r['R_new']['upper'],{'path':str(p),'sha256':sha(p),'kind':'original091_R_new','x':r['x']})
 add('0','1','1',{'identity':'R(0)=1'});missing=sorted(x for x in xs if str(x) not in catalogs);gz(H/'nodes.json.gz',rows);gz(H/'geometries.json.gz',profiles);gz(H/'plans.json.gz',plans);(H/'inherited_ratio_catalog.json').write_text(ser(catalogs)+'\n');selected=sorted(st['frontier'],key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:256];reg={'status':'all97972_geometry73time8panels_registered_before_new_R','source103_manifest_sha256':sha(B/'upper_adaptive103/manifest.json'),'source103_state_path':str(sp),'source103_state_sha256':sha(sp),'source111_manifest_sha256':sha(A/'upper_pathmean111/manifest.json'),'pathmean_sha256':sha(H/'pathmean.py'),'secant_sha256':sha(H/'secant.py'),'prepare_sha256':sha(__file__),'nodes_sha256':sha(H/'nodes.json.gz'),'geometries_sha256':sha(H/'geometries.json.gz'),'plans_sha256':sha(H/'plans.json.gz'),'ratio_sources':sources,'inherited_catalog_sha256':sha(H/'inherited_ratio_catalog.json'),'all_used_x':list(map(str,sorted(xs))),'missing_x':list(map(str,missing)),'bits192_new_only':True,'new_R_registered':len(missing),'old_R_reused':len(xs)-len(missing),'panels_per_CDF_piece':8,'time_endpoints':times,'pilot_selected_physical_leaves':selected,'pilot_count':256,'nodes':97972,'physical_frontier':94580,'actual_geometries':len(profiles),'geometry_time_plans':len(plans),'positive_panels':sum(p['panels'] for p in plans),'new_H_calls':0,'new_R_calls_so_far':0,'seconds':time.monotonic()-t};(H/'scalar_registration.json').write_text(json.dumps(reg,indent=2)+'\n');print(json.dumps({k:reg[k] for k in ('nodes','physical_frontier','actual_geometries','geometry_time_plans','positive_panels','new_R_registered','old_R_reused','seconds')},indent=2))
if __name__=='__main__':main()
