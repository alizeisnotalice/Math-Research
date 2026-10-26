from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import lcm,factorial,prod
from collections import Counter
import gzip,hashlib,json,datetime
H=Path(__file__).resolve().parent
M=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
C=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-cover/cube-dim-5-bounds')
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-lower/cube-dim-5-bounds')
U=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert not (H/'registration.json').exists()
files={'table017':M/'astra_vertex017/boundary_triangulation.json.gz','code040':M/'upper_cutoff040/cutoff_box.py','report040':M/'upper_cutoff040/report.md','report029':M/'upper_alpha1_029/report.md','report067':C/'upper_integral067/report.md','stencil067':C/'upper_integral067/zero_stencil.json','report072':C/'upper_bandmean072/report.md','report074':C/'upper_bandlip074/report.md','report075':A/'audit_frame075/report.md','manifest075':A/'audit_frame075/manifest.json','manifest076':U/'upper_shorttime076/manifest.json','registration063':M/'upper_integral063/registration.json'}
src={k:{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for k,p in files.items()}
assert src['report075']['sha256']=='3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225'
assert src['report074']['sha256']=='492b8e23cd3ebd03be9874a16ec19a4d221448ef6a33fd5f396218c22c1d8baf'
t=json.loads(gzip.decompress(files['table017'].read_bytes()));verts=t['vertices']
nom=[('ordinary',[F(31,37),F(19,37),F(11,37),F(5,37)]),('low',[F(v,104729) for v in (104115,103354,71094,6327)]),('nearS1',[F(255,256),F(127,128),F(253,256),F(251,256)]),('crossing',[F(1),F(1),F(2,3),F(1,16)])]
offs=sorted(product((-1,0,1),repeat=4),key=lambda v:(sum(map(abs,v)),v))
jitters=[(0,0,0,0),(-1,-2,-4,-8),(-1,-3,-7,-11),(-2,-5,-11,-19),(-3,-7,-17,-29),(-2,-11,-23,-43),(-5,-13,-31,-47),(-7,-19,-37,-49)]
old=json.loads(files['registration063'].read_text());oldpoints={tuple(r['candidates_in_order'][r['selected_candidate']]['point']) for r in old['anchors']}
rows=[];seen=set()
for name,p in nom:
 base=tuple(min(63,int(64*x)) for x in p);picked=[]
 for off in offs:
  cell=tuple(sorted((max(0,min(63,b+d)) for b,d in zip(base,off)),reverse=True))
  if cell in seen:continue
  seen.add(cell);picked.append(cell)
  if len(picked)==8:break
 assert len(picked)==8
 for n,cell in enumerate(picked):
  center=[F(1)]+[F(2*i+1,128) for i in cell];candidates=[];selected=None
  for ji,jit in enumerate(jitters):
   c=[F(1)]+[x+F(v,3351328) for x,v in zip(center[1:],jit)];D=lcm(*(x.denominator for x in c));nodes=[sum(int(x*D)*z for x,z in zip(c,v)) for v in verts]
   generic=len(set(nodes))==280;inside=all(F(i,64)<x<F(i+1,64) for i,x in zip(cell,c[1:]));new=tuple(map(str,c)) not in oldpoints
   candidates.append({'jitter':jit,'point':list(map(str,c)),'D':D,'distinct':len(set(nodes)),'min_abs_projection':min(map(abs,nodes)),'inside':inside,'new_vs063':new})
   if selected is None and generic and inside and new:selected=ji
  assert selected is not None,(name,cell)
  rows.append({'id':f'{name}{n:02d}','class':name,'nominal':list(map(str,p)),'indices':cell,'weight':factorial(4)//prod(factorial(v) for v in Counter(cell).values()),'candidates':candidates,'selected':selected})
reg={'id':'upper_integral079','registered_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'model':'n5 alpha1 body16 unnormalized inverse Fourier (2pi)^-5','eps':'1/512','q':'5/4','level':6,'box_width':'1/64','precision_bits':384,'worker_count':12,'soft_checkpoint_seconds':540,'hard_alarm_seconds':590,'max_wall_seconds':600,'time_endpoints':[str(F(2**i,512)*(1+F(j,8))) for i in range(9) for j in range(8)]+['1'],'selection_rule':'for class in ordinary,low,nearS1,crossing: base floor(64 nominal) clamp 63; offsets {-1,0,1}^4 sorted by (L1 norm, tuple lex); coordinate clamp0..63, sort nonincreasing, take first8 globally unseen; anchors first generic of fixed8 jitter vectors about true box center','cells':rows,'sources':src,'mechanism_sha256':sha(H/'mechanism.md'),'registration_code_sha256':sha(Path(__file__)),'no_mesh_no_oldjet':True}
(H/'registration.json').write_text(json.dumps(reg,indent=2)+'\n')
print('REGISTERED',len(rows),sha(H/'registration.json'));print([(r['id'],r['indices'],r['selected']) for r in rows])
