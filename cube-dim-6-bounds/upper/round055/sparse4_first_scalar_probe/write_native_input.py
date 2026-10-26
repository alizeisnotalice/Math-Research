"""New driver input construction; reuses frozen053 receipts, not its formal batch."""
import pathlib,json,gzip,hashlib,itertools,math,time
from fractions import Fraction as Q
from flint import fmpz_mat
P=pathlib.Path(__file__).resolve().parent;OLD=P.parents[1]/'round053/sparse4_orbit_scalar_engine'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if (P/'native_input.txt').exists():raise SystemExit('Refuse input regeneration')
 start=time.monotonic();reg=json.loads((P/'registration.json').read_text())
 for f,h in reg['source_sha256'].items():assert sha(pathlib.Path(f))==h
 s2=json.loads((P/'S2_certificate.json').read_text());assert s2['status']=='PASS_STRICT_OFF_TRUE_S2' and Q(s2['pair_min_squared'])>0
 # New input test2: independent 3x3 Gram determinant witness.
 x=reg['point_numerators'];v,w=s2['pair_min_witnesses'][0]['directions'];dot=lambda a,b:sum(q*r for q,r in zip(a,b))
 gram=[[dot(a,b) for b in [x,v,w]] for a in [x,v,w]];small=[[dot(a,b) for b in [v,w]] for a in [v,w]]
 witness=Q(int(fmpz_mat(gram).det()),2500*int(fmpz_mat(small).det()));assert witness==Q(s2['pair_min_squared'])
 src=next(pathlib.Path(f) for f in reg['source_sha256'] if f.endswith('global_fan.json.gz'));data=json.loads(gzip.decompress(src.read_bytes()))
 rays=[]
 for r in data['rays']:
  q=Q(r['q_true_152']);rays.append(dict(R=[q.denominator*z for z in r['primitive']],Q=q.numerator,scale=q.denominator))
 cones=[[c['abs_primitive_determinant']*math.prod(rays[i]['scale'] for i in c['ray_ids']),c['ray_ids']] for c in data['cones']]
 # Independent distinct multiset permutation generation, not46080 group replay.
 images=set()
 for a in set(itertools.permutations(x)):
  nz=[i for i,z in enumerate(a) if z]
  for signs in itertools.product([-1,1],repeat=5):
   z=list(a)
   for i,s in zip(nz,signs):z[i]*=s
   images.add(tuple(z))
 assert len(images)==3840
 reps=sorted({min(z,tuple(-v for v in z)) for z in images});assert len(reps)==1920
 packed=dict(rays=rays,cones=cones,point_orbit_representatives=reps,orbit_weight=24,point_denominator=50)
 digest=hashlib.sha256(json.dumps(packed,separators=(',',':'),sort_keys=True).encode()).hexdigest();old=json.loads((OLD/'adapter_results.json').read_text());assert digest==old['normalized_input_canonical_sha256']
 lines=['C055_SPARSE4_BETA1_80_X050_T32_825','478 13229 1920','50 24 32 825']
 lines += [' '.join(map(str,[r['Q'],*r['R']])) for r in rays]
 lines += [' '.join(map(str,[d,*ids])) for d,ids in cones]
 lines += [' '.join(map(str,x)) for x in reps]
 (P/'native_input.txt').write_text('\n'.join(lines)+'\n')
 out=dict(status='PASS_NEW_DRIVER_INPUT_AND_WITNESS',input_structure_tests=[dict(id=1,name='actual new152/11476S2 preflight',status=s2['status']),dict(id=2,name='independent exact Gram determinant witness',distance_squared=str(witness)),dict(id=3,name='new complete driver input and distinct multiset orbit construction',canonical_digest_matches053=True)],input_sha256=sha(P/'native_input.txt'),canonical_data_sha256=digest,input_bytes=(P/'native_input.txt').stat().st_size,rays=478,cones=13229,representatives=1920,weight=24,paired_terms=25399680,logical_terms_each_DD_field=609592320,joint_DD_field_evaluations=50799360,positive_time_evaluations=0,full_kernel_calls=0,seconds=time.monotonic()-start)
 (P/'input_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
