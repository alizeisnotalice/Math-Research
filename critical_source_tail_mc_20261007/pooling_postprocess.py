from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,math
B=Path('/Users/zhengzhihao/Desktop/T/output/cube_general_20261003/geom_20261006/critical_source_tail_mc_20261007')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
regpath=B/'pooling_postprocess_registration.json'
outpath=B/'pooling_postprocess_results.json'
assert regpath.exists() and not outpath.exists()
reg=json.loads(regpath.read_text())
for name,h in reg['inputs_sha256'].items():assert sha(B/name)==h
original=json.loads((B/'results.json').read_text())
rounds=[]
for old in original['rounds']:
 n=old['n'];A=16*n;Ns=old['source_samples'];nr=old['receivers_per_replica']
 with gzip.open(B/f'n{n}_source_profiles.jsonl.gz','rt') as f:rows=[json.loads(line) for line in f]
 assert len(rows)==Ns==256 and nr==8
 selected=[];rest=[];allz=[]
 old_class_sum=F(0);rest_sum=F(0);new_first=F(0)
 for row in rows:
  y=list(map(F,row['source_y']))
  inclass=all(v.denominator==1 and abs(v)<=A-2 for v in y)
  scores=[[F.from_float(z['score']) for z in rep] for rep in row['receiver_replicas']]
  assert all(len(v)==nr for v in scores)
  s1,s2=[sum(rep,F(0))/nr for rep in scores]
  new_first+=(s1+s2)/2
  if inclass:
   selected.append(row['source_index']);old_class_sum+=s1*s2
   allz.extend(scores[0]+scores[1])
  else:
   rest.append(row['source_index']);rest_sum+=s1*s2
 Nc=len(selected);N=len(allz);assert N==2*nr*Nc
 if Nc:
  total=sum(allz,F(0));square=sum((z*z for z in allz),F(0))
  uc=(total*total-square)/(N*(N-1));assert uc>=0
  pool_class=F(Nc,Ns)*uc
 else:
  uc=F(0);pool_class=F(0)
 old_exact=(old_class_sum+rest_sum)/Ns
 new_exact=pool_class+rest_sum/Ns
 inp=json.loads((B/'registration.json').read_text())['source_inputs'][str(n)]['input']
 pclass=F((2*A-3)**n)*sum((F(w)/q**n for w,q in zip(inp['weights'],inp['qs'])),F(0))
 assert 0<pclass<1
 assert abs(float(old_exact)-old['I2_over_W']['mean'])<1e-12
 assert abs(float(new_first/Ns)-old['I1_over_W']['mean'])<1e-12
 # Exact arithmetic here certifies only the saved floating scores' postprocessing.
 raw=str(new_exact.numerator)+'/'+str(new_exact.denominator)
 rounds.append(dict(n=n,source_count=Ns,class_source_count=Nc,class_draw_count=N,
  class_source_indices=selected,class_probability_exact=str(pclass),class_probability=float(pclass),
  I1_original=old['I1_over_W']['mean'],I1_saved_score_exact_reconstruction=float(new_first/Ns),
  I2_original=old['I2_over_W']['mean'],I2_saved_score_exact_reconstruction=float(old_exact),
  I2_pooled_point=float(new_exact),class_old_contribution=float(old_class_sum/Ns),
  class_pooled_contribution=float(pool_class),nonclass_unchanged_contribution=float(rest_sum/Ns),
  pooled_class_U_point=float(uc),pooled_total_fraction_sha256=hashlib.sha256(raw.encode()).hexdigest(),
  pooled_total_fraction_bits=[new_exact.numerator.bit_length(),new_exact.denominator.bit_length()],
  original_moment_difference=old['I2_over_W']['mean']-old['I1_over_W']['mean']**2,
  pooled_moment_difference=float(new_exact)-old['I1_over_W']['mean']**2,
  pooled_I2_over_I1_point=float(new_exact)/old['I1_over_W']['mean'],
  no_new_confidence_interval=True,no_projection_to_moment_cone=True))
result=dict(status='complete',scope=reg['scope'],registration_sha256=sha(regpath),
 script_sha256=sha(Path(__file__)),rounds=rounds,new_samples=False,oracle_rerun=False,
 primary_files_modified=False,statistical_intervals_not_reused=True)
outpath.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps([dict(n=x['n'],Nc=x['class_source_count'],old_I2=x['I2_original'],pooled_I2=x['I2_pooled_point'],old_gap=x['original_moment_difference'],pooled_gap=x['pooled_moment_difference']) for x in rounds]))
