import json,pathlib,itertools,collections,time,datetime,hashlib
from fractions import Fraction as F
P=pathlib.Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes');D=pathlib.Path(__file__).parent;start=time.perf_counter(); inp=json.loads((P/'inputs.json').read_text());cones=json.loads((P/'cone_table.json').read_text())['cones'];assert (len(cones),inp['body_count'],inp['pair_count'])==(21,32,496); groups={};total=F(0);n=0
for cone in cones:
 a=F(cone['alpha']);nodes=[[F(x) for x in row] for row in cone['nodes_by_body']];assert len(nodes)==32 and all(len(row)==6 for row in nodes)
 for i,j in itertools.combinations(range(32),2):
  cols=tuple(sorted(zip(nodes[i],nodes[j])));swap=tuple(sorted((v,u) for u,v in cols));key=min(cols,swap);g=groups.setdefault(key,{'weight':F(0),'sources':[]});g['weight']+=a;g['sources'].append([cone['id'],i,j]);total+=a;n+=1
assert n==10416 and total==496*sum((F(c['alpha']) for c in cones),F(0));rows=[];ranks=collections.Counter()
for key,g in sorted(groups.items()):
 nz=[(u,v) for u,v in key if u or v];rank=0 if not nz else (2 if any(a*d-b*c for (a,b),(c,d) in itertools.combinations(nz,2)) else 1);ranks[rank]+=1;counts=collections.Counter(key);rows.append({'columns':[[str(a),str(b)] for a,b in key],'rank':rank,'zero_columns':6-len(nz),'repeated_column_excess':sum(v-1 for v in counts.values()),'weight':str(g['weight']),'source_count':len(g['sources']),'sources':g['sources']})
assert sum((F(x['weight']) for x in rows),F(0))==total
(D/'patterns.json').write_text(json.dumps({'indexing':'zero-based body indices as input row order','patterns':rows},separators=(',',':'))+'\n');receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start,'items':n,'pattern_count':len(rows),'rank_counts':dict(ranks),'total_alpha_weight':str(total),'expected_496_sum_alpha':str(496*sum((F(c['alpha']) for c in cones),F(0))),'checks':'item count, exact total weight and grouped exact total passed','calls':{'erfc':0,'kernel':0,'native':0,'exp':0},'patterns_sha256':hashlib.sha256((D/'patterns.json').read_bytes()).hexdigest()};(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
