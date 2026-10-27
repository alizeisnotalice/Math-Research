import json,hashlib,ast,itertools,collections,sys
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round134/fixed_atlas_local_bisection');P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();canon=lambda x:(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode();key=lambda k,p:hashlib.sha256(canon({'semantic_key':k,'point':p})).hexdigest()
H=rd(A/'input_hashes.json');cfg=rd(A/'fixed_roots_and_atlas.json');entries=rd(A/'fixed_source_entries.json')['entries'];initial=rd(A/'initial_point_cache.json');cache=initial['entries'];req=rd(A/'release_requirements.json');packet=rd(A/'prepared_packet.json')
checks={};checks['all_frozen_source_hashes']=all(sha(p)==h for group in ['frozen_inputs','fixed_selected_source_artifacts']for p,h in H[group].items())
s132=A.parents[1]/'round132/cross_center_band_envelopes';s133=A.parents[1]/'round133/fixed_corner_envelope_diagnostic';oldcfg=rd(s132/'selected_boxes_and_centers.json');oldentries=rd(s132/'selected_source_cache_entries.json')['entries'];summary=rd(s133/'summary.json');chosen=summary['targets']['1/63500']['local_subdivision_candidate_node_ids'];roots=cfg['roots'];ids=[r['original_node_id']for r in roots];allkeys={k for r in roots for b in r['bands']for k in b['candidate_keys']};byid={x['node_id']:x for x in oldcfg['selected_boxes']}
checks['exact7_C133_candidates_in_sorted_order']=ids==sorted(chosen)and len(ids)==len(set(ids))==7
checks['per_original_atlas_unexpanded_exact_windows']=all(r['box']==byid[r['original_node_id']]['box']and[b['window']for b in r['bands']]==oldcfg['windows']and[b['candidate_keys']for b in r['bands']]==[[q['semantic_key']for q in b['candidate_keys']]for b in byid[r['original_node_id']]['bands']]and[len(b['candidate_keys'])for b in r['bands']]==[16,16,0,0,0]and all(len(b['candidate_keys'])==len(set(b['candidate_keys']))for b in r['bands'])for r in roots)
checks['82_union_cache_only_original_source_bytes']=len(entries)==len(allkeys)==82 and set(entries)==allkeys and all(entries[k]==oldentries[k]for k in allkeys)
checks['all_five_bands_positive_disjoint']=len(oldcfg['windows'])==5 and all(0<Q(b['a'])<Q(b['T'])for b in oldcfg['windows'])and all(Q(oldcfg['windows'][j+1]['T'])==Q(oldcfg['windows'][j]['a'])for j in range(4))
expected={};records=0;files=list(sorted((s132/'box_results').glob('*.json')))
for path in files:
 out=rd(path)
 for band in out['band_envelopes']:
  for center in band['center_bounds']:
   k=center['semantic_key']
   if k not in allkeys:continue
   for c in center['corners']:
    point=['1']+[out['box'][i][b]for i,b in enumerate(c['corner_bits'])];ck=key(k,point);fields={z:c[z]for z in ['displacement','gradient_dot_lower','clamped_true_value','frame32_quadratic','exponent']};value=(c['enclosure_lower'],c['enclosure_upper'],fields);records+=1
    if ck in expected:assert expected[ck]==value,'all duplicate original sourcepoint values agree'
    else:expected[ck]=value
checks['all93_C132_files_only82key_cache_complete']=len(files)==93 and records==47808 and len(expected)==16929 and set(expected)==set(cache)
checks['every_cache_identity_canonical_full6point']=all(k==key(r['semantic_key'],r['point'])and r['semantic_key']in allkeys and len(r['point'])==6 and r['point'][0]=='1'and all(str(Q(v))==v for v in r['point'])for k,r in cache.items())
checks['all16929_cached_values_and_formulas_unchanged']=all((r['lower'],r['upper'],r['exact_formula_fields'])==expected[k]and Q(0)<Q(r['lower'])<=Q(r['upper'])for k,r in cache.items())
checks['all_cache_original_provenance_hashes']=all(r['origin']=='C132_ORIGINAL_CORNER_NO_RECOMPUTATION'and str(s132/'box_results')in r['source_certificate_path']and sha(r['source_certificate_path'])==r['source_certificate_sha256']for r in cache.values())
# Reproduce only rational cache accounting; no original endpoint evaluator is called.
oldrecords=0;oldbounds_ok=True;allcorners_ok=True
for root in roots:
 points=[['1']+[root['box'][i][b]for i,b in enumerate(bits)]for bits in itertools.product((0,1),repeat=5)];boxlower=Q(0)
 for band in root['bands']:
  sourcebounds=[]
  for k in band['candidate_keys']:
   sourcebounds.append(min(Q(cache[key(k,p)]['lower'])for p in points));oldrecords+=len(points)
  boxlower+=max(sourcebounds,default=Q(0))
 source=rd(root['original_corner_certificate_path']);oldbounds_ok &= source['summed_disjoint_band_lower']==str(boxlower)and boxlower<Q(cfg['target'])
 for p in points:
  lo=sum(max((Q(cache[key(k,p)]['lower'])for k in b['candidate_keys']),default=Q(0))for b in root['bands']);allcorners_ok &= lo>=Q(cfg['target'])
checks['7_original_bounds_reused_no_exp']=oldbounds_ok and oldrecords==7168
checks['224_original_point_envelopes_paid']=allcorners_ok
# Pure rational geometry from the exact frozen function, not the scientific run/preflight.
tree=ast.parse((A/'run_batch.py').read_text());names={'require','points','volume','split_geometry','ordered'};ns={'Q':Q,'itertools':itertools};exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name in names],type_ignores=[]),'<frozen_geometry_only>','exec'),ns)
geometry_ok=True;nodes=leaves=splits=0
for root in roots:
 queue=collections.deque([('',root['box'])]);leafvol=Q(0);leafwords=[]
 while queue:
  word,box=queue.popleft();nodes+=1
  if len(word)==6:leaves+=1;leafvol+=ns['volume'](box);leafwords.append(word);continue
  g=ns['split_geometry'](box);splits+=1;axis=max(range(5),key=lambda j:(Q(box[j][1])-Q(box[j][0]),-j));mid=(Q(box[axis][0])+Q(box[axis][1]))/2
  lower=[b.copy()for b in box];upper=[b.copy()for b in box];lower[axis][1]=str(mid);upper[axis][0]=str(mid)
  geometry_ok &= g['axis_free_index']==axis and Q(g['midpoint'])==mid and g['lower_box']==lower and g['upper_box']==upper
  old={tuple(p)for p in ns['points'](box)};child={tuple(p)for childbox in [lower,upper]for p in ns['points'](childbox)};new={tuple(p)for p in g['midface_points']}
  geometry_ok &= len(new)==16 and child==old|new and old.isdisjoint(new)and all(all(str(Q(v))==v for v in p)for p in new)
  queue.append((word+'0',lower));queue.append((word+'1',upper))
 geometry_ok &= leafvol==ns['volume'](root['box'])and set(leafwords)=={''.join(x)for x in itertools.product('01',repeat=6)}
checks['closed_geometry_full_depth6']=geometry_ok and(nodes,leaves,splits)==(889,448,441)
checks['budget_889_448_441_225792']=cfg['max_depth']==6 and cfg['max_nodes_total']==889 and cfg['max_split_attempts_total']==441 and cfg['max_leaves_per_original_box']==64 and cfg['max_new_exp']==7*63*16*32==225792
paths={'run_batch_sha256':'run_batch.py','detached_runner_sha256':'detached_runner.py','detached_launcher_sha256':'launch_detached.py','fixed_roots_and_atlas_sha256':'fixed_roots_and_atlas.json','fixed_source_entries_sha256':'fixed_source_entries.json','initial_point_cache_sha256':'initial_point_cache.json','input_hashes_sha256':'input_hashes.json','proof_sha256':'PROOF.md','prepared_packet_sha256':'prepared_packet.json'}
checks['all_current_release_pins']=all(req[k]==sha(A/n)for k,n in paths.items())
checks['frozen_source_proof_requirements_match_requested']=sha(A/'run_batch.py')=='919dbcb0394862265b3afe49a8ef705aa438c525263d9d4c6f51063f92bd7f1d'and sha(A/'PROOF.md')=='ee6fa5fef6f58430f28fae5ea758ec92dd2b8936957f7e7b7d5a657833dbf88e'and sha(A/'release_requirements.json')=='b82763a7879fe4c8a14dbeee4bad15e63b5bb80620d796014d4a4effc8dc8cb7'
checks['AST_all_and_UNRUN']=all(ast.parse(p.read_text())is not None for p in A.glob('*.py'))and not(A/'started.json').exists()
result={'status':'PASS_C134_STATIC_FIXED_ATLAS_CACHE_AND_CLOSED_GEOMETRY_GATE'if all(checks.values())else'FAIL','checks':checks,'old_saved_cache_entries_checked':16929,'old_saved_records_checked':records,'original_source_union':82,'pure_geometry_fixture':{'nodes':nodes,'leaves':leaves,'splits':splits},'new_exp_native_old_baseline_evaluator_calls':0,'frozen_pins':{n:sha(A/n)for n in paths.values()},'scope':'seven full closed rectangles, own unchanged per-band16 centers, no global certificate or trueG upper'}
(P/'static_review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
