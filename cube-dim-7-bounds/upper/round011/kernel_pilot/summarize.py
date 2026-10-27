#!/usr/bin/env python3
import json
from pathlib import Path
p=Path(__file__).resolve().parent;d=json.loads((p/'results.json').read_text());keys=['top_requests','unique_top','unique_subdd','dd_cache_hits','log_evaluations','hermite_evaluations','zero_ray_knots','repeated_top_nodes'];out={}
for group in ['symmetry','generic_chunks']:
 phases=[y for x in d[group] for y in x['phase_records']];out[group]=dict(phases=len(phases),counts={k:sum(x['counts'][k] for x in phases) for k in keys})
out['scalar_counts']={x['name']:x['counts'] for x in d['scalar_cases']};out['origin_counts']={x['epsilon']:x['counts'] for x in d['origin']};out['real_scalar_main_top_requests']=len(d['scalar_cases']);out['origin_scalar_main_top_requests']=len(d['origin']);out['complex_scalar_main_and_checker_each']={'top_requests':4,'unique_subdd':sum(x['counts']['unique_subdd'] for x in d['scalar_cases']),'log_evaluations':sum(x['counts']['log_evaluations'] for x in d['scalar_cases']),'hermite_evaluations':sum(x['counts']['hermite_evaluations'] for x in d['scalar_cases']),'count_method':'reconstructed from identical cached recurrence over exact same node keys','near_expansion_terms':33};out['separate_checker_work']='four fixed scalar complex DD rechecks and 3 origin constants; no additional phases';(p/'counts.json').write_text(json.dumps(out,indent=2)+'\n')
