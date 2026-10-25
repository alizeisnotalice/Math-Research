#!/usr/bin/env python3
"""Certify at most two exact-reconstructed candidates from registered float outputs."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,importlib.util,json,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
ENGINE=ROOT/'lower_opt013/event_candidates.py'
spec=importlib.util.spec_from_file_location('lower013_engine',ENGINE)
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();F,u,v,bounds,roots,den=engine.field()
 source=ROOT/'lower_round003/a2_equals_a3_results.json'
 C=engine.decode(F,json.loads(source.read_text())['coefficients'])
 screen=json.loads((HERE/'float_results.json').read_text());assert screen['status']=='completed_bounded_numerical_screen'
 def reconstruct(q):
  if 'rational'in q:return F(Q(q['rational']))
  def line(data):
   n,sign,slope=data
   return (F(0) if n==0 else sign*engine.root(n,u,v)/2),Q(slope,2)
  left,ls=line(q['line_left']);right,rs=line(q['line_right'])
  t=(right-left)/(ls-rs)
  assert left+t*ls==right+t*rs
  return t
 chosen=[];seen=set()
 for row in sorted(screen['rows'],key=lambda x:-x['ratio_display_only']):
  b,c=reconstruct(row['b']),reconstruct(row['c']);a=(b,b,c,c,c)
  if a in seen:continue
  seen.add(a);chosen.append((row,a))
  if len(chosen)==2:break
 assert len(chosen)<=2
 out={'status':'partial','C_source_sha256':sha(source),'C_coefficients':engine.coefficients(C),'source_sha256':sha(Path(__file__)),'engine_sha256':sha(ENGINE),
  'field_source_sha256':sha(ROOT/'lower/certify_rectangular.py'),'registration_sha256':sha(HERE/'registration.json'),'float_results_sha256':sha(HERE/'float_results.json'),
  'root_lower_numerators':roots,'root_denominator':str(den),'candidates':[]}
 for row,a in chosen:
  result={'selected_float_row':row}
  try:result.update(engine.calculate(F,u,v,bounds,C,a,(0,0,0,0,0),F(0)))
  except (RuntimeError,TimeoutError) as e:result.update(status='COMPUTE_LIMIT',reason=str(e))
  out['candidates'].append(result)
  (HERE/'exact_results.json').write_text(json.dumps(out,indent=2)+'\n')
  print(json.dumps({k:result[k] for k in ('status','display_only','difference_display_only','raw_boxes','retained_boxes','recursive_nodes','wall_seconds') if k in result}),flush=True)
 out.update(status='completed_selected_exact_certification',total_wall_seconds=time.monotonic()-start)
 (HERE/'exact_results.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
