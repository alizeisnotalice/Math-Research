"""One fixed positive T, one existing C010 H execution; refuses repeat."""
import json,subprocess,time,sys
sys.dont_write_bytecode=True
from pf_integral import P,ROOT,check_inputs,formal_checks,evaluate,sha,arb,ctx,ball,endpoints,bounds
def compare(pf,native):
 assert native['precision_bits']==512 and native['full_group_terms']==967680
 assert native['zero_node_terms']==88320 and native['confluent_terms']==201408
 assert native['computed_DD_evaluations']==483840 and native['exact_global_sign_pairing']
 ctx.prec=640
 wide=ball(pf['G_from_frozen_Gamma'])-ball(native['dyadic_interval']);wd=endpoints(wide)
 pl,pu=bounds(pf['K_T']);wl,wu=bounds(wd)
 passed=wl<pl<=pu<wu
 return dict(status='PASS' if passed else 'UNRESOLVED_NO_SEPARATION_STOP',PF640_strictly_inside_G_minus_H512=passed,independent_G_minus_H_reference=wd,p_T_crosscheck='Exact PF rational value; formal derivative recurrence, no independent full density call',native_full_group_terms=native['full_group_terms'],native_zero_node_terms=native['zero_node_terms'],native_confluent_terms=native['confluent_terms'],new_math_C6_bound=False)
def main():
 if (P/'execution_record.json').exists():raise SystemExit('Refuse any rerun')
 check_inputs();start=time.monotonic();record=dict(start_unix=time.time(),native_executions=0,positive_T_values=0,new_G_full=0)
 (P/'execution_record.json').write_text(json.dumps(record,indent=2)+'\n')
 s=formal_checks();(P/'scalar_checks.json').write_text(json.dumps(s,indent=2)+'\n')
 tick=time.monotonic();pf=evaluate();pf['PF_seconds']=time.monotonic()-tick
 (P/'pf_result.json').write_text(json.dumps(pf,indent=2)+'\n');record['positive_T_values']=1
 binary=ROOT/'Math-Research-cube6-audit/cube-dim-6-bounds/upper/round010/arb_engine/scratch/arb_kernel'
 fan=ROOT/'Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel/fan_integer_input.txt'
 cmd=[str(binary),'full',str(fan),'1','25','50','50','50','50','45','27','0']
 # Canonical literal tokens and SHA-pinned existing fan gate the legacy parser.
 record.update(native_executions=1,command=cmd,binary_sha256=sha(binary));(P/'execution_record.json').write_text(json.dumps(record,indent=2)+'\n')
 tick=time.monotonic();r=subprocess.run(cmd,capture_output=True,text=True,timeout=120);record.update(native_wall_seconds=time.monotonic()-tick,native_exit_code=r.returncode,native_stderr=r.stderr,elapsed_seconds=time.monotonic()-start)
 (P/'execution_record.json').write_text(json.dumps(record,indent=2)+'\n');(P/'native_H.json').write_text(r.stdout)
 if r.returncode:raise SystemExit('Native failure; stop without rerun')
 cross=compare(pf,json.loads(r.stdout));(P/'comparison.json').write_text(json.dumps(cross,indent=2)+'\n')
 print(json.dumps(dict(status=cross['status'],PF_seconds=pf['PF_seconds'],native_seconds=json.loads(r.stdout)['total_kernel_seconds'],elapsed_seconds=record['elapsed_seconds'])))
 if cross['status']!='PASS':raise SystemExit(2)
if __name__=='__main__':main()
