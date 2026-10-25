from pathlib import Path
from fractions import Fraction as Q
import json,time,hashlib,os
P=Path(__file__).resolve().parent

def bd(x):
 e=Q(2)**int(x['common_binary_exponent']);return int(x['lower_mantissa'])*e,int(x['upper_mantissa'])*e

def main():
 assert not(P/'comparison_results.json').exists()
 clock=json.loads((P/'input_manifest.json').read_text());S=Path(clock['source_directory']);ref=json.loads((P/'reference_results.json').read_text());native=json.loads((S/'native_hessian.json').read_text());old=json.loads((S.parents[1]/'round016/band_engine/native_result.json').read_text())
 assert ref['logical_terms_each']==ref['computed_terms_each']==967680 and ref['zero_base_terms']==ref['repeated_base_terms']==0
 assert native['full_group_terms']==967680 and native['computed_sign_pairs']==483840 and native['confluent_terms']==native['zero_node_terms']==0
 assert len(native['hessian_upper_lex'])==21 and all(bd(v)[0]<=bd(v)[1] for v in native['hessian_upper_lex'])
 for name in ['band','gradient','p_a','p_T']:assert native[name]==old[name],name
 assert bd(native['Euler_xgrad_plus5B_plusTpT_minusapa'])[0]<=0<=bd(native['Euler_xgrad_plus5B_plusTpT_minusapa'])[1]
 checks=[]
 for name,index in [('x2x2',6),('x5x6',19)]:
  a,b=bd(native['hessian_upper_lex'][index]);lo,hi=bd(ref[name]);assert a<=lo<=hi<=b,(name,'containment')
  checks.append({'component':name,'upper_lex_index_zero_based':index,'reference_lower_exact':str(lo),'reference_upper_exact':str(hi),'native_lower_exact':str(a),'native_upper_exact':str(b),'native_contains_independent_reference':True,'reference_display_midpoint_only':format(float((lo+hi)/2),'.17g')})
 for f,h in clock['input_sha256'].items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h,f
 for f,h in clock['source_files_sha256'].items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h,f
 elapsed=time.time()-clock['started_unix'];assert elapsed<1200
 out={'status':'PASS_TWO_FULL_COMPONENTS_AND_SAVED_CONSISTENCY','checks':checks,'each_full_component_logical_terms':967680,'reference_computed_all_64_signs':True,'native_computed_32_signs_weight2':True,'old_nine_fields_byte_identical':True,'saved_21_fields_finite_ordered':True,'independent_full_component_count':2,'full_21_independent_reference':False,'no_box_remainder':True,'no_full_G_Hessian_replacement':True,'new_C6_bound':False,'pid':os.getpid(),'audit_elapsed_seconds':elapsed}
 (P/'comparison_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'components':[(v['component'],v['reference_display_midpoint_only'])for v in checks],'elapsed':elapsed}))
if __name__=='__main__':main()
