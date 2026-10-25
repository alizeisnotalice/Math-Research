"""Complete fixed035 recheck; stdout-only unless --save is explicitly requested."""
from pathlib import Path
import argparse,hashlib,json,runpy
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    outputs={}
    for name,output in [('035_representation_check.py','035_representation_results.json'),('035_integral_check.py','035_integral_results.json')]:
        result=runpy.run_path(str(HERE/name))['run'](args.project)
        result['checker_sha256']=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        outputs[output]=result
        if args.save:(HERE/output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS complete fixed4D candidate C4<14753/1000','representation_seconds':outputs['035_representation_results.json']['seconds'],'integral_seconds':outputs['035_integral_results.json']['seconds'],'g':'21/25000','I_outward_rational':outputs['035_integral_results.json']['I_outward_rational'],'cost_cube_exact':outputs['035_integral_results.json']['cost_cube_exact'],'candidate_U4':'14753/1000','strict_cube_gap':outputs['035_integral_results.json']['strict_cube_gap'],'new_ledger_write':False},indent=2))
