"""Armed four-point finite qualification; a failure stops without restart."""
from pathlib import Path
from fractions import Fraction as F
import json,time,os,sys
from flint import ctx
from zero_kernel import ZeroGreenKernel,sha,save,TABLE
P=Path(__file__).resolve().parent;CACHE=Path('/tmp/astra_green019')
POINTS=[(1,1,1,F(2,3),F(1,16)),(1,1,1,F(3,4),F(1,8)),(1,1,1,F(1,2),F(1,16)),(1,1,F(3,4),F(1,2),F(1,8))]
def main():
    sources={str(x):sha(x) for x in [Path(__file__),P/'zero_kernel.py',P/'registration.md',TABLE]}
    if '--arm' in sys.argv:
        print(json.dumps({'phase':'armed_before_math','pid':os.getpid(),'source_hashes':sources,'checkpoint_directory':str(CACHE)}),flush=True)
        input()
    began=time.monotonic();ctx.prec=256;engine=ZeroGreenKernel();records=[]
    for index,point in enumerate(POINTS):
        result=engine.evaluate(point,CACHE/('point_'+str(index)+'_checkpoint.json'),120)
        path=P/('point_'+str(index)+'.json');save(path,result)
        record={'index':index,'path':path.name,'sha256':sha(path),'status':result['status'],'point':result['point'],'seconds':result['seconds']}
        if 'G_ball' in result:record.update(G_ball=result['G_ball'],G_upper_dyadic=result['G_upper_dyadic'],zero_projection_vertices=result['zero_projection_vertices'],principal_parts=result['principal_coefficients_inside_common_negative_factor'])
        records.append(record)
        out={'status':'complete' if len(records)==4 else 'in_progress','records':records,'source_hashes':sources,
             'new_global_g':False,'new_C5_bound':False,'pid':os.getpid(),'seconds':time.monotonic()-began,'python':sys.executable}
        if not result['status'].startswith('PASS_'):out['status']='stopped_after_recorded_failure'
        save(P/'batch_result.json',out);print(json.dumps(record),flush=True)
        if not result['status'].startswith('PASS_'):return
    print(json.dumps({'status':'complete','pid':os.getpid(),'points':len(records),'seconds':time.monotonic()-began}),flush=True)
if __name__=='__main__':main()
