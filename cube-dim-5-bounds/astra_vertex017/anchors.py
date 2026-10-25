"""Exactly two new-versus-locked-old positive-time Arb qualification anchors."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,sys,time,os,json
from flint import arb,acb,fmpq,ctx
from exact import read,save,sha
P=Path(__file__).resolve().parent
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/numerics/body_fan_arb_anchor.py')
def q(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))
def new_heat(table,x,t):
    inverse=[]
    for z in table['vertices']:
        dot=sum((a*b for a,b in zip(x,z)),F(0))
        inverse.append(1/acb(q(144*t),-q(dot)))
    total=acb(0);absolute_sum=arb(0)
    for s in table['simplices']:
        term=acb(s['det_numerator'])
        for j in s['vertices']:term*=inverse[j]
        total+=term;absolute_sum+=abs(term)
    normalization=(2*acb.pi())**5;value=total/normalization
    return value,absolute_sum/(2*arb.pi())**5
def main():
    began=time.monotonic();ctx.prec=256
    table_path=P/'boundary_triangulation.json.gz';table=read(table_path)
    assert table['denominator']==144 and len(table['vertices'])==280 and len(table['simplices'])==3504
    sys.path.insert(0,str(OLD.parent))
    spec=importlib.util.spec_from_file_location('locked017oldanchor',OLD);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    inputs=[((F(1),F(7,8),F(3,4),F(5,8),F(1,2)),F(1,3)),((F(1),F(1),F(1),F(1),F(0)),F(1,8))]
    records=[]
    for x,t in inputs:
        start=time.monotonic();value,l1=new_heat(table,x,t);new_seconds=time.monotonic()-start
        start=time.monotonic();prior=old.heat(x,t);old_seconds=time.monotonic()-start
        difference=value.real-prior
        assert value.imag.contains(0) and difference.contains(0) and value.real.overlaps(prior) and value.real>0 and prior>0
        records.append({'point':[str(z) for z in x],'time':str(t),'new3504_real':str(value.real),'new3504_imag':str(value.imag),
                        'old26880_real':str(prior),'new_minus_old':str(difference),'new_term_absolute_sum':str(l1),
                        'new_cancellation_ratio':str(l1/value.real),'new_seconds':new_seconds,'old_seconds':old_seconds,
                        'checks':{'overlap':True,'difference_contains_zero':True,'imaginary_contains_zero':True,'both_strictly_positive':True}})
        save(Path('/tmp/astra_vertex017/anchors_checkpoint.json'),{'completed':len(records),'records':records})
        if time.monotonic()-began>120:raise RuntimeError('120-second complete-anchor checkpoint reached; no automatic repeat')
    out={'status':'two_registered_qualification_crosschecks_pass','precision_bits':ctx.prec,'records':records,
         'new_global_g':False,'new_C5_bound':False,'scope':'Geometry is proved separately; these two points only crosscheck implementation and normalization.',
         'source_hashes':{str(z):sha(z) for z in [Path(__file__),P/'exact.py',P/'registration.md',table_path,OLD,OLD.with_name('body_fan_geometry.py')]},
         'pid':os.getpid(),'seconds':time.monotonic()-began,'python':sys.executable}
    save(P/'anchors_result.json',out);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
