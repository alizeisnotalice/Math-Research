"""Eight finite exact arithmetic/symbolic checks; no maximal-set sampling."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,time,datetime
from flint import fmpq as Q
import sympy as S
H=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def checks():
    d=Q(1,10**6);a=Q(19,20);b=1-d;r0=a-d
    theta=Q(43,4000000);height=1-theta;mass=Q(2080)
    small=Q(1,2)+Q(7,20)/r0**6
    big=Q(1,2)/r0**6+Q(1,4)+Q(1,80)/r0
    volume=(b-a)*(2*b)**5;allowance=Q(43,62500)*mass
    rows=[]
    def add(name,ok,**data):assert ok,name;rows.append({'name':name,'PASS':True,**data})
    add('mass and finite L1 source',Q(1,2)*4**6+32==mass,mass=str(mass),spike_density=str(32/(2*d)**6))
    add('exact L6 allowance',theta==Q(43,15625)/256 and allowance==Q(17888,12500),theta=str(theta),allowance=str(allowance))
    add('all small components below h',small<height,upper=str(small),height=str(height),margin=str(height-small))
    add('unique big component below h',big<height,upper=str(big),height=str(height),margin=str(height-big))
    add('slab exceeds L6 allowance',volume>allowance,volume=str(volume),margin=str(volume-allowance),leakage_ratio_lower=str(volume/mass))
    r=S.symbols('r',positive=True)
    ds=S.factor(S.diff(S.Rational(1,2)+S.Rational(7,20)/r**6,r))
    add('small radial bound decreasing',ds==-S.Rational(21,10)/r**7,derivative=str(ds))
    db=S.factor(S.diff(S.Rational(1,2)/r**6+S.Rational(1,4)+S.Rational(1,80)/r,r))
    add('big radial bound decreasing',S.cancel(db+3/r**7+1/(80*r**2))==0,derivative=str(db))
    e=S.Rational(1,2)+S.Rational(1,2)/r**6-1
    factor=(1-r)*sum(r**i for i in range(6))/(2*r**6)
    add('strict full-input witness on 0<r<1',S.cancel(e-factor)==0,residual='0',positive_factorization=str(factor))
    assert len(rows)==8
    return {'status':'PASS_EIGHT_EXACT_CHECKS_ONLY','checks':rows,'kernel_calls':0,'union_or_grid_runs':0,'new_C6_bound':False}
def main():
    reg=json.loads((H/'registration.json').read_text())
    for path,h in reg['source_sha256'].items():assert sha(Path(path))==h
    if '--record' in sys.argv:
        assert not(H/'results.json').exists(),'one registered arithmetic run only'
        tic=time.perf_counter();result=checks();result.update(pid=os.getpid(),elapsed_seconds=time.perf_counter()-tic,python=sys.version,sympy=S.__version__,registration_sha256=sha(H/'registration.json'),exit_code=0)
        (H/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        result=checks();saved=json.loads((H/'results.json').read_text());assert result['checks']==saved['checks']
    print(result['status'],'0kernel/0union; independent full proof review separate')
if __name__=='__main__':main()
