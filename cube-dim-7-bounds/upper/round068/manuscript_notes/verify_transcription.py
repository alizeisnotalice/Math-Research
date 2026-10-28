"""New R068 paper-tuple transcription checks only; no old producer or baseline."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json,re,sys

rows=[]
def test(name,ok,left=None,right=None):
    assert ok,name
    r={'name':name,'status':'PASS'}
    if left is not None:r['left']=str(left)
    if right is not None:r['right']=str(right)
    rows.append(r)
def lt(name,a,b):test(name,a<b,a,b)
def le(name,a,b):test(name,a<=b,a,b)
def eq(name,a,b):test(name,a==b,a,b)

s=F(3,2);qlo=F(2,3);qhi=F(917,1000);delta=F(1,2**30)
rlo=1+delta;rhi=s*qhi-delta
sigma=F(83,1000);chi=F(9,200);hc=F(1,2**80)
e=F(1,2**20);epsA=F(1,65536);epsRet=F(1,4096)
P0=1-e-F(1,2**60);Ulo=P0*(1-e)/(1+e);Uhi=(1+e)/(P0*(1-e))
M=2*(1+e)*Uhi;b=2*(1-epsRet)*(1-e)*Ulo-2*e
m0=F(637,2000);Gamma=m0*b/M;gamma=F(317,1000)
Q5=F(34187,187500);eta=(1-qhi)/(2*qhi);Csrc=(qhi-qlo)/qlo
mode=sys.argv[1] if len(sys.argv)>1 else 'core'
if mode=='core':
    test('relative endpoint geometry',F(1,2)<qlo<qhi<1 and 1<rlo<rhi<F(3,2))
    lt('degree above sharp column',2,M);lt('positive row budget',0,b)
    test('reserved energy deficit',0<gamma<(1-e)*Gamma<Gamma<1)
    u=(rhi-rlo)/(rhi+rlo)
    mass=(1-epsA)*2*sum((u**(2*j+1)/F(2*j+1) for j in range(4)),F(0))
    lt('four-term rational target mass',m0,mass)
    S6=sum(((-1)**(j+1)*F(1,5)**j/F(j) for j in range(1,7)),F(0))
    eq('target mean lower-log sum',S6,F(2279,12500))
    eq('target mean direction',F(1,2)-F(5,2)*S6,F(221,5000))
    eq('odd-five pin upper',sum(((-1)**(j+1)*F(1,5)**j/F(j) for j in range(1,6)),F(0)),Q5)
    lt('source-face mean',F(2211,50000)/(1-e),F(223,5000))
    lt('source-face variance',F(2,3)*s*Q5**3/(1-e),F(1,150))
    eq('retention gap',chi-F(223,5000)-F(1,5000),F(1,5000))
    eq('whole-q log bridge',512*hc,F(1,2**71));lt('bridge reserved',512*hc,F(1,2**60))
    eq('full source mass coefficient',Csrc,F(751,2000))
    eq('cap exponent gap',eta-chi,F(47,183400))
    eq('tau exponent gap',sigma/2-chi*(1-sigma),F(47,200000))
    u=F(83,917);S7=sum(((-1)**(j+1)*u**j/F(j) for j in range(1,8)),F(0))
    lt('same-Phi exact price',S7,F(4333,50000))
    eq('exact-price display slack',gamma-F(4333,50000)-F(23,100),F(17,50000))
    # Abstract direct-sum overlap example; no original operator lower bound.
    eq('shared-source missing mass',1-F(1,4),F(3,4))
    lt('shared-source sum failure',F(3,4),2*(1-F(1,4)))
elif mode=='threshold':
    # Only the new finite-L30 threshold is checked in this run; core not replayed.
    L=F(2**30);rootL=F(2**15);x=(1-sigma)*L;h=hc/L;N=L**8/factorial(8)
    ML=F(factorial(8))*(F(3,4)/(qlo**2*L**6)+F(1,2)/(qlo*L**7))
    VL=F(3,2)*factorial(8)/(qlo**2*L**6)
    le('N lower',400,N);lt('t lower',16,N/L);lt('source radius',4*rhi,qlo*N/L)
    lt('trim width',h,delta/2);lt('local width',h,F(1,1000))
    lt('epsilon positive domain',1/(x*N**2),F(1,2))
    lt('relative gap drift',qhi/(x*N**2),2*h/3)
    lt('individual J precision',L/(3*qlo*N),e)
    lt('decreasing-pin count drift',2*h+1/N,3*e/4)
    lt('per-face cap fraction',F(1,2),1/rhi-L/N)
    lt('baseline mean',ML,e/2);lt('baseline Chebyshev',4*VL/e**2,epsA/8)
    lt('endpoint-count failure',4/(e**2*N),epsA/8)
    lt('target cap failure',F(3,4)**64,epsA/8)
    le('profile pin',200000,N*qlo)
    lt('profile Chebyshev',F(90000000000)/N,epsA/8)
    lt('profile bridge width',h,F(1,1600000))
    le('pin over x',1000,x)
    lt('perface retention',F(5000**2,150)/x,epsRet)
    lt('radial logwidth',1/(1-(h/rlo)**2),1+e)
    lt('finite n cap halfgap',eta/400,(eta-chi)/2)
    lt('finite n cap safe rate',F(1,10000),eta-chi-eta/400)
    lt('full source prefactor',Csrc*rhi,1)
    cS=F(917,10000000);nu=F(47,200000)
    lt('cap prefactor power',Csrc*rhi*L/hc,2**110)
    lt('cap exponential power',2**16,cS*L)
    lt('cap coarse Taylor10',F(factorial(10),2**50),e)
    lt('cap actual Taylor10',F(factorial(10))*Csrc*rhi/(hc*cS**10*L**9),e)
    lt('tau prefactor power',rhi*L*rootL/hc,2**126)
    lt('tau exponential power',2**17,nu*L)
    lt('tau coarse Taylor10',F(factorial(10),2**44),e)
    lt('tau actual Taylor10',F(factorial(10))*rhi/(hc*nu**10*L**8*rootL),e)
    lt('root lower domain',0,L/4-4/L-1)
    lt('exact-price exponential1',2**26,sigma*L)
    lt('exact-price exponential2',2**25,sigma*L/2)
    R4=2/L**2+F(factorial(4))/(sigma*L)**4+2*rootL*factorial(4)/(sigma*L/2)**4
    R4coarse=F(1,2**59)+F(factorial(4),2**104)+F(2**16*factorial(4),2**100)
    lt('full scalar power envelope',R4,R4coarse)
    lt('full scalar strict margin',R4coarse,F(17,50000))
    # Full functions compared via their positive B components and retained price.
    oldgamma56=F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208)
    lt('Phi56 full-cost buffer',oldgamma56-F(7,100)+F(4,2**20),F(23,100))
    lt('Phi64 full-cost buffer',F(61,200)-F(91,1000)+F(4,2**20),F(23,100))
else:
    raise SystemExit('mode must be core or threshold')
source=Path(__file__).resolve().parents[3]/'paper/main.tex'
tex=source.read_text();labels=re.findall(r'\\label\{([^}]+)\}',tex)
refs=re.findall(r'\\(?:ref|eqref)\{([^}]+)\}',tex)
test('unique labels',len(labels)==len(set(labels)));test('defined refs',set(refs)<=set(labels))
test('no drafting placeholders','@@' not in tex)
import hashlib
out={'actual_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'actual_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'new R068 fixed paper tuple/source checks, not a search or old baseline','mode':mode,'passed':len(rows),'failed':0,'labels':len(labels),'refs':len(refs),'checks':rows}
Path(__file__).with_name(mode+'_receipt.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(mode,len(rows),'checks PASS;',len(labels),'labels,',len(refs),'refs')
