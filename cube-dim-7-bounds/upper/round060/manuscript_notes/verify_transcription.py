"""R060: exact checks of the newly transcribed R056 fixed theorem tuple.
No search, old baseline replay, or dimension-seven computation is performed.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json

records = []
def check(name, statement, left=None, right=None):
    assert statement, name
    entry = {"name": name, "status": "PASS"}
    if left is not None: entry["left"] = str(left)
    if right is not None: entry["right"] = str(right)
    records.append(entry)
def le(name, a, b): check(name, a <= b, a, b)
def lt(name, a, b): check(name, a < b, a, b)
def eq(name, a, b): check(name, a == b, a, b)

L = F(2**20); rootL = F(2**10); fact = F(factorial(8))
d = F(3,10); sigma = F(7,100); chi = F(7,200)
x = (1-sigma)*L; eps0 = F(1,2**24); h = eps0/(L*d)
e = F(1,1024); epsA = F(1,64)
qmin = F(31,40); qmax = F(37,40)
rlo = F(43,40)+F(1,2**20); rhi = F(49,40)-F(1,2**20)
N = L**8/fact
P0 = 1-F(1,256)-1280*eps0
Ulo = P0*(1-e)/(1+e); Uhi = (1+e)/(P0*(1-e)); M=F(21,8)
zA = (rhi-rlo)/rlo
mA = (1-epsA)*(zA-zA**2/2+zA**3/3-zA**4/4)
bstar = (1-e)*2*Ulo/qmax-2*e
Gamma = mA*bstar/M; gamma = (1-e)*Gamma
slack = gamma-F(7,93)-F(1,40)
Q = F(500,3083)-F(500,3083)**2/2+F(500,3083)**3/3
ML = 6*d**3*L**2/(qmin**2*N)+d*L/(qmin*N)
VL = 12*d**3*L**2/(qmin**2*N)
eq('product floor',P0,F(65275,65536))
eq('matching bridge',1280*L*d*h,F(5,65536))
eq('mass coefficient',mA,F(259465051078193863093554363,2018111068113594386771608802))
eq('gamma',gamma,F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208))
eq('scalar slack',slack,F(3593475707213855800847480839340601290687,764281719784595327319928601545123005726720))
lt('positive row excess',0,bstar); lt('row degree',2*Uhi/qmin,M)
lt('retained deficit',gamma,Gamma); lt('deficit below one',Gamma,1)
eq('pin polynomial',Q,F(13226208500,87910718361))
lt('pin bound',Q,F(1,6)); lt('mean rate',F(43,31)*Q**2/(1-e),F(4,125))
lt('variance rate',F(2,3)*F(43,31)*Q**3/(1-e),F(1,300))
eq('Chebyshev gap',chi-F(4,125)-F(1,2000),F(1,400))
eq('cap exponent',F(1,25)-chi,F(1,200))
eq('tau exponent',sigma/2-chi*(1-sigma),F(49,20000))
le('N lower bound',75,N); le('x lower bound',1000,x)
le('t lower bound',4,N/L); le('source radius',4*rhi,qmin*N/L)
for label,bound in [('trim',F(1,2**20)),('local width',d/8),('R floor',F(1,10000))]: le(label,h,bound)
le('shrinkage versus h',1/(x*N**2),h/2)
le('J uniform precision',2*d*L/(qmin*N),e)
le('baseline mean',ML,F(1,512)); le('baseline variance',VL,F(1,2**26))
le('endpoint count bad probability',4/(e**2*N),epsA/4)
le('count interpolation',3*h+1/N,e/2)
le('per-face target cap',(F(3,4))**74,epsA/4)
le('source-face Chebyshev',F(1600,3)/x,e)
le('pinned mean over x',Q/x,F(1,2000))
le('cap Taylor8',fact*d**2/(eps0*F(93,20000)**8*L**7),e)
le('tau Taylor8',fact*d/(eps0*F(49,20000)**8*L**6*rootL),e)
P = fact/(sigma*L)**8+2*rootL*fact/(chi*L)**8+fact/L**8
lt('complete root scalar',P,slack)
# The endpoint computations below check new expanded-domain constants, not J's old baseline.
eq('local derivative coarse',1+F(1,3)+3*F(1,3)**2,F(5,3))
eq('local numerator ratio',F(6,5)*F(1,3)**2+F(5,2)*F(1,3),F(29,30))
eq('J endpoint floor',1-F(5,32)*F(1,3),F(91,96))
for label,g in [('simple',F(9775559871,281956214374400)),('raw',F(735063024285513,21076775872652902400))]: lt('old '+label+' gamma',g,F(1,1000))
lt('full-cost comparison margin',F(1,1000)+4/L,F(1,40))
# R059 wide-band tuple as transcribed into the paper (no eventual N computed).
dw=F(391,1905); ab=F(283,2032); bb=F(63,64)
ew=F(1,2**20); hc=F(1,2**48); epsAw=F(1,65536); ret=F(1,4096)
qwmin=1-bb*dw; qwmax=1-ab*dw; R0=qwmin+dw
rwlo=R0+ew; rwhi=qwmax+dw-ew
Pw=1-ew-1280*hc; Ulw=Pw*(1-ew)/(1+ew); Uhw=(1+ew)/(Pw*(1-ew))
Mw=2*Uhw/qwmin; bw=(1-ret)*2*Ulw/qwmax-2*ew
massw=F(159477,1000000); Gw=massw*bw/Mw; gw=F(1637,12500)
check('wide local and band domains',0<dw<F(1,4) and F(1,2)<qwmin<qwmax<1 and 1<R0<rwlo<rwhi<F(5,4))
lt('wide inset width',2*ew,qwmax-qwmin)
lt('wide column versus degree',2,Mw); lt('wide row excess',0,bw)
check('wide deficit chain',0<gw<(1-ew)*Gw<Gw<1)
def slog(u): return 2*sum((u**(2*j+1)/F(2*j+1) for j in range(64)),F(0))
def elog(u): return 2*u**129/(129*(1-u*u))
uA=(rwhi-rwlo)/(rwhi+rwlo); uQ=dw/(4*R0-dw)
lt('wide target log mass',massw,(1-epsAw)*slog(uA))
Qw=F(26979,250000)
lt('wide pin logarithm',slog(uQ)+elog(uQ),Qw)
lt('wide pin level',Qw,F(1,9))
lt('wide inverse mean rate',R0*Qw**2/(qwmin*(1-ew)),F(293,20000))
lt('wide inverse variance rate',2*R0*Qw**3/(3*qwmin*(1-ew)),F(1,900))
etaw=ab*dw/(2*qwmax); chiw=F(147,10000); sigmaw=F(3,100)
eq('wide cap gap',etaw-chiw,F(499871,37603070000))
lt('wide cap positive exponent',F(1,100000),etaw-chiw)
eq('wide tau positive exponent',sigmaw/2-chiw*(1-sigmaw),F(741,1000000))
eq('wide scalar coefficient',gw-sigmaw/(1-sigmaw),F(121289,1212500))
eq('wide strict scalar slack',gw-sigmaw/(1-sigmaw)-F(1,10),F(39,1212500))
lt('wide versus explicit complete root',gamma-F(7,100)+4/L,F(1,10))

# Scan identifiers/cross references; this is a source check, not a theorem check.
import re
source = Path(__file__).resolve().parents[3]/'paper'/'main.tex'
s = source.read_text()
labels = re.findall(r'\\label\{([^}]+)\}',s)
refs = re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',s)
check('unique labels',len(labels)==len(set(labels)))
check('defined references',set(refs)<=set(labels))
check('no draft placeholders','@@' not in s)
output = {'kind':'exact fixed-tuple transcription checks; not new mathematical search',
          'passed':len(records),'failed':0,'checks':records,
          'labels':len(labels),'references':len(refs)}
outfile = Path(__file__).with_name('transcription_results.json')
outfile.write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
print(f'{len(records)} exact/source checks PASS; {len(labels)} labels, {len(refs)} references.')
