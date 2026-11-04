"""Exact incompatibility certificates and high-precision geometry regressions.
The universal asymptotic theorem is analytic, not certified by floating searches.
"""
from fractions import Fraction as F
from pathlib import Path
import sys,json,time,platform
from witness_capacity import atoms
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]/'results'

def box(x,y,l):return (x-l/2,x+l/2,y-l/2,y+l/2)
def geometry(b,r,lo,hi):
 a=1-b
 boxes=[box(a*r,b*r,lo),box(-a*r,-b*r,lo),box(b*r,a*r,hi),box(-b*r,-a*r,hi)]
 cap=atoms(boxes)
 IL=sum(v for m,v in cap.items() if m&3==3);IH=sum(v for m,v in cap.items() if m&12==12)
 bad=sum(v for m,v in cap.items() if m&12==12 and m&3==0)
 return IL,IH,bad

def run():
 st=time.monotonic();cert=[]
 for ts in ['1/1000','1/10000']:
  t=F(ts);alpha=1-t;b=F(9,20);a=1-b;r=t+(F(1,2)+2*a*b-F(1,100))*t*t;lo=alpha+t*t;hi=1/alpha-t*t
  IL,IH,bad=geometry(b,r,lo,hi);el=1+IL-2*alpha*lo;eh=1+IH-2*alpha*hi
  assert alpha<lo<hi<1/alpha
  assert min(lo*lo,1)>alpha*lo and min(hi*hi,1)>alpha*hi and el>0 and eh>0
  assert bad>el+eh
  cert.append(dict(t=ts,b='9/20',r=str(r),low_scale=str(lo),high_scale=str(hi),low_pair_slack=str(el),high_pair_slack=str(eh),incompatible_corner_area=str(bad),strict_contradiction_margin=str(bad-el-eh)))
 mp.mp.dps=80;rows=[]
 for bs in ['.05','.1','.2','.25','.35','.45']:
  b=mp.mpf(bs);a=1-b;d=a-b
  for ts in ['.01','.001','.0001']:
   t=mp.mpf(ts);al=1-t;w=1/al-al;r0=(1-al*al)/(al+mp.sqrt(al*al-4*a*b*(1-al*al)))
   IL,IH,bad=geometry(F(bs),F(mp.nstr(r0,75)),F(1)-F(ts),1/(F(1)-F(ts)))
   formula=2*((w/2)**2-d*d*r0*r0);assert abs(mp.mpf(bad.numerator)/bad.denominator-formula)<mp.mpf('1e-70')
   rows.append(dict(b=bs,t=ts,corner_area_over_t2=mp.nstr(formula/t**2,35),limit=mp.nstr(8*a*b,35)))
 # Convex support chart and triangle inclusion regression on T itself.
 b=mp.mpf('0.45');a=1-b;eta=mp.mpf('1e-4');s0=mp.sqrt(eta)/16;charts=[]
 for ts in ['.01','.001','.0001']:
  t=mp.mpf(ts);al=1-t;c=2*al*al-1;r0=(1-al*al)/(al+mp.sqrt(al*al-4*a*b*(1-al*al)));p=mp.matrix([a*r0,b*r0]);nn=mp.matrix([al-2*p[1],al-2*p[0]]);nn=nn/mp.norm(nn);vv=mp.matrix([-nn[1],nn[0]])
  def ht(w):return al*(w[0]+w[1])/2-mp.sqrt(c*w[0]*w[1])
  def chart(s):return (ht(nn+s*t*vv)-(p.T*(nn+s*t*vv))[0])/t**2
  rho=mp.sqrt(c)/(4*((nn[0]*nn[1])**mp.mpf('1.5')))
  for s in [-s0,s0]:
   f=chart(s);fp=mp.diff(chart,s);assert abs(f-rho*s*s/2)<mp.mpf('0.1')*eta
   charts.append(dict(t=ts,s=mp.nstr(s,25),chart=mp.nstr(f,35),slope=mp.nstr(fp,35),boundary_normal=mp.nstr(f-s*fp,35)))
 out=dict(python=platform.python_version(),mpmath=mp.__version__,precision=80,exact_pair_compatibility_certificates=cert,corner_regressions=rows,support_chart_regressions=charts,seconds=time.monotonic()-st,scope='exact fixed-scale four-point certificates; numerical regression only for universal cubic theorem')
 (ROOT/'independent_cubic_python.json').write_text(json.dumps(out,indent=2)+'\n');print({'exact_certificates':len(cert),'geometry_regressions':len(rows),'chart_regressions':len(charts),'seconds':out['seconds']})
if __name__=='__main__':run()
