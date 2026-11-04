import mpmath as m, json, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def task(arg):
 m.mp.dps=70
 al,fr=map(m.mpf,arg); w=1-(1-al)*fr; W=1/w; s=(w+W)/2; A=s-al; X=(1/al-W)/2; Y=(1/al-w)/2
 def val(x,y):
  lengths=[W-2*x,w-2*y,W+2*x,w+2*y]
  return max(max(0,min(W/2,x+l/2)-max(-W/2,x-l/2))*max(0,min(w/2,y+l/2)-max(-w/2,y-l/2))/l for l in lengths if l>0)
 err=0
 for j in range(1,40):
  x=A*m.mpf(j)/40
  y=(A-x)/(2*al/w-1) if x>=X else A-(2*al/W-1)*x
  err=max(err,abs(val(x,y)-al))
 # independent high precision midpoint boundary quadrature, split at actual vertex
 area=0
 for l,r in [(m.mpf(0),X),(X,A)]:
  for j in range(8):
   x=l+(r-l)*(m.mpf(j)+m.mpf('.5'))/8; low=0; high=w/2
   for _ in range(210):
    mid=(low+high)/2
    if val(x,mid)>al: low=mid
    else: high=mid
   area+=(r-l)/8*(low+high)/2
 exact=2*(s-al)*(1/al-s)
 return [str(al),str(w),str(err),str(abs(4*area-exact))]
if __name__=='__main__':
 start=time.monotonic(); jobs=[(a,f) for a in ['.707107','.71','.8','.9','.99','.99999'] for f in ['0','.1','.5','.9','.99','.9999']]
 with ProcessPoolExecutor(max_workers=8) as pool: rows=list(pool.map(task,jobs))
 m.mp.dps=70
 out={'mpmath_version':m.__version__,'precision_digits':70,'workers':8,'bisection_steps':210,'cases':len(rows),'rows':rows,'max_boundary_residual':str(max(m.mpf(r[2]) for r in rows)),'max_area_error':str(max(m.mpf(r[3]) for r in rows)),'elapsed_seconds':time.monotonic()-start,'scope':'finite high precision evidence; split midpoint integration exact for predicted piecewise linear boundary, boundary obtained from original overlap'}
 (Path(__file__).resolve().parents[1]/'results'/'rectangle_python.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k!='rows'})
