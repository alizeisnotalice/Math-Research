"""Exact finite-input margins and two actual witness cubes per dimension.
Uses the archived root brackets, not a root search. Intersection averages use
an independent exact interval-intersection sum; the geometric proof excludes
the remaining lattice indices and transfers the witnesses into E_N.
"""
from fractions import Fraction as F
from pathlib import Path
from math import prod
import argparse,json
from compute import evaluate

def frac(x):return str(x)
def exact(x):return {'numerator_hex':hex(x.numerator),'denominator_hex':hex(x.denominator)}
def overlap(c,x,r,d):return max(F(0),min(x+r,c+d)-max(x-r,c-d))
def make(record):
 n=record['n'];k=F(record['kappa']);u=F(record['u']);delta=u/(1+2*u);r0=1/(1+2*u);e=k*u**4;t=u**n;a=(2*u)**n;v=2*u;m=(2*delta)**n;det=1-e**2
 N=(u.denominator//u.numerator)**9
 assert N==u**-9 and N>3
 fresh,(gl,gu)=evaluate((n,k,u.denominator),record['root']);assert fresh['status']=='certified'
 denom=(det-m)*(1-m)*a*v**8;bu=(m*e**2-gl*denom)/(1-m);assert 0<bu<m
 lower=1+(F(2*N-5,2*N))**n*(m-bu)/(det-m)
 finite_margin=lower-1/(1-m);assert finite_margin>0
 hx=e**2/(2*n*u**2*r0);assert 0<hx<e/2<delta
 points=[[F(0)]*n,[delta-hx,-delta+e/2]+[F(0)]*(n-2)]
 witnesses=[]
 for x in points:
  s=delta-max(map(abs,x));r=r0+s;assert r+delta==1-max(map(abs,x)) and r>2*delta
  total=F(0);positive=[]
  for i in range(-2,3):
   for j in range(-2,3):
    c=[F(i)+e*j,e*i+F(j)];mass=overlap(c[0],x[0],r,delta)*overlap(c[1],x[1],r,delta)*(2*delta)**(n-2)
    total+=mass
    if mass:positive.append([i,j])
  Pp=max(F(0),e+s-(delta-x[0]))*max(F(0),e+s-(delta+x[1]))
  Pm=max(F(0),e+s-(delta+x[0]))*max(F(0),e+s-(delta-x[1]))
  assert total==m+(2*delta)**(n-2)*(Pp+Pm)
  defect=total/(2*r)**n;margin=t-defect;assert margin>0
  # In index coordinates this proves the same actual cube lies in Omega_N.
  assert (delta+r)/(1-e)<2<N-3
  witnesses.append({'center':[frac(z) for z in x],'radius':frac(r),'positive_hole_indices_first_two':positive,'defect_average':frac(defect),'input_average':frac(1-defect),'average_over_threshold_margin':frac(margin)})
 return {'n':n,'kappa':str(k),'u':str(u),'finite_input':{'definition':'E_N=(R^n minus periodic T_e lattice cube pores) intersect T_e[-N,N]^n','N':N,'pore_halfwidth':str(delta),'shear':str(e)},'finite_ratio_minus_Ln_lower_exact':exact(finite_margin),'finite_normalized_gap_lower_exact':exact(finite_margin/(a*v**8)),'witnesses':witnesses}

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args();records=json.loads(args.input.read_text())['cases']
 selected=[r for r in records if F(r['u'])==F(1,1024) and F(r['kappa'])==F(16*r['n'],3)]
 assert sorted(r['n'] for r in selected)==list(range(2,17))
 rows=[make(r) for r in selected];out={'scope':'15 explicit finite inputs and 30 exact sample witness averages; not all-input optimality','all_finite_margins_positive':True,'all_sample_witness_margins_positive':True,'cases':rows}
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'finite_inputs':len(rows),'sample_witnesses':sum(len(r['witnesses']) for r in rows),'finite_margins_positive':True,'witness_margins_positive':True}))
if __name__=='__main__':main()
