"""Exact arithmetic exploration of sparse versus dense fixed-window budgets."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
HERE=Path(__file__).resolve().parent

def main():
 target=HERE/'SCALE_BUDGET.json';assert not target.exists();rows=[]
 for Q in [1,2,4,16,64,256]:
  s=sum(F(1,4*q*(q+1))for q in range(1,Q+1));assert s==F(Q,4*(Q+1))
  rows.append(dict(Q=Q,sparse_sum_exact=str(s),tail_to_quarter_exact=str(F(1,4)-s),sparse_sum_float=float(s)))
 dense=[]
 for J in [4,16,64,256,1024]:
  L=[F(1)+F(j+1,J+1)for j in range(J)];s=sum(F(1,4*int(v)*(int(v)+1))for v in L);assert s==F(J,8)
  dense.append(dict(J=J,all_scales_in='(1,2)',termwise_sum_exact=str(s)))
 out=dict(sparse=rows,dense=dense,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Exact budget arithmetic after actual joint integral probes; sparse-tail sum finite, dense within-unit-interval sum diverges. No C_n estimate.')
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
