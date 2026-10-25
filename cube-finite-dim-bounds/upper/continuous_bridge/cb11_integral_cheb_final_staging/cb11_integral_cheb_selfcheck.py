"""Author algebra/output consistency check; not independent audit or new nodes."""
from pathlib import Path
from fractions import Fraction as F
import json,runpy,hashlib,argparse
from flint import arb,arb_mat
OUT=Path(__file__).resolve().parent
ns=runpy.run_path(str(OUT/'cb11_integral_cheb_check.py'));N=ns['N'];aq=ns['aq'];unpack=ns['unpack'];pack=ns['pack'];exact=ns['exact']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--parts',type=int,default=0);args=ap.parse_args()
 BF=ns['BF'];MU=ns['MU'];A=ns['A'];P=ns['P'];degree=N-1
 # Exact Bernstein mean agrees with exact shifted-power integral for every T_k.
 for k in range(N):assert sum((BF[j][k] for j in range(N)),F(0))/N==MU[k]
 # Arb encloses all exact DCT inverse identities.
 E=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos() for k in range(N)] for j in range(N)]);I=A*E
 for j in range(N):
  for k in range(N):assert (I[j,k]-(1 if j==k else 0)).contains(0) and I[j,k].rad()<arb(2)**(-300)
 rows=[];files=[]
 if args.parts:
  for part in range(args.parts):
   path=OUT/f'cb11_integral_cheb_part_{part:02}.json';r=json.loads(path.read_text());assert r['config']==ns['CONFIG'] and r['checker_sha256']==hashlib.sha256((OUT/'cb11_integral_cheb_check.py').read_bytes()).hexdigest()
   expected={(i,j) for i in range(128) for j in range(i+1) if (i*(i+1)//2+j)%args.parts==part}
   assert len(r['patches'])==len(expected) and {tuple(v['patch']) for v in r['patches']}==expected
   rows+=r['patches'];files.append(path)
  assert len(rows)==8256
 else:
  path=OUT/'cb11_integral_cheb_four_result.json';r=json.loads(path.read_text());assert r['config']==ns['CONFIG'];rows=r['patches'];files=[path]
  assert [tuple(v['patch']) for v in rows]==ns['FOUR']
 for r in rows:
  i,j=r['patch'];L=unpack(r['L']);U=unpack(r['U']);m=unpack(r['mean_H_upper']);assert 0<=L<=m<=U and U>0
  slope=(U+(U*L).sqrt()+L)/(U.sqrt()+L.sqrt());c=L*L.sqrt()+slope*(m-L);w=1 if i==j else 2
  bound=aq(F(r['IH_contribution_exact']));assert bound>=(c*w/128**2).upper() and bound==unpack(r['IH_contribution_upper'])
 out={'status':'PASS author exact means/DCT/output chord consistency','new_point_evaluations':0,'patches_checked':len(rows),'full_coverage':bool(args.parts),'I_claim':'only main full aggregate can claim I, conditional on its audited source/time inputs','input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
 (OUT/('cb11_integral_cheb_all_selfcheck.json' if args.parts else 'cb11_integral_cheb_four_selfcheck.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
