"""Exactly certify the already-scanned complement/quarter-mass structural point."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
from scipy.stats import qmc
import round006_search as eng
from round006_certify import verify
ROOT=Path(__file__).resolve().parent

def main():
 start=time.monotonic();seed=eng.SEED;eps=F(seed['weights'][1]);theta=list(map(F,seed['phases'][1]));theta2=[1-x for x in theta];z=np.r_[np.zeros(3),np.array(list(map(float,theta))),np.array(list(map(float,theta2))),np.log(float(eps)),np.log(float(eps/4))];keys=sorted(set(eng.OLD)|set(eng.CLUSTERS)|set(eng.capture(z,qmc.Sobol(3,scramble=True,seed=20261025602).random_base2(12))));b,active=eng.boxes(z,keys);saved=sorted(set(map(tuple,active))|set(eng.OLD))
 inp={'periods':seed['periods'],'phases':[['0','0','0'],seed['phases'][1],[str(x) for x in theta2]],'weights':['1',str(eps),str(eps/4)],'keys':[list(k) for k in saved]};fp=ROOT/'round006_structural_input.json';fp.write_text(json.dumps(inp,indent=2)+'\n');out=verify(inp,json.loads((ROOT/'round006_target.json').read_text())['ratio_interval'][1]);out['input_sha256']=hashlib.sha256(fp.read_bytes()).hexdigest();out['attribution']='original P, theta1, eps1 fixed exactly; only second satellite theta2=1-theta1, eps2=eps1/4 is added';out['elapsed_seconds']=time.monotonic()-start
 assert out['lower_above_full_007c_upper'];out.pop('radius_records');(ROOT/'round006_structural_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
if __name__=='__main__':main()
