"""Enclose whole dyadic query cells to connect ideal 53-bit sampling to continuous queries."""
from pathlib import Path
from fractions import Fraction as F
import json,time,sys
import numpy as np
import interval_recheck as ir
import probe as core
HERE=Path(__file__).resolve().parent

def cell_overlap(x,L,left,right):
 # Generator.random uses the 53-bit grid. x+2^-53 is exact on that grid.
 assert np.all(np.ldexp(x,53)==np.floor(np.ldexp(x,53)))
 xhigh=x+2.**-53
 h=ir.rational(L/2)
 wl=(ir.dn(x-h[1]),ir.up(xhigh-h[0]));wr=(ir.dn(x+h[0]),ir.up(xhigh+h[1]))
 lo=np.zeros_like(x);hi=lo.copy()
 for shift in range(-2,3):
  al,au=ir.rational(left+shift);bl,bu=ir.rational(right+shift)
  ol=np.maximum(0.,ir.dn(np.minimum(wr[0],bl)-np.maximum(wl[1],au)))
  oh=np.maximum(0.,ir.up(np.minimum(wr[1],bu)-np.maximum(wl[0],al)))
  lo,hi=ir.add((lo,hi),(ol,oh))
 return lo,hi

def main():
 outpath=HERE/'CELL_RESULTS.json';pf=HERE/'CELL_PLAN.json';assert not outpath.exists()and not pf.exists();start=time.time()
 plan=dict(script_SHA256=core.sha(__file__),parent_SHA256=core.sha(ir.__file__),core_SHA256=core.sha(core.__file__),
   query_cells='[x_i,x_i+2^-53]; x_i on exact53bit grid',alpha_each_exact='1/720',
   prior_interval_results_SHA256=core.sha(HERE/'INTERVAL_RESULTS.json'),
   scope='Same frozen samples/selectors; whole query-cell enclosure, no new RNG draws beyond replay; ideal iid uniform53bits assumption explicit, not a proof of PRNG randomness')
 pf.write_text(json.dumps(plan,indent=2)+'\n');ir.overlap=cell_overlap;cases=[]
 for folder in [HERE,HERE/'balanced/n7',HERE/'balanced/n15']:
  cases.append(ir.run_case(folder));print(json.dumps(cases[-1]),flush=True)
 out=dict(cases=cases,plan_SHA256=core.sha(pf),query_cell_width_exact='1/9007199254740992',
   same_frozen_queries=True,ideal_iid_uniform_bits_assumed=True,PRNG_randomness_proved=False,
   deterministic_continuum_integral_certificate=False,Cn_log_gain_proved=False,
   python=sys.version,numpy=np.__version__,seconds=time.time()-start)
 outpath.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
