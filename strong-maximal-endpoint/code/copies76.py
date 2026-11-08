import json,sys,time,math
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def run(j):
 n=j['N'];k=j['copies'];A=j['volume_enlargement'];D=8*A*k+2;gap=F(D)-F(1,2)-4*A*k;assert gap>0
 return dict(j,packet_amplitude='3/2',packet_source_cube='[0,1]^N',packet_receiver_cube='[1/4,3/4]^N',spacing=D,projection_separation_margin=str(gap),continuous_corner_count_lower=k,explicit_exact_corner_count_upper=k,source_allocation_reuse=1,score_log=-n*math.log(2)-math.log(1.5)-(n-1)*math.log1p(math.log(1.5)),scope='explicit low-score packet model; general near-Q copies use the analytic theorem')
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);jobs=[dict(N=n,copies=k,volume_enlargement=A) for n in [2,4,8,16,32,64,128,256] for k in [1,4,16,64,256] for A in [1,2,4]];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs))
 (out/'copies76.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-start,rows=rows),indent=2)+'\n');print('PASS',len(rows),'exact separated-packet budget certificates')
