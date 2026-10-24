#!/usr/bin/env python3
"""Fixed 3 addRows probes, no solves or kernel evaluations, 10s hard wrapper."""
import json,sys,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
if '--worker' not in sys.argv:
    out=ROOT/'interface_probe.json'
    if out.exists():raise SystemExit('No repeat of existing probe.')
    start=time.monotonic()
    p=subprocess.run([sys.executable,__file__,'--worker'],capture_output=True,text=True,timeout=10)
    result={'budget_seconds':10,'wall_seconds':time.monotonic()-start,'exit_code':p.returncode,'stderr':p.stderr,'worker':json.loads(p.stdout)}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));raise SystemExit
import highspy,numpy as np
records=[]
for label,indices,values in [('ordinary',[0,1],[1.,.25]),('tiny_entry',[0,1],[1.,1e-12]),('invalid_column',[0,2],[1.,.25])]:
    h=highspy.Highs();h.setOptionValue('output_flag',False)
    vs=h.addVars(2,np.array([-1.,-1.]),np.array([1.,1.]))
    rs=h.addRows(1,np.array([0.]),np.array([highspy.kHighsInf]),2,np.array([0],dtype=np.int32),np.array(indices,dtype=np.int32),np.array(values))
    lp=h.getLp();a=lp.a_matrix_
    records.append({'case':label,'addVars':str(vs),'addRows':str(rs),'num_rows':lp.num_row_,'num_cols':lp.num_col_,
                    'matrix_format':str(a.format_),'start':list(a.start_),'indices':list(a.index_),'values':list(a.value_),
                    'row_lower':list(lp.row_lower_),'row_upper':list(lp.row_upper_)})
print(json.dumps({'version':highspy.Highs().version(),'addRows_calls':3,'solver_run_calls':0,'kernel_evaluations':0,'cases':records}))
