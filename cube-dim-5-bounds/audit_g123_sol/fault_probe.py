#!/usr/bin/env python3
"""Before any H: bad pilot header; prove idle-worker cleanup returns."""
import sys
sys.dont_write_bytecode=True
import multiprocessing as mp,os,time
import common as C
from numeric import worker,release

def main():
    began=time.monotonic();ctx=mp.get_context('spawn');tq=ctx.Queue();rq=ctx.Queue();stop=ctx.Event();acks=[ctx.Queue() for _ in range(4)]
    pool=[ctx.Process(target=worker,args=(i,tq,rq,acks[i],stop)) for i in range(4)]
    try:
        for p in pool:p.start()
        ready=[rq.get(timeout=20) for _ in pool];assert all(m['kind']=='ready' for m in ready)
        # raw_gate fails on the header before P.values or any logarithm.
        tq.put((0,[dict(index=0,entry={},task={})]));msg=rq.get(timeout=20)
        assert msg['kind']=='fault' and "raw_gate" in msg['error'] and "KeyError" in msg['error']
    finally:exits=release(pool,tq,acks,stop,True)
    assert len(exits)==4 and all(not x['alive'] and x['exitcode']==0 for x in exits)
    result=dict(status='bad_pilot_before_H_releases_all_four_idle_workers',pid=os.getpid(),ready=ready,pool_exits=exits,seconds=time.monotonic()-began,new_H=0,new_gradient=0,new_log=0,shutdown_code_sha256=C.sha(C.HERE/'numeric.py'))
    C.write(C.HERE/'fault_probe.json',result);print(__import__('json').dumps(result),flush=True)

if __name__=='__main__':main()
