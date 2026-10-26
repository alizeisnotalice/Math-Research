#!/usr/bin/env python3
"""Resume the immutable213-family fine prefix, with transactional slab checkpoints."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import prod
import sys,os,time,json,hashlib,resource,threading,sqlite3,gzip
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent;OLD=ROOT/'lower_round025'
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
ZERO=Q(0)
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,separators=(',',':'))+'\n');t.replace(p)
def rt(start):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def decode(d):
    z=[ZERO]*25
    for ij,c in d.items():i,j=map(int,ij.split(','));z[5*i+j]=Q(c)
    return z
def encode(z):return {f'{i//5},{i%5}':str(c) for i,c in enumerate(z) if c}
def mul(x,y):
    z=[ZERO]*25
    for i,c in enumerate(x):
        if not c:continue
        ii,ij=divmod(i,5)
        for j,d in enumerate(y):
            if not d:continue
            ji,jj=divmod(j,5);a,b=ii+ji,ij+jj
            z[(a%5)*5+b%5]+=c*d*2**(a//5)*3**(b//5)
    return z
def addto(target,x):
    for i,c in enumerate(x):
        if c:target[i]+=c

def main(stage):
    start=time.monotonic();reg=load(P/'registration.json');assert sha(Path(__file__))==reg['source_sha256']
    for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
    dbpath=Path(reg['scratch_database']);db=sqlite3.connect(dbpath);db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=FULL');db.execute('PRAGMA cache_size=-8192')
    db.execute('CREATE TABLE IF NOT EXISTS hist(mask TEXT PRIMARY KEY, coefficients TEXT NOT NULL)')
    db.execute('CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL)')
    def meta():return dict(db.execute('SELECT key,value FROM meta'))
    def put(k,v):db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)',(k,str(v)))
    state=meta()
    if not state:
        assert db.execute('SELECT COUNT(*) FROM hist').fetchone()[0]==0
        old=load(OLD/'tail_checkpoint.json')
        assert old['completed_first_tail_slabs']==11 and old['total_first_tail_slabs']==26 and old['nodes_this_process']==193336
        with db:
            db.executemany('INSERT INTO hist VALUES (?,?)',((m,json.dumps(c,separators=(',',':'))) for m,c in old['histogram'].items()))
            put('source_sha256',reg['source_sha256']);put('old_prefix_sha256',sha(OLD/'tail_checkpoint.json'));put('completed_slabs',11);put('cells_done',193336);put('total_slabs',26)
        del old;state=meta()
    assert state['source_sha256']==reg['source_sha256'] and state['old_prefix_sha256']==reg['old_prefix_sha256']
    F,u,v,_,_,_=field();geom=load(OLD/'tail_geometry.json');assert len(geom['families'])==213
    if stage=='resume':
        axes=[[F(a=decode(x)) for x in axis] for axis in geom['axes']]
        boxes=[(r['family_id'],[F(a=decode(z)) for z in r['tail_box']]) for r in geom['families']]
        segments=[]
        for k,axis in enumerate(axes):
            segments.append([((hi-lo).a,sum(1<<fid for fid,b in boxes if b[2*k]<=lo and hi<=b[2*k+1])) for lo,hi in zip(axis,axis[1:])])
        assert list(map(len,segments))==[26]*4
        triples=[]
        for (dy,my),(dz,mz),(dw,mw) in product(*segments[1:]):triples.append((mul(mul(dy,dz),dw),my&mz&mw))
        prefix=int(state['completed_slabs']);newcells=0
        print({'stage':'resume_original213_fine_unfinished_slabs','pid':os.getpid(),'initial_prefix_slabs':prefix,'remaining_cells':(26-prefix)*17576,'new_cell_cap':200000,'old_prefix_sha256':state['old_prefix_sha256'],'runtime':rt(start)},flush=True)
        for ix in range(prefix,26):
            if newcells+17576>200000:break
            dx,mx=segments[0][ix];delta={}
            for yzww,myzw in triples:
                mask=str(mx&myzw);value=mul(dx,yzww)
                if mask not in delta:delta[mask]=value
                else:addto(delta[mask],value)
            # Atomic transaction: either the entire new xslab and its prefix
            # commit together, or SQLite leaves the preceding exact prefix.
            with db:
                for m,z in delta.items():
                    oldrow=db.execute('SELECT coefficients FROM hist WHERE mask=?',(m,)).fetchone()
                    if oldrow:
                        oldz=decode(json.loads(oldrow[0]));addto(oldz,z);z=oldz
                    db.execute('INSERT OR REPLACE INTO hist VALUES (?,?)',(m,json.dumps(encode(z),separators=(',',':'))))
                put('completed_slabs',ix+1);put('cells_done',(ix+1)*17576)
            newcells+=17576;prefix=ix+1
            save(P/'checkpoint.json',{'status':'COMPLETE_full213_fine_histogram' if prefix==26 else 'EXACT_fine_slab_prefix_saved','completed_slabs':prefix,'total_slabs':26,'cells_done':prefix*17576,'new_cells_this_process':newcells,'old_cells_recomputed':0,'source_sha256':reg['source_sha256'],'scratch_database':str(dbpath),'runtime':rt(start)})
        receipt={'status':'COMPLETE_full213_fine_histogram' if prefix==26 else 'FINITE_NEW_CELL_CAP_prefix_saved','initial_slabs':int(state['completed_slabs']),'completed_slabs':prefix,'new_cells':newcells,'total_cells_done':prefix*17576,'old_cells_recomputed':0,'histogram_signatures':db.execute('SELECT COUNT(*) FROM hist').fetchone()[0],'runtime':rt(start),'source_sha256':reg['source_sha256']}
        save(P/f'batch_from_slab{state["completed_slabs"]}.json',receipt);print(receipt,flush=True)
    elif stage=='export':
        assert int(state['completed_slabs'])==26
        total=[ZERO]*25
        for coeff, in db.execute('SELECT coefficients FROM hist'):addto(total,decode(json.loads(coeff)))
        assert F(a=total)==((1+u)/4)**4
        artifact=P/'full_tail_histogram.json.gz';tmp=artifact.with_suffix('.gz.tmp')
        with tmp.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as out:
                out.write(b'{"status":"COMPLETE_full213_fine_static_histogram","histogram":{');first=True
                for m,coeff in db.execute('SELECT mask,coefficients FROM hist ORDER BY length(mask),mask'):
                    if not first:out.write(b',')
                    first=False;out.write(json.dumps(m).encode()+b':'+coeff.encode())
                out.write(b'},"total_physical_tail_volume":'+json.dumps(encode(total),separators=(',',':')).encode()+b'}\n')
        tmp.replace(artifact)
        out={'status':'COMPLETE_full213_fine_static_histogram','families':213,'cells':26**4,'old_successful_cells_reused':193336,'new_cells':263640,'old_cells_recomputed':0,'signatures':db.execute('SELECT COUNT(*) FROM hist').fetchone()[0],'total_physical_tail_volume':encode(total),'artifact':artifact.name,'artifact_bytes':artifact.stat().st_size,'artifact_sha256':sha(artifact),'old_prefix_sha256':state['old_prefix_sha256'],'geometry_sha256':sha(OLD/'tail_geometry.json'),'source_sha256':reg['source_sha256'],'runtime':rt(start)}
        save(P/'static_result.json',out);print(out,flush=True)
    else:raise ValueError('stage must be resume or export')
    db.close()
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    if len(sys.argv)!=2:raise SystemExit('usage: full_static.py resume|export')
    main(sys.argv[1])
