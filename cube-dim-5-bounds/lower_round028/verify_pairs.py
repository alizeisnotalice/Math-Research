#!/usr/bin/env python3
"""Independent2D+2D mask factorization,422500 actual pairs in finite batches."""
from pathlib import Path
from fractions import Fraction as Q
import sys,os,time,json,gzip,sqlite3
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from full_static import ENV,load,save,sha,rt,decode,encode,mul,addto,ZERO

def stream_histogram(path):
    decoder=json.JSONDecoder();marker='"histogram":{'
    with gzip.open(path,'rt') as f:
        buffer=''
        while marker not in buffer:
            block=f.read(65536)
            if not block:raise ValueError('missing histogram marker')
            buffer+=block
        buffer=buffer.split(marker,1)[1]
        while True:
            buffer=buffer.lstrip()
            if not buffer:
                buffer=f.read(65536)
                if not buffer:raise ValueError('truncated histogram')
                continue
            if buffer[0]=='}':break
            try:
                mask,pos=decoder.raw_decode(buffer)
                while pos<len(buffer) and buffer[pos].isspace():pos+=1
                if pos==len(buffer):raise IndexError()
                assert buffer[pos]==':'
                coeff,end=decoder.raw_decode(buffer,pos+1)
            except (json.JSONDecodeError,IndexError):
                block=f.read(65536)
                if not block:raise ValueError('truncated mask coefficient')
                buffer+=block;continue
            yield mask,coeff
            buffer=buffer[end:].lstrip()
            while not buffer:
                buffer=f.read(65536)
                if not buffer:raise ValueError('truncated item separator')
            if buffer[0]==',':buffer=buffer[1:]
            elif buffer[0]=='}':break
            else:raise ValueError('invalid histogram separator')

def main(stage):
    start=time.monotonic();reg=load(P/'pair_verification_registration.json');assert sha(Path(__file__))==reg['source_sha256']
    for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
    plan=load(P/'pair_plan.json');left,right=[sorted((int(m),decode(c)) for m,c in h.items()) for h in plan['histograms']]
    assert len(left)==len(right)==650 and len(left)*len(right)==422500
    db=sqlite3.connect(reg['scratch_database']);db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=FULL');db.execute('PRAGMA cache_size=-8192')
    db.execute('CREATE TABLE IF NOT EXISTS hist(mask TEXT PRIMARY KEY, coefficients TEXT NOT NULL)');db.execute('CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL)')
    def meta():return dict(db.execute('SELECT key,value FROM meta'))
    def put(k,v):db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)',(k,str(v)))
    state=meta()
    if not state:
        assert db.execute('SELECT COUNT(*) FROM hist').fetchone()[0]==0
        with db:put('source_sha256',reg['source_sha256']);put('completed_left_rows',0)
        state=meta()
    assert state['source_sha256']==reg['source_sha256']
    initial=int(state['completed_left_rows']);prefix=initial;newpairs=0
    if stage=='resume':
        print({'stage':'full213_independent_pair_factorization','pid':os.getpid(),'initial_left_rows':initial,'actual_total_pairs':422500,'per_process_cap':reg['max_new_pairs_per_process']},flush=True)
        for i in range(initial,len(left)):
            if newpairs+len(right)>reg['max_new_pairs_per_process']:break
            lm,lw=left[i];delta={}
            for rm,rw in right:
                m=str(lm&rm);value=mul(lw,rw)
                if m not in delta:delta[m]=value
                else:addto(delta[m],value)
            with db:
                for m,z in delta.items():
                    row=db.execute('SELECT coefficients FROM hist WHERE mask=?',(m,)).fetchone()
                    if row:old=decode(json.loads(row[0]));addto(old,z);z=old
                    db.execute('INSERT OR REPLACE INTO hist VALUES (?,?)',(m,json.dumps(encode(z),separators=(',',':'))))
                put('completed_left_rows',i+1)
            newpairs+=len(right);prefix=i+1
            save(P/'pair_verification_checkpoint.json',{'status':'COMPLETE_pair_histogram' if prefix==650 else 'EXACT_left_signature_row_prefix_saved','completed_left_rows':prefix,'total_left_rows':650,'completed_pairs':prefix*650,'new_pairs_this_process':newpairs,'source_sha256':reg['source_sha256'],'runtime':rt(start)})
        receipt={'status':'COMPLETE_pair_histogram' if prefix==650 else 'FINITE_NEW_PAIR_CAP_prefix_saved','initial_left_rows':initial,'completed_left_rows':prefix,'new_pairs_this_process':newpairs,'completed_pairs':prefix*650,'source_sha256':reg['source_sha256'],'runtime':rt(start)}
        save(P/f'pair_batch_from_{initial}.json',receipt);print(receipt,flush=True)
    elif stage=='compare':
        assert int(state['completed_left_rows'])==650
        total=[ZERO]*25;count=0
        for m,c in stream_histogram(P/'full_tail_histogram.json.gz'):
            row=db.execute('SELECT coefficients FROM hist WHERE mask=?',(m,)).fetchone();assert row is not None
            assert json.loads(row[0])==c,m
            addto(total,decode(c));count+=1
        assert count==db.execute('SELECT COUNT(*) FROM hist').fetchone()[0]==139843
        static=load(P/'static_result.json');assert encode(total)==static['total_physical_tail_volume']
        result={'status':'PASS_entire_full213_histogram_independent_pair_factorization','actual_pair_count':422500,'signature_counts':[650,650],'all_exact_mask_coefficients_equal':139843,'full_physical_tail_volume_equal':True,'new_parameter_or5D_union_evaluations':0,'arithmetic_boundary':'pairfactorizationandfineslabproducer shareFractionandQ(u,v)convolutionkernel;geometricintegrationgroupingsandmaskintersectionsaredifferent','source_sha256':reg['source_sha256'],'static_artifact_sha256':sha(P/'full_tail_histogram.json.gz'),'runtime':rt(start)}
        save(P/'pair_verification_result.json',result);print(result,flush=True)
    else:raise ValueError('resume or compare')
    db.close()
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    if len(sys.argv)!=2:raise SystemExit('usage: verify_pairs.py resume|compare')
    main(sys.argv[1])
