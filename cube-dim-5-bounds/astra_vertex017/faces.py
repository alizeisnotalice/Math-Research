"""All faces by exact facet-mask intersections, then codimension-one incidences."""
from pathlib import Path
from collections import Counter
import hashlib,time,os,json
from exact import save,read,sha,dumps,mask_of,ids,rank
P=Path(__file__).resolve().parent;CACHE=Path('/tmp/astra_vertex017')

def main():
    began=time.monotonic();input_path=CACHE/'contacts.json.gz';data=read(input_path)
    source={str(p):sha(p) for p in [Path(__file__),P/'exact.py',input_path]}
    fingerprint=hashlib.sha256(dumps(source)).hexdigest();cp=CACHE/'faces_checkpoint.json.gz'
    vertices=data['vertices'];full=(1<<len(vertices))-1;facets=[mask_of(f['vertices']) for f in data['facets']]
    state=read(cp) if cp.exists() else {'fingerprint':fingerprint,'phase':'closure','next':0,'queue':[hex(full)],'dimensions':[],'incidences':[]}
    assert state['fingerprint']==fingerprint
    def checkpoint():save(cp,state)
    if state['phase']=='closure':
        queue=list(map(lambda x:int(x,16),state['queue']));known=set(queue)|{0};index=state['next']
        while index<len(queue):
            F=queue[index]
            for G in sorted({F&H for H in facets}-{0,F}):
                if G not in known:known.add(G);queue.append(G)
            index+=1
            if index%512==0 or time.monotonic()-began>120:
                state.update(next=index,queue=[hex(m) for m in queue]);checkpoint()
                print(json.dumps({'pid':os.getpid(),'phase':'closure','processed':index,'discovered':len(known),'seconds':time.monotonic()-began}),flush=True)
            if time.monotonic()-began>120:return
        ordered=sorted(known,key=lambda m:(m.bit_count(),m));state.update(phase='ranks',next=0,queue=[hex(m) for m in ordered]);checkpoint()
    masks=[int(x,16) for x in state['queue']]
    if state['phase']=='ranks':
        for i in range(state['next'],len(masks)):
            F=masks[i];d=5 if F==full else (-1 if F==0 else rank([vertices[j] for j in ids(F)])-1)
            assert -1<=d<=5;state['dimensions'].append(d);state['next']=i+1
            if (i+1)%512==0 or time.monotonic()-began>120:checkpoint()
            if time.monotonic()-began>120:return
        state.update(phase='incidences',next=0);checkpoint()
    index_of={m:i for i,m in enumerate(masks)};dims=state['dimensions']
    if state['phase']=='incidences':
        for i in range(state['next'],len(masks)):
            F=masks[i];d=dims[i]
            if d==-1:children=[]
            elif d==0:children=[index_of[0]]
            else:
                intersections={F&H for H in facets}-{F}
                assert all(G in index_of for G in intersections)
                children=sorted(index_of[G] for G in intersections if dims[index_of[G]]==d-1)
                assert children
            state['incidences'].append(children);state['next']=i+1
            if (i+1)%512==0 or time.monotonic()-began>120:checkpoint()
            if time.monotonic()-began>120:return
        state.update(phase='complete');checkpoint()
    assert state['phase']=='complete'
    # The complete face lattice is transient; its fingerprint authenticates pulling input.
    result={'fingerprint':fingerprint,'contacts_sha256':sha(input_path),'masks':state['queue'],'dimensions':dims,'children':state['incidences']}
    path=CACHE/'faces.json.gz';save(path,result)
    bydim=Counter(dims);fvector=[bydim[d] for d in range(5)];assert bydim[-1]==bydim[5]==1
    assert sum((-1)**d*fvector[d] for d in range(5))==2
    summary={'status':'complete_exact_face_lattice','faces_including_empty_and_polytope':len(masks),'boundary_f_vector':fvector,
             'incidence_counts_by_parent_dimension':{str(d):sum(len(ch) for dim,ch in zip(dims,state['incidences']) if dim==d) for d in range(6)},
             'Euler_characteristic_boundary':2,'faces_path':str(path),'faces_sha256':sha(path),'source_hashes':source,
             'pid':os.getpid(),'seconds_this_process':time.monotonic()-began}
    save(P/'faces_result.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
