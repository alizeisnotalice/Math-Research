"""Exhaustive private-cell subset search on finite rectangular tori."""
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path
import argparse,json,math,os,subprocess,tempfile,time


def catalog(shape,weights):
    points=list(product(*(range(m) for m in shape)))
    axes=[{tuple(sorted((o+i)%m for i in range(s))) for o in range(m) for s in range(1,m+1)} for m in shape]
    boxes={sum(1<<i for i,point in enumerate(points) if all(x in a for x,a in zip(point,rect))) for rect in product(*axes)}
    mass=[math.prod(axis[c] for axis,c in zip(weights,point)) for point in points]
    return {b:sum(v for i,v in enumerate(mass) if b>>i&1) for b in boxes},mass


def main():
    ap=argparse.ArgumentParser();ap.add_argument('outdir',type=Path);ap.add_argument('--shapes',nargs='*');args=ap.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    cases=[(2,2),(2,3),(3,3),(2,2,2),(3,4),(4,4),(2,3,3),(2,2,2,2),(2,3,4),(2,2,2,3),(3,3,3)]
    if args.shapes:cases=[tuple(map(int,x.split('x'))) for x in args.shapes]
    with tempfile.TemporaryDirectory() as tmp:
        exe=Path(tmp)/'private101';subprocess.run(['c++','-std=c++20','-O3',str(Path(__file__).with_name('private101.cpp')),'-o',str(exe)],check=True)
        for shape in cases:
            weights=[[1]*m for m in shape];boxes,mass=catalog(shape,weights);n=len(mass);end=1<<n
            data=f'{n} {len(boxes)}\n'+''.join(f'{k} {v}\n' for k,v in sorted(boxes.items()));jobs=[(end*i//12,end*(i+1)//12) for i in range(12)]
            def search(pair):
                return list(map(int,subprocess.check_output([str(exe),*map(str,pair)],input=data,text=True).split()))
            start=time.monotonic()
            with ThreadPoolExecutor(max_workers=12) as pool:parts=list(pool.map(search,jobs))
            assert [tuple(x[:2]) for x in parts]==jobs
            best=max(x[2] for x in parts);arg=next(x[3] for x in parts if x[2]==best)
            chosen=[max((b for b in boxes if b&arg==1<<i),key=boxes.get) for i in range(n) if arg>>i&1];union=0
            for b in chosen:union|=b
            assert union==end-1 and sum(boxes[b] for b in chosen)==best
            if n<=12:
                oracle=max(sum(max(v for b,v in boxes.items() if b&p==1<<i) for i in range(n) if p>>i&1) for p in range(end));assert oracle==best
            record=dict(shape=shape,cell_count=n,number_of_masks=end,rectangles=len(boxes),score_numerator=best,score_denominator=sum(mass),private_mask=arg,max_boxes=chosen,chunks=parts,seconds=time.monotonic()-start,workers=12)
            path=args.outdir/('grid_'+'x'.join(map(str,shape))+'.json');path.write_text(json.dumps(record,indent=2)+'\n');print(shape,best,sum(mass),flush=True)
if __name__=='__main__':main()
