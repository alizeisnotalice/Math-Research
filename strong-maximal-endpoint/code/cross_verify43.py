"""Independent integer verification of the finite real-rectangle witnesses."""
import hashlib,itertools,json,sys
from fractions import Fraction
from pathlib import Path
edges=(-8,-4,-1,0,1,2,3,4,5,8,12)
pairs=list(itertools.combinations(edges,2))
rectangles=[(xl,xr,yl,yr) for (xl,xr),(yl,yr) in itertools.product(pairs,pairs)]
receiver=[(i,j) for i in range(10) for j in range(10)]
def contained(rect,cell):
    xl,xr,yl,yr=rect;i,j=cell
    return xl<=edges[i] and edges[i+1]<=xr and yl<=edges[j] and edges[j+1]<=yr
def area(cells):
    return sum((edges[i+1]-edges[i])*(edges[j+1]-edges[j]) for i,j in cells)
def source_count(rect,cells):
    xl,xr,yl,yr=rect
    return sum(xl<=i-1 and i<=xr and yl<=j-1 and j<=yr for i,j in cells)
def witness(cells,amplitude):
    coverage={};margins=[];ngood=0
    for rect in rectangles:
        ar=(rect[1]-rect[0])*(rect[3]-rect[2]);mass=amplitude*source_count(rect,cells)
        if 1000*mass<=1001*ar:continue
        ngood+=1;margins.append(Fraction(mass,ar)-Fraction(1001,1000))
        for cell in receiver:
            if cell not in coverage and contained(rect,cell):coverage[cell]=(rect,mass,ar)
    return coverage,ngood,min(margins)
if __name__=='__main__':
    repo=Path(__file__).resolve().parent.parent
    source=Path(sys.argv[1]) if len(sys.argv)>1 else repo/'rounds/r043/cross_rectangles43.json'
    output=Path(sys.argv[2]) if len(sys.argv)>2 else Path('/tmp/cross_certificate43.json')
    data=json.loads(source.read_text());summaries=[];certificates=[]
    assert len(rectangles)==3025 and len(data['rows'])==20
    for row in data['rows']:
        components=[{tuple(x) for x in part} for part in row['components']]
        cells=set.union(*components);assert sum(map(len,components))==len(cells)
        amp=row['amplitude'];full,nrect,margin=witness(cells,amp)
        comp_results=[witness(c,amp) for c in components]
        comp=set().union(*(set(v[0]) for v in comp_results));cross=set(full)-comp
        assert set(full)>=comp
        assert (area(full),area(comp),area(cross),nrect)==(row['full_witness_union_area'],row['component_witness_union_area'],row['cross_component_witness_area'],row['full_good_rectangle_count'])
        assert [v[1] for v in comp_results]==row['component_good_rectangle_counts']
        assert margin==Fraction(row['least_full_rectangle_margin_over_tau'])
        assert len(cells)==row['input_cell_count']
        summaries.append({'pattern':row['pattern'],'amplitude':amp,'full_area':area(full),'component_area':area(comp),'cross_area':area(cross),'denominator_multiplier':len(cells)*amp,'denominator_log_argument':amp,'least_margin_over_tau':str(margin)})
        if amp==8 and row['pattern'] in ('horizontal_two_1x2','offset_two_1x3'):
            witnesses=[]
            for cell,(rect,mass,ar) in sorted(full.items()):
                i,j=cell;assert 1000*mass>1001*ar and contained(rect,cell)
                witnesses.append({'cell_bounds':[edges[i],edges[i+1],edges[j],edges[j+1]],'rectangle_bounds':list(rect),'source_integral':mass,'rectangle_area':ar,'exact_average':str(Fraction(mass,ar))})
            certificates.append({'pattern':row['pattern'],'amplitude':amp,'source_components':row['components'],'full_area':area(full),'witnesses':witnesses})
    result={'status':'PASS: independent integer enumeration; explicit strict continuous witnesses for two candidates','input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rectangle_count':len(rectangles),'verified_rows':len(summaries),'summaries':summaries,'certificates':certificates}
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('PASS: 20 exact unions and denominators; explicit certificates',[(c['pattern'],c['full_area'],len(c['witnesses'])) for c in certificates])
