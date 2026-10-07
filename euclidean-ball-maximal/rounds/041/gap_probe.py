#!/usr/bin/env python3
"""Round 41 exact full-band pair counterexample. Fraction only, no source sampling."""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import argparse
import json
import sys
import time
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'runtime/work/general_input_20261003'))
from stage24_clock_replacement_probe import observer


def overlap(lo, hi, center, radius):
    return max(F(0), min(hi, center + radius) - max(lo, center - radius))


def evaluate(k):
    assert k >= 10
    r = F(2) ** (2-k)
    count = 2 ** (k-8)
    micro = [(-8*r*(i+F(1,2)), F(3,2)*r) for i in range(count)]
    atoms = sorted(micro + [(F(9,32), F(27,64)), (F(10), F(71,128))])
    assert sum(w for y,w in atoms) == 1
    saved, _, radii, _ = observer(atoms, F(1), 'band', 4)
    assert F(saved['entrance_a']) == F(1,8)
    assert F(saved['terminal_b']) == F(1,2)
    cells = [(F(v['lo']),F(v['hi']),v['J'],v['K']) for v in saved['observer_cells']]
    by_layer = defaultdict(list)
    for lo,hi,J,K in cells:
        assert J >= 2
        by_layer[K].append((lo,hi))
    def h(y, layer):
        return sum((overlap(lo,hi,y,radii[layer])/(2*radii[layer])
                    for lo,hi in by_layer[layer]), F(0))
    def window_is_actual(lo,hi,layer):
        length = sum((max(F(0),min(hi,b)-max(lo,a))
                      for a,b,J,K in cells if J == 2 and K == layer),F(0))
        assert length == hi-lo
    window_is_actual(F(3,32),F(1,8),4)
    for y,w in micro:
        window_is_actual(y-r/2,y-F(3,8)*r,k)
        window_is_actual(y+F(3,8)*r,y+r/2,k)
        assert h(y,4) >= F(1,16) and h(y,k) >= F(1,8)
    pair = sum((w*h(y,4)*h(y,k) for y,w in atoms),F(0))
    X = F(saved['eligible_volume'])
    X4 = sum((b-a for a,b in by_layer[4]),F(0))
    Xk = sum((b-a for a,b in by_layer[k]),F(0))
    Q = sum((w*sum((h(y,l) for l in by_layer),F(0))**2 for y,w in atoms),F(0))
    assert pair >= F(3,16384)
    assert 0 < X4 <= X <= F(77,32) and 0 < Xk <= X
    assert Q <= 4*X
    return dict(k=k,gap=k-4,micro_count=count,r=str(r),D=saved['D'],
        original_volume=saved['original_volume'],X=str(X),X4=str(X4),Xk=str(Xk),
        I4k=str(pair),I4k_over_X=str(pair/X),I4k_over_X4=str(pair/X4),
        I4k_over_Xk=str(pair/Xk),Q=str(Q),Q_over_X=str(Q/X),
        min_micro_h4=str(min(h(y,4) for y,w in micro)),
        min_micro_hk=str(min(h(y,k) for y,w in micro)),
        eligible_cell_count=len(cells),
        certified_pair_lower='3/16384',certified_total_X_upper='77/32',
        actual_complete_band_and_original_J_K=True,
        source=[dict(location=str(y),mass=str(w)) for y,w in atoms])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('k',type=int,nargs='*')
    args = parser.parse_args()
    values = args.k or [10,12,14]
    rows = []
    start = time.monotonic()
    for k in values:
        row = evaluate(k)
        rows.append(row)
        print(json.dumps({key:value for key,value in row.items() if key != 'source'},ensure_ascii=False),flush=True)
    result = dict(status='passed',arithmetic='Fraction exact',dimension=1,
        actual_source_weighting=True,stopping_simulated=False,
        uniform_Q_bound_disproved=False,per_pair_gap_decay_disproved=True,
        elapsed_seconds=time.monotonic()-start,cases=rows)
    output = args.output
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('Saved '+str(output),flush=True)

if __name__ == '__main__':
    main()
