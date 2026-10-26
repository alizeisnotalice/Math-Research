"""Exact unchanged-continuation comparison, retaining full algebraic sources."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import os
H=Path(__file__).resolve().parent
P=H.parent


def main():
    ledger=json.loads((P/'results/bounds_ledger.json').read_text())
    source=P/ledger['lower']['certificate']
    assert hashlib.sha256(source.read_bytes()).hexdigest()=='f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2'
    Llo,Lhi=map(F,ledger['lower']['exact_interval'])
    intervals={}
    for name,old in [('checked_public',ledger['public_best_verified']['lower']['interval']),
                     ('opening_internal',ledger['opening_thresholds']['lower']['exact_interval'])]:
        lo,hi=map(F,old)
        assert Llo-hi>0
        intervals[name]=list(map(str,[Llo-hi,Lhi-lo]))
    out={'status':'formal_C5_bounds_unchanged_in_round013',
         'full_lower_certificate':ledger['lower']['certificate'],
         'full_lower_certificate_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'same_full_algebraic_lower_as_continuation':True,
         'lower_exact_difference_from_continuation':'0',
         'lower_minus_old_bounds_outward_intervals':intervals,
         'current_upper':'32','upper_minus_checked_public':'0',
         'upper_minus_opening_internal':'0','upper_minus_continuation':'0',
         'checked_public_scope':'best among sources actually checked; not global record claim',
         'internal_not_publicly_published_claim':True,
         'new_C5_lower':False,'new_C5_upper':False,
         'dual_breakthrough_time':None,'post_dual_rounds':0,
         'pid':os.getpid(),'new_baseline_or_union_calls':0}
    path=H/'formal_comparison.json';assert not path.exists()
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','lower_exact_difference_from_continuation','current_upper','pid']}))


if __name__=='__main__':main()
