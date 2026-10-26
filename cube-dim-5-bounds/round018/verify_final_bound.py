"""Minimal portable exact C5 comparison check; no old H/R/log/corner replay."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parents[1]
def read(n):return json.loads((P/n).read_bytes())
def main():
    a=read('round018/primary_integral156_acceptance.json')
    b=read('round018/primary_curvature155_acceptance.json')
    c=read('round018/primary_integral150_acceptance.json')
    old=read('round017/primary_global_acceptance.json')
    audit=read('audit_integral161/run002/results.json')
    assert audit['all_source_hashes_unchanged'] and audit['worker_exitcodes']==[0]*15
    assert audit['qualified_total_mean_checks']==97907 and audit['qualified_total_band_products']==7147211
    assert Q(a['new_UH'])==Q(audit['new_UH'])
    m=Q(b['minimum_margin110']);term=Q(1);poly=term
    for k in range(1,9):term=term*m/k;poly+=term
    g=Q(a['qualified_g']);exact=Q(11,100000)*poly;N=1<<192
    assert 0<=exact-g<Q(1,N) and (g*N).denominator==1 and g>Q(old['qualified_g'])
    I=Q(a['complete_I_upper']);assert I==Q(a['new_UH'])+Q(a['same_complete_J139_upper'])<Q(c['complete_I_upper'])
    fourth=3125*I**4/g**5;assert fourth==Q(a['best_upper']['fourth_power'])
    lo,hi=map(Q,a['best_upper']['interval']);assert 0<lo<hi and lo**4<=fourth<=hi**4
    assert hi<=Q(a['simple_rational_upper'])<Q(32)
    assert hashlib.sha256((P/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest()=='f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2'
    print('PASS:complete g192/I156/C5 exact certificate, source lower unchanged')
if __name__=='__main__':main()
