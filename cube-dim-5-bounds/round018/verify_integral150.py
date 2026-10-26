"""Portable exact 150 finite arithmetic replay; no H/R/old-power evaluation.

Uses original compressed input bytes and frozen qualified scalar intervals.
This is a reader entry point; its execution is separate from the historical
actual producer69275, independent73479 and primary76240 receipts.
"""
from fractions import Fraction as Q
from pathlib import Path
import base64, gzip, hashlib, json

P=Path(__file__).resolve().parents[1]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    b=p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def ceilq(q): return -((-q.numerator)//q.denominator)

def main():
    S=P/'upper_integral150'; pre=read(S/'preflight_results.json')
    for name in ('factor_catalog.json','ratio_catalog.json','geometries.json.gz'):
        assert sha(P/'upper_path112'/name)==pre['source_locks'][name]['sha256']
    old=read(P/'upper_path112/factor_catalog.json')
    ratios=read(P/'upper_path112/ratio_catalog.json')
    geo=read(P/'upper_path112/geometries.json.gz')['g002']['profile']
    h,mu=Q(geo['h']),Q(geo['mu_infinity']); times=pre['time_endpoints']
    expected=read(S/'actual/factors.json'); factors=[]
    def ends(t):
        if t is None:return Q(1),Q(1)
        a=Q(t)
        return Q(ratios[str(h/a)]['upper']),1+mu/a+Q(25,24)*Q(old['g002'][t]['J_upper'])/a**2
    for j,t in enumerate(times):
        bt=times[j+1] if j+1<len(times) else None
        za,zb=1/Q(t),Q(0) if bt is None else 1/Q(bt)
        # Exact z^3 and z^4 moments, both endpoints including infinity.
        ez=Q(4,5)*(za**5-zb**5)/(za**4-zb**4)
        theta=(ez-zb)/(za-zb);assert 0<theta<1
        ra,ta=ends(t);rb,tb=ends(bt)
        cr=theta*ra+(1-theta)*min(rb,ra)
        ct=theta*ta+(1-theta)*min(tb,ta)
        factor=Q(ceilq(min(Q(old['g002'][t]['factor_upper']),cr,ct)*(1<<192)),1<<192)
        e=expected[j]
        assert (Q(e['theta']),Q(e['range_time_upper']),Q(e['taylor_time_upper']),Q(e['new_upper']))==(theta,cr,ct,factor)
        factors.append(factor)
    receipt=read(P/'round018/input150/receipt.json')
    packets_path=P/'round018/input150/original_packets64.jsonl.gz'
    assert sha(packets_path)==receipt['sha256']
    payload=gzip.decompress(packets_path.read_bytes())
    assert hashlib.sha256(payload).hexdigest()==receipt['uncompressed_sha256']
    packets={}
    for line in payload.splitlines():
        row=json.loads(line);b=base64.b64decode(row['gzip_base64'])
        assert len(b)==row['original_bytes'] and hashlib.sha256(b).hexdigest()==row['original_sha256']
        packets[row['original_sha256']]=json.loads(gzip.decompress(b))
    selected=read(S/'selected64.json');results=read(S/'actual/nodes.json')
    assert len(selected)==len(results)==len(packets)==64
    gain_total=0;N=1<<160
    for item,result in zip(selected,results):
        row=item['source_node'];cert=item['source112_certificate']
        assert row['level']==5 and row['geometry']=='g002' and row['old_U_i']>=1<<144
        assert cert['new_selfU']==item['saved_selfU']<cert['old_selfU']
        raw=packets[row['packet_sha256']];endpoints=raw['endpoints']
        assert raw['point']==row['point'] and [x['t'] for x in endpoints]==times
        total=Q(0)
        for j,e in enumerate(endpoints):
            nxt=endpoints[j+1]['H'] if j+1<len(endpoints) else {'lower':'0','upper':'0'}
            upper=Q(e['H']['upper'])-Q(nxt['lower'])
            lower=Q(e['H']['lower'])-Q(nxt['upper']);assert 0<lower<=upper
            total+=upper*factors[j]
        newmean=min(cert['new_M_i'],ceilq(N*total))
        gain=row['weight']*(cert['new_M_i']-newmean)//(16*(1<<(4*row['level'])))
        assert result['key']==item['key'] and result['new_mean_i']==newmean and result['gain_units']==gain
        assert result['new_selfU']==item['saved_selfU']-gain>=item['saved_selfL']
        gain_total+=gain
    prod=read(S/'actual/results.json');assert Q(gain_total,N)==Q(prod['exact_UH_gain'])>0
    c=read(P/'round018/primary_integral150_acceptance.json')
    assert Q(c['complete_I_upper'])==Q(read(P/'round017/primary_global_acceptance.json')['complete_I_upper'])-Q(gain_total,N)
    fourth=3125*Q(c['complete_I_upper'])**4/Q(c['qualified_g'])**5
    lo,hi=map(Q,c['best_upper']['interval'])
    assert lo**4<=fourth<=hi**4 and hi<=Q(c['simple_rational_upper'])<Q(32)
    print('PASS:73 factors,64 roots,4672 cached products,complete I and C5 comparison')

if __name__=='__main__':main()
