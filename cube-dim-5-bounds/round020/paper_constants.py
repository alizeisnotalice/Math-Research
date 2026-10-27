"""Readable outward scalar certificate for current audited175 upper;not new bound."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,os,time
P=Path(__file__).resolve().parents[1];r=json.loads((P/'round020/paper_constants_registration.json').read_bytes());f=P/r['source'];assert hashlib.sha256(f.read_bytes()).hexdigest()==r['source_sha256'];assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==r['code_sha256'];out=P/'round020/paper_constants.json';assert not out.exists();t=time.monotonic();s=json.loads(f.read_bytes());g,I=Q(s['qualified_g']),Q(s['complete_I_upper']);N=10**20;G=g.numerator*N//g.denominator;H=(I.numerator*N+I.denominator-1)//I.denominator;K=28920678;gap=K**4*G**5-3125*H**4*10**44;assert Q(G,N)<=g and Q(H,N)>=I and gap>0
z={'status':'PASS_readable_outward175_paper_certificate','pid':os.getpid(),'seconds':time.monotonic()-t,'N':N,'G':G,'H':H,'K':K,'strict_integer_gap':str(gap),'g_lower':str(Q(G,N)),'I_upper':str(Q(H,N)),'strict_portable_upper':str(Q(K,10**6)),'source_sha256':r['source_sha256'],'new_stronger_bound':False};out.write_text(json.dumps(z,indent=2)+'\n');print(json.dumps(z))
