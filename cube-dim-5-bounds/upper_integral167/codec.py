"""Exact binary projection only. No numerical kernel or envelope evaluation."""
from fractions import Fraction as F
import hashlib,json,struct
MAGIC=b'I167H001'
def dyad(s):
    x=F(s);d=x.denominator;assert d>0 and d&(d-1)==0
    return x.numerator,d.bit_length()-1
def norm(s):return str(F(s))
def projection(point,times,values):
    x={'point':[norm(s) for s in point],'endpoints':[{'t':t,'H':{'lower':str(l),'upper':str(u)}} for t,(l,u) in zip(times,values)]}
    return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def uint(n):
    assert n>=0;b=n.to_bytes(max(1,(n.bit_length()+7)//8),'big');assert len(b)<65536
    return struct.pack('>H',len(b))+b
def take(f,n):
    b=f.read(n);assert len(b)==n,'truncated167';return b
def readuint(f):
    n=struct.unpack('>H',take(f,2))[0];assert 0<n<65536
    b=take(f,n);assert n==1 or b[0]!=0;return int.from_bytes(b,'big')
def encode(row,raw,times):
    point=raw.get('point',raw.get('anchor'));assert point==row['point']
    ends=raw['endpoints'];assert len(ends)==73 and [e['t'] for e in ends]==times
    b=bytearray(struct.pack('>I',row['id']));vals=[]
    for e in ends:
        l,u=map(dyad,[e['H']['lower'],e['H']['upper']]);E=max(l[1],u[1]);L=l[0]<<(E-l[1]);U=u[0]<<(E-u[1]);assert 0<L<=U and E<2147483648
        b.extend(struct.pack('>i',E));b.extend(uint(L));b.extend(uint(U-L));vals.append((F(L,1<<E),F(U,1<<E)))
    psha=projection(point,times,vals);b.extend(bytes.fromhex(psha));b.extend(bytes.fromhex(row['packet_sha256']))
    original=raw.get('original_raw_sha256',row['packet_sha256']);b.extend(bytes.fromhex(original))
    return bytes(b),psha,original
def decode(f,row,times):
    assert struct.unpack('>I',take(f,4))[0]==row['id'];vals=[]
    for _ in times:
        E=struct.unpack('>i',take(f,4))[0];assert 0<=E<=100000
        L,W=readuint(f),readuint(f);assert L>0;vals.append((F(L,1<<E),F(L+W,1<<E)))
    psha,rawsha,original=(take(f,32).hex() for _ in range(3))
    assert rawsha==row['packet_sha256'] and psha==projection(row['point'],times,vals)
    return vals,psha,original
