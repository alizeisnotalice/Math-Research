"""Small integer polytope helpers, no floating geometry."""
from pathlib import Path
from math import gcd
import hashlib,json,gzip,os

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dumps(obj):return json.dumps(obj,sort_keys=True,separators=(',',':')).encode()
def save(path,obj):
    path=Path(path);tmp=Path(str(path)+'.partial');data=dumps(obj)
    if path.suffix=='.gz':data=gzip.compress(data,mtime=0)
    tmp.write_bytes(data);os.replace(tmp,path)
def read(path):
    path=Path(path);data=path.read_bytes()
    if path.suffix=='.gz':data=gzip.decompress(data)
    return json.loads(data)
def ids(mask):
    while mask:
        bit=mask&-mask;yield bit.bit_length()-1;mask-=bit
def mask_of(indices):
    m=0
    for i in indices:m|=1<<i
    return m
def rank(rows):
    a=[list(map(int,r)) for r in rows if any(r)];r=0
    for j in range(5):
        p=next((i for i in range(r,len(a)) if a[i][j]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];pivot=a[r][j]
        for i in range(r+1,len(a)):
            if a[i][j]:
                z=[pivot*x-a[i][j]*y for x,y in zip(a[i],a[r])];g=0
                for x in z:g=gcd(g,abs(x))
                a[i]=[x//g for x in z] if g else z
        r+=1
        if r==5:return r
    return r
def determinant(rows):
    a=[list(map(int,r)) for r in rows];n=len(a);previous=1;sign=1
    for j in range(n-1):
        p=next((i for i in range(j,n) if a[i][j]),None)
        if p is None:return 0
        if p!=j:a[j],a[p]=a[p],a[j];sign=-sign
        pivot=a[j][j]
        for i in range(j+1,n):
            for k in range(j+1,n):
                value=pivot*a[i][k]-a[i][j]*a[j][k]
                assert value%previous==0;a[i][k]=value//previous
            a[i][j]=0
        previous=pivot
    return sign*a[-1][-1]
