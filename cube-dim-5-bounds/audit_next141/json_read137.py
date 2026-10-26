"""Read-only mmap JSON selections; never imports previous numerical checkers."""
import json
import mmap
from pathlib import Path
import re

tokens=re.compile(rb'"(?:[^"\\]|\\.)*"|[{}\[\]]')
quoted=re.compile(rb'"(?:[^"\\]|\\.)*"')


class Document:
    def __init__(self,path):
        self.handle=Path(path).open('rb')
        self.bytes=mmap.mmap(self.handle.fileno(),0,access=mmap.ACCESS_READ)

    def skip(self,i):
        while i<len(self.bytes) and self.bytes[i] in b' \t\r\n':i+=1
        return i

    def end(self,start):
        if self.bytes[start]==34:return quoted.match(self.bytes,start).end()
        if self.bytes[start] in (91,123):
            balance=0
            for m in tokens.finditer(self.bytes,start):
                c=self.bytes[m.start()]
                if c in (91,123):balance+=1
                elif c in (93,125):
                    balance-=1
                    if balance==0:return m.end()
            raise AssertionError('Incomplete JSON')
        i=start
        while self.bytes[i] not in b',}] \r\n\t':i+=1
        return i

    def object_spans(self,start=0):
        i=self.skip(start);assert self.bytes[i]==123;i=self.skip(i+1)
        while self.bytes[i]!=125:
            q=quoted.match(self.bytes,i);assert q
            k=json.loads(self.bytes[i:q.end()]);i=self.skip(q.end());assert self.bytes[i]==58
            a=self.skip(i+1);b=self.end(a);yield k,(a,b)
            i=self.skip(b)
            if self.bytes[i]==44:i=self.skip(i+1)
            else:assert self.bytes[i]==125

    def load_span(self,span):return json.loads(self.bytes[span[0]:span[1]])

    def array(self,span):
        i=self.skip(span[0]);assert self.bytes[i]==91;i=self.skip(i+1)
        while self.bytes[i]!=93:
            e=self.end(i);yield json.loads(self.bytes[i:e]);i=self.skip(e)
            if self.bytes[i]==44:i=self.skip(i+1)
            else:assert self.bytes[i]==93

    def __enter__(self):return self
    def __exit__(self,*args):self.bytes.close();self.handle.close()
