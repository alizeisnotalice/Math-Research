"""Private framed IPC for this run's own authorized local children only."""
import hashlib,pickle,struct
MAX_FRAME_BYTES=256*1024*1024
def send_frame(stream,obj):
 data=pickle.dumps(obj,protocol=5)
 if len(data)>MAX_FRAME_BYTES:raise RuntimeError('finite IPC frame cap exceeded')
 stream.write(struct.pack('!Q',len(data)));stream.write(data);stream.flush()
 return hashlib.sha256(data).hexdigest()
def _exact(stream,n):
 out=bytearray()
 while len(out)<n:
  part=stream.read(n-len(out))
  if not part:raise EOFError('truncated owned-child frame')
  out.extend(part)
 return bytes(out)
def read_frame(stream):
 n=struct.unpack('!Q',_exact(stream,8))[0]
 if n>MAX_FRAME_BYTES:raise RuntimeError('finite IPC frame cap exceeded')
 data=_exact(stream,n)
 # Never read from external/untrusted pipes or arbitrary pickle files.
 return pickle.loads(data),hashlib.sha256(data).hexdigest()
