"""Python ctypes wrapper for the numerical C++ full4D union kernel."""
from pathlib import Path
import ctypes,subprocess,time
import numpy as np
P=Path(__file__).resolve().parent
SOURCE=P/'round007_union_kernel.cpp';LIBRARY=P/'round007_union_kernel.dylib'
if not LIBRARY.exists() or LIBRARY.stat().st_mtime<SOURCE.stat().st_mtime:
    subprocess.run(['/usr/bin/c++','-O3','-std=c++17','-dynamiclib',str(SOURCE),'-o',str(LIBRARY)],check=True)
_lib=ctypes.CDLL(str(LIBRARY));_call=_lib.union4_volume;_call.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.POINTER(ctypes.c_longlong)];_call.restype=ctypes.c_double

def volume(boxes):
    st=time.monotonic();b=np.ascontiguousarray(boxes,dtype=np.float64).reshape(-1,8);states=ctypes.c_longlong()
    value=_call(b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),len(b),ctypes.byref(states))
    return value,dict(backend='C++ endpoint event sweep',cache_states=states.value,seconds=time.monotonic()-st)
