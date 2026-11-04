import ctypes,time
lib=ctypes.CDLL('/usr/lib/libSystem.B.dylib')
lib.mach_host_self.restype=ctypes.c_uint
lib.mach_task_self.restype=ctypes.c_uint
lib.host_processor_info.argtypes=[ctypes.c_uint,ctypes.c_int,ctypes.POINTER(ctypes.c_uint),ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),ctypes.POINTER(ctypes.c_uint)]
lib.host_processor_info.restype=ctypes.c_int
lib.vm_deallocate.argtypes=[ctypes.c_uint,ctypes.c_uint64,ctypes.c_uint64]
def read():
 n=ctypes.c_uint();count=ctypes.c_uint();p=ctypes.POINTER(ctypes.c_int)()
 rc=lib.host_processor_info(lib.mach_host_self(),2,ctypes.byref(n),ctypes.byref(p),ctypes.byref(count))
 if rc: raise RuntimeError('host_processor_info status '+str(rc))
 data=[p[i] for i in range(count.value)]
 lib.vm_deallocate(lib.mach_task_self(),ctypes.cast(p,ctypes.c_void_p).value,count.value*ctypes.sizeof(ctypes.c_int))
 return [data[i*4:(i+1)*4] for i in range(n.value)]
a=read();time.sleep(2);b=read()
for i,(x,y) in enumerate(zip(a,b)):
 d=[(v-u)%(2**32) for u,v in zip(x,y)];total=sum(d)
 print('CPU',i,'busy_percent',round(100*(total-d[2])/total,1) if total else 'no ticks','ticks',d)
