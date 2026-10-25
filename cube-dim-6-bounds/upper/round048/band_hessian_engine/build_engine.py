from pathlib import Path
import json,hashlib,tarfile,subprocess,flint
HERE=Path(__file__).resolve().parent;OLD=HERE.parents[1]/'round010/arb_engine';deps=json.loads((OLD/'dependencies.json').read_text());archive=OLD/'headers_min.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest()==deps['minimal_header_archive_sha256']and flint.__FLINT_VERSION__=='3.6.0'
target=HERE/'scratch';target.mkdir(exist_ok=True)
with tarfile.open(archive,'r:gz')as t:t.extractall(target,filter='data')
lib=Path(flint.__file__).resolve().parent/'.dylibs'
for fn,sha in deps['runtime_library_sha256'].items():assert hashlib.sha256((lib/Path(fn).name).read_bytes()).hexdigest()==sha
out=target/'band_hessian';cmd=['/usr/bin/clang++','-O3','-std=c++20','-I',str(target/'include'),str(HERE/'band_hessian.cpp'),str(lib/'libflint.24.0.dylib'),'-Wl,-rpath,'+str(lib),'-o',str(out)];subprocess.run(cmd,check=True)
patch=['/usr/bin/install_name_tool','-change','/DLC/flint/.dylibs/libflint.24.0.dylib','@rpath/libflint.24.0.dylib',str(out)];subprocess.run(patch,check=True)
(HERE/'build_record.json').write_text(json.dumps({'compile':cmd,'local_binary_patch':patch,'binary_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/'band_hessian.cpp').read_bytes()).hexdigest(),'headers_archive':str(archive),'headers_sha256':deps['minimal_header_archive_sha256'],'runtime_dylib_sha256':deps['runtime_library_sha256']},indent=2)+'\n');print(out)
