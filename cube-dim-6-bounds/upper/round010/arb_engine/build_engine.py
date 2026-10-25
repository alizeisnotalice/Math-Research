"""Build against the existing python-flint3.6.0 wheel; no installation or PATH edit."""
from pathlib import Path
import flint,hashlib,json,subprocess,tarfile

HERE=Path(__file__).resolve().parent
deps=json.loads((HERE/'dependencies.json').read_text())
assert flint.__FLINT_VERSION__=='3.6.0'
archive=HERE/'headers_min.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest()==deps['minimal_header_archive_sha256']
target=HERE/'scratch'
target.mkdir(exist_ok=True)
with tarfile.open(archive,'r:gz')as t:t.extractall(target,filter='data')
for name,digest in deps['minimal_header_members'].items():
    assert hashlib.sha256((target/'include'/name).read_bytes()).hexdigest()==digest
lib=Path(flint.__file__).resolve().parent/'.dylibs'
for path,digest in deps['runtime_library_sha256'].items():
    candidate=lib/Path(path).name
    assert hashlib.sha256(candidate.read_bytes()).hexdigest()==digest
output=target/'arb_kernel'
cmd=['/usr/bin/clang++','-O3','-std=c++20','-I',str(target/'include'),str(HERE/'arb_kernel.cpp'),
     str(lib/'libflint.24.0.dylib'),'-Wl,-rpath,'+str(lib),'-o',str(output)]
subprocess.run(cmd,check=True)
patch=['/usr/bin/install_name_tool','-change','/DLC/flint/.dylibs/libflint.24.0.dylib',
       '@rpath/libflint.24.0.dylib',str(output)]
subprocess.run(patch,check=True)
(HERE/'build_record.json').write_text(json.dumps({'compile_command':cmd,'local_binary_install_name_command':patch,
    'binary_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'kernel_calls':0,
    'headers_archive_sha256':deps['minimal_header_archive_sha256']},indent=2)+'\n')
print(output)
