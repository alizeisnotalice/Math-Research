from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from PIL import Image,ImageDraw
import subprocess,os,json,time,hashlib,resource,shutil
out=Path(__file__).resolve().parent;pdf=out/'input.pdf';renderer=shutil.which('pdftoppm')
expected='5f84aaa0a3d18fad3e9a86eea5a275bf9c5222a13d8171f2d10a7e47d6ff98a1'
def sha():return hashlib.sha256(pdf.read_bytes()).hexdigest()
def render(page):
 before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.monotonic()
 cmd=[renderer,'-f',str(page),'-l',str(page),'-singlefile','-r','100','-png',str(pdf),str(out/f'page-{page:03d}')]
 p=subprocess.run(cmd,capture_output=True,text=True)
 if p.returncode:raise RuntimeError((page,p.returncode,p.stderr))
 after=resource.getrusage(resource.RUSAGE_CHILDREN);image=out/f'page-{page:03d}.png'
 with Image.open(image) as im:size=list(im.size)
 return {'page':page,'worker_pid':os.getpid(),'seconds':time.monotonic()-t,'cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'size':size,'bytes':image.stat().st_size}
if __name__=='__main__':
 t=time.monotonic();assert sha()==expected
 info=subprocess.check_output(['pdfinfo',str(pdf)],text=True);pages=int(next(s.split(':',1)[1] for s in info.splitlines() if s.startswith('Pages:')));assert pages==73
 with ProcessPoolExecutor(max_workers=15) as pool:rows=list(pool.map(render,range(1,pages+1)))
 for group in range((pages+15)//16):
  sheet=Image.new('RGB',(1200,1720),'#dddddd');draw=ImageDraw.Draw(sheet)
  for j,page in enumerate(range(group*16+1,min(pages+1,group*16+17))):
   with Image.open(out/f'page-{page:03d}.png') as im:
    im.thumbnail((280,395));x=(j%4)*300+(300-im.width)//2;y=(j//4)*430+23;sheet.paste(im,(x,y));draw.text((x,y-18),str(page),fill='black')
  sheet.save(out/f'contact-{group+1}.png')
 assert sha()==expected
 report={'input_sha256':expected,'pages':pages,'dpi':100,'geometry_sizes':sorted({tuple(r['size']) for r in rows}),'requested_workers':15,'actual_workers':len({r['worker_pid'] for r in rows}),'wall_seconds':time.monotonic()-t,'sum_cpu_seconds':sum(r['cpu_seconds'] for r in rows),'rows':rows}
 (out/'render-result.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
