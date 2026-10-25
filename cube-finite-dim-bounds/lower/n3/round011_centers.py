"""Generate <=128 deterministic receiver-gap/exposed-face centers, no captures."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,heapq,json,time
import numpy as np
import round011_geometry as geo
ROOT=Path(__file__).resolve().parent

class Picker:
 def __init__(self):self.octants=[[] for _ in range(8)];self.globalheap=[];self.count=0
 def add(self,quality,octant,record):
  self.count+=1
  for heap,cap in [(self.octants[octant],8),(self.globalheap,64)]:
   item=(quality,record)
   if len(heap)<cap:heapq.heappush(heap,item)
   elif quality>heap[0][0]:heapq.heapreplace(heap,item)
 def selected(self):
  ranked=[sorted(h,key=lambda t:t[0],reverse=True) for h in self.octants];out=[];seen=set()
  for q in range(8):
   for octant in range(8):
    if len(ranked[octant])>q:
     record=ranked[octant][q][1];identifier=record['candidate_id']
     if identifier not in seen:out.append(record);seen.add(identifier)
  for quality,record in sorted(self.globalheap,key=lambda t:t[0],reverse=True):
   if len(out)==64:break
   if record['candidate_id'] not in seen:out.append(record);seen.add(record['candidate_id'])
  return out

def main():
 start=time.monotonic();boxes,oldids=geo.old_geometry();mid=[[[geo.midpoint_integer(e) for e in pair] for pair in b] for b in boxes];edges=[]
 for p in geo.P:
  q=p*geo.Q;assert q.denominator==1;edges.append(int(q))
 for b in mid:
  for j in range(3):b[j][0]=max(0,b[j][0]);b[j][1]=min(edges[j],b[j][1])
 assert all(all(b[j][1]>b[j][0] for j in range(3)) for b in mid)
 cuts=[sorted({0,edges[j]}|{b[j][k] for b in mid for k in range(2)}) for j in range(3)];lookup=[{x:i for i,x in enumerate(c)} for c in cuts];idx=np.array([[[lookup[j][b[j][k]] for k in range(2)] for j in range(3)] for b in mid],dtype=np.int32);shape=[len(c)-1 for c in cuts];widths=[[v-u for u,v in zip(c,c[1:])] for c in cuts];eligible=[np.array([w>4 for w in ww],dtype=bool) for ww in widths];octflags=[[int(c[i]+c[i+1]>=edges[j]) for i in range(len(c)-1)] for j,c in enumerate(cuts)]
 def plane(axis,i):
  tangents=[j for j in range(3) if j!=axis];a,b=tangents;chosen=idx[(idx[:,axis,0]<=i)&(idx[:,axis,1]>i)];delta=np.zeros((len(cuts[a]),len(cuts[b])),dtype=np.int32)
  if len(chosen):
   aa,bb=chosen[:,a,0],chosen[:,a,1];cc,dd=chosen[:,b,0],chosen[:,b,1]
   np.add.at(delta,(aa,cc),1);np.add.at(delta,(bb,cc),-1);np.add.at(delta,(aa,dd),-1);np.add.at(delta,(bb,dd),1)
  return delta.cumsum(axis=0,dtype=np.int32).cumsum(axis=1,dtype=np.int32)[:-1,:-1]>0
 empty_picker=Picker();first=plane(0,0);last=plane(0,shape[0]-1);prev=last;cur=first;nxt=plane(0,1)
 for i in range(shape[0]):
  if eligible[0][i]:
   ym=np.roll(cur,1,axis=0);yp=np.roll(cur,-1,axis=0);zm=np.roll(cur,1,axis=1);zp=np.roll(cur,-1,axis=1);axes=(prev|nxt).astype(np.uint8)+(ym|yp).astype(np.uint8)+(zm|zp).astype(np.uint8);faces=prev.astype(np.uint8)+nxt.astype(np.uint8)+ym.astype(np.uint8)+yp.astype(np.uint8)+zm.astype(np.uint8)+zp.astype(np.uint8);valid=(~cur)&(axes>=2)&eligible[1][:,None]&eligible[2][None,:]
   for j,k in np.argwhere(valid):
    j=int(j);k=int(k);degree=int(axes[j,k]);facecount=int(faces[j,k]);vol=widths[0][i]*widths[1][j]*widths[2][k];octant=octflags[0][i]+2*octflags[1][j]+4*octflags[2][k];record={'candidate_id':f'A:{i}:{j}:{k}','type':'empty_gap_cell','cell_indices':[i,j,k],'contact_axes':degree,'covered_neighbor_faces':facecount,'volume_lattice_units':str(vol),'octant':octant};empty_picker.add((degree,facecount,-vol,-i,-j,-k),octant,record)
  prev,cur=cur,nxt
  if i+2<shape[0]:nxt=last if i+2==shape[0]-1 else plane(0,i+2)
  else:nxt=first
 face_picker=Picker()
 for axis in range(3):
  tangents=[j for j in range(3) if j!=axis];a,b=tangents;loevents={};hievents={}
  for bi,bb in enumerate(idx):loevents.setdefault(int(bb[axis,0]),[]).append(bi);hievents.setdefault(int(bb[axis,1]),[]).append(bi)
  left=plane(axis,0)
  for f in range(1,shape[axis]):
   right=plane(axis,f)
   if eligible[axis][f-1] and eligible[axis][f]:
    codes=left.astype(np.int8)-right.astype(np.int8);owners={}
    for j,k in np.argwhere(codes!=0):
     j=int(j);k=int(k);side=int(codes[j,k]);mask=0
     for bi in (hievents.get(f,[]) if side>0 else loevents.get(f,[])):
      bb=idx[bi]
      if bb[a,0]<=j<bb[a,1] and bb[b,0]<=k<bb[b,1]:mask|=1<<bi
     assert mask;owners[(j,k)]=mask
    for (j,k),mask in owners.items():
     if not eligible[a][j] or not eligible[b][k]:continue
     code=int(codes[j,k]);switch=0
     for dj,dk in [(-1,0),(1,0),(0,-1),(0,1)]:
      jj,kk=j+dj,k+dk
      if jj<0 or jj>=shape[a] or kk<0 or kk>=shape[b] or int(codes[jj,kk])!=code or owners.get((jj,kk),0)!=mask:switch+=1
     area=widths[a][j]*widths[b][k];mid2=[0,0,0];mid2[axis]=2*cuts[axis][f];mid2[a]=cuts[a][j]+cuts[a][j+1];mid2[b]=cuts[b][k]+cuts[b][k+1];octant=sum((1<<d)*int(mid2[d]>=edges[d]) for d in range(3));ownerids=[];m=mask
     while m:
      bit=m&-m;bi=bit.bit_length()-1;ownerids.append(oldids[bi]);m^=bit
     record={'candidate_id':f'B:{axis}:{f}:{j}:{k}','type':'exposed_face_patch','normal_axis':axis,'plane_index':f,'tangent_axes':[a,b],'tangent_cell_indices':[j,k],'inside_side':'left' if code>0 else 'right','owner_switch_degree':switch,'owner_old_key_ids':ownerids,'area_lattice_units':str(area),'octant':octant};face_picker.add((switch,len(ownerids),-area,-axis,-f,-j,-k),octant,record)
   left=right
 selected=empty_picker.selected()+face_picker.selected();assert len(selected)<=128;centers=[]
 for q,record in enumerate(selected):
  if record['type']=='empty_gap_cell':
   ijk=record['cell_indices'];point=[F(cuts[d][ijk[d]]+cuts[d][ijk[d]+1],2*geo.Q) for d in range(3)]
  else:
   axis=record['normal_axis'];f=record['plane_index'];a,b=record['tangent_axes'];j,k=record['tangent_cell_indices'];point=[None]*3;point[axis]=F(cuts[axis][f],geo.Q);point[a]=F(cuts[a][j]+cuts[a][j+1],2*geo.Q);point[b]=F(cuts[b][k]+cuts[b][k+1],2*geo.Q)
  assert all(0<x<geo.P[d] for d,x in enumerate(point));record=dict(record);record['center_id']=q;record['point']=[str(x) for x in point];record['normalized_point']=[str(point[d]/geo.P[d]) for d in range(3)];centers.append(record)
 tail=geo.W
 for p in geo.P:tail*=1/p+F(1,6)
 assert tail<1 and min(geo.P)>1
 output={'status':'geometry center bank locked before event capture','seed_input_sha256':hashlib.sha256((ROOT/'round008_certificate_input.json').read_bytes()).hexdigest(),'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'geometry_helper_sha256':hashlib.sha256((ROOT/'round011_geometry.py').read_bytes()).hexdigest(),'root_denominator':geo.D,'coordinate_lattice_denominator':geo.Q,'face_midpoint_error_upper':str(F(1,geo.Q)),'minimum_candidate_cell_width_units_exclusive':4,'old_nonempty_boxes':len(boxes),'axis_coordinate_cut_counts':[len(c) for c in cuts],'eligible_empty_cells':empty_picker.count,'eligible_face_patches':face_picker.count,'chosen_empty_centers':len(empty_picker.selected()),'chosen_face_centers':len(face_picker.selected()),'center_count':len(centers),'tail_radius':3,'exact_uniform_density_tail_bound':str(tail),'tail_bound_below1':True,'centers':centers,'elapsed_seconds':time.monotonic()-start};(ROOT/'round011_centers.json').write_text(json.dumps(output,indent=2)+'\n');print({k:v for k,v in output.items() if k!='centers'});print('FIRST_FIXED_CENTERS',json.dumps(centers[:3]))
if __name__=='__main__':main()
