#!/usr/bin/env python3
"""Restore the two exact123 incremental states; byte checks only, no math rerun."""
from pathlib import Path
import argparse,gzip,hashlib,json
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,default=Path(__file__).resolve().parent/'packet');ap.add_argument('--out-dir',type=Path);args=ap.parse_args();p=args.packet;m=json.load(open(p/'manifest.json'));decoded={}
 for name,row in m['compressed_files'].items():
  data=(p/name).read_bytes();assert sha(data)==row['gzip_sha256'];raw=gzip.decompress(data);assert sha(raw)==row['uncompressed_sha256'];decoded[row['role']]=json.loads(raw)
 final=decoded['state_batch0002'];tasks1=decoded['tasks_batch0001'];tasks2=decoded['tasks_batch0002'];points1=decoded['points_batch0001'];points2=decoded['points_batch0002']
 assert len(tasks1)==len(points1['records'])==4096 and len(tasks2)==len(points2['records'])==2141
 for task,row in zip(tasks1+tasks2,points1['records']+points2['records']):
  assert task['point']==row['point'] and task['node']==row['node'] and row['record']==final['new_cache_by_epsilon']['1/128'][task['pair_cache_key'][1]]
 nodes={tuple(t['node']) for t in tasks1};header=json.load(open(p/'state_header_batch0001.json'));first=dict(header)
 first['new_cache_by_epsilon']={'1/128':{t['pair_cache_key'][1]:final['new_cache_by_epsilon']['1/128'][t['pair_cache_key'][1]] for t in tasks1}}
 first['new_passed']=[r for r in final['new_passed'] if tuple(r['node']) in nodes];first['new_splits']=[r for r in final['new_splits'] if tuple(r['node']) in nodes]
 first['queue']=[t['node'] for t in tasks2]+[n for s in first['new_splits'] for n in s['children']]
 enc=lambda x:(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode();raw1=enc(first);raw2=enc(final)
 assert sha(raw1)==m['state_batch0001_uncompressed_sha256'] and sha(raw2)==m['state_batch0002_uncompressed_sha256']
 if args.out_dir:
  args.out_dir.mkdir(exist_ok=False,parents=True);(args.out_dir/'state_batch0001.json').write_bytes(raw1);(args.out_dir/'state_batch0002.json').write_bytes(raw2)
 print(json.dumps({'status':'exact_two_increment_states_and6237_raw_records_restored','state_batch0001_sha256':sha(raw1),'state_batch0002_sha256':sha(raw2),'new_H_gradient_log_corner_calls':0,'source_seed_copied':False}))
if __name__=='__main__':main()
