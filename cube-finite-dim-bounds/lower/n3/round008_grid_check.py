"""Independent integer-grid enclosure of the new three-orbit candidate."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
from round006_grid_check import build,union_integer,D,G
ROOT=Path(__file__).resolve().parent

def main():
 start=time.monotonic();fp=ROOT/'round008_certificate_input.json';inp=json.loads(fp.read_text());inner,outer,W=build(inp);vi,ci=union_integer(inner);vo,co=union_integer(outer);lo=F(vi,G**3)/W;hi=F(vo,G**3)/W;target=F(json.loads((ROOT/'round008_target.json').read_text())['ratio_interval'][1]);out={'status':'independent three-orbit integer-grid enclosure','input_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'root_denominator':D,'coordinate_grid_denominator':G,'ratio_interval':[str(lo),str(hi)],'decimal_display_only':[float(lo),float(hi)],'target_exact_upper':str(target),'lower_above_full_015_upper':lo>target,'rigorous_gain_lower':str(lo-target),'mass_per_index_triple':str(W),'inner_boxes':len(inner),'outer_boxes':len(outer),'inner_cuts':ci,'outer_cuts':co,'elapsed_seconds':time.monotonic()-start};(ROOT/'round008_grid_result.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
if __name__=='__main__':main()
