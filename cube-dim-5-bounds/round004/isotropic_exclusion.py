from fractions import Fraction as F
from pathlib import Path
import hashlib,json
pi_lower=F(333,106)
lower_sixth=pi_lower**12*5**25/(60**6*2**10)
diff=lower_sixth-32**6
assert diff>0
out={'status':'exact_exclusion_of_isotropic_positive_part_cost_not_C5_lower_bound','alpha_range':'0<alpha<=2','face_shape_integral':'pi^2/60','minimum_cost':'(5/2)^(5/3)*5^(5/2)*pi^2/60','pi_lower':'333/106','minimum_cost_sixth_power_strict_lower':str(lower_sixth),'difference_from_32_sixth_power':str(diff),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
