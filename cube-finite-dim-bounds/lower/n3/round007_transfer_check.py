"""Reuse the verified three-orbit finite-block proof on the new exact input."""
import sys
from round006_transfer_check import main
if __name__=='__main__':
 sys.argv=[sys.argv[0],'--input','round007_certificate_input.json','--grid','round007_grid_result.json','--output','round007_finite_result.json']
 main()
