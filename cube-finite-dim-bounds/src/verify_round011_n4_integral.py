"""Portable entry point for the locked CI4B certificate (no optimizer).

Only input locations are relocated. The locked checker, arithmetic, certificate,
and source hashes remain unchanged. --inputs-only does not certify the integral.
"""
from pathlib import Path
import argparse
import hashlib
import sys

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / 'upper/continuous_integral_n4'
LOCKS = {
    'ci4b_affine.py': '5c54cac11da5989d65cf35f1cff03729be47b5d10d4ceb9181a1ec51789e9e34',
    'ci4b_certificate.py': 'babf044645845ad0a6bc5f372e9a007f2348601fb25565182c7bc6f736dbaead',
    'ci4b_certificate.json': '4dbb544ba34762b239adea8e4dc152987ce33c15ef855f459efb63bffa1281b0',
}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs-only', action='store_true')
parser.add_argument('--save', help='Optional output result path; default stdout only')
args = parser.parse_args()
for name, expected in LOCKS.items():
    assert hashlib.sha256((DIRECTORY / name).read_bytes()).hexdigest() == expected, name
sys.path.insert(0, str(DIRECTORY))
import ci4b_certificate as verifier
verifier.engine.SOURCE = ROOT / 'upper/u11_minimum_global_results.json'
verifier.engine.INTERFACE = ROOT / 'upper/u11_minimum_global_interface.py'
if args.inputs_only:
    assert len(verifier.engine.models()) == 3
    assert len(verifier.cells()) == 2440
    print('PASS portable input hashes, models, cell geometry; integral NOT rechecked')
else:
    verifier.check(str(DIRECTORY / 'ci4b_certificate.json'), args.save)
