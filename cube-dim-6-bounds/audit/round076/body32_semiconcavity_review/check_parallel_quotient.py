#!/usr/bin/env python3
"""One exact, fixed 3D quotient fixture. No 6D or body32 kernel calls."""
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import platform
import sys


class QI:
    """Gaussian rationals, used solely for two simple residues."""
    def __init__(self, re=0, im=0):
        self.re, self.im = F(re), F(im)

    @staticmethod
    def cast(z):
        return z if isinstance(z, QI) else QI(z)

    def __add__(self, z):
        z = self.cast(z)
        return QI(self.re + z.re, self.im + z.im)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.re, -self.im)

    def __sub__(self, z):
        return self + -self.cast(z)

    def __rsub__(self, z):
        return self.cast(z) + -self

    def __mul__(self, z):
        z = self.cast(z)
        return QI(self.re*z.re-self.im*z.im, self.re*z.im+self.im*z.re)

    __rmul__ = __mul__

    def __truediv__(self, z):
        z = self.cast(z)
        d = z.re*z.re+z.im*z.im
        return self*QI(z.re/d, -z.im/d)

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = QI(1)
        for _ in range(n):
            out = out*self
        return out


def det3(columns):
    a, b, c = columns
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            - b[0]*(a[1]*c[2]-a[2]*c[1])
            + c[0]*(a[1]*b[2]-a[2]*b[1]))


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def main():
    root = Path(__file__).resolve().parent
    registration = json.loads((root/'registration.json').read_text())
    cfg = registration['configurations'][0]
    require(len(registration['configurations']) == 1, 'Exactly one registered configuration')
    B = [tuple(F(x) for x in b) for b in cfg['columns']]
    t = F(cfg['t'])
    s, y2, y3 = map(F, cfg['basis_coordinates'])
    basis = [B[j] for j in (0, 2, 3)]
    det = abs(det3(basis))
    require(det == 3, 'Fixed slanted basis determinant')
    require(all(B[1][i] == 2*B[0][i] for i in range(3)), 'Parallel columns')
    require(all(det3([B[j] for j in J]) == 0
                for J in itertools.combinations(range(4), 3) if 0 in J and 1 in J),
            'Parallel columns cannot be in one quotient basis')
    basis_dets = {str(j+1): str(abs(det3([B[k] for k in J])))
                  for j, J in [(0, (0, 2, 3)), (1, (1, 2, 3))]}

    # Fixed coarea coordinates: w1=s-2w2, w3=y2, w4=y3.
    # c_t(w2)c_t(s-2w2)=pi^-2*t^2/D(w2).
    # Poles in the upper half plane are it and (s+it)/2, both simple.
    # Degree <=2 weights give O(|w2|^-2) at infinity, so the arc vanishes.
    poles = [QI(0, t), QI(s/2, t/2)]
    tail_factor = (t/(t*t+y2*y2))*(t/(t*t+y3*y3))/det
    c3_coeff = 3*t/(9*t*t+s*s)  # pi*c_(3t)(s)
    results = []
    for j in (0, 1):
        for k in (0, 1, 2):
            residues = QI()
            for w in poles:
                A, D = t*t+(s-2*w)**2, t*t+w**2
                derivative = -4*(s-2*w)*D + A*2*w
                weight = (s-2*w if j == 0 else w)**k
                residues = residues + t*t*weight/derivative
            raw = QI(0, 2)*residues  # pi times one-dimensional moment
            require(raw.im == 0, 'Residue sum must be exactly real')
            expected = (c3_coeff if k == 0 else
                        s*c3_coeff/3 if k == 1 else
                        t-t*t*c3_coeff if j == 0 else
                        t/2-t*t*c3_coeff)
            require(raw.re == expected, 'Exact parallel conditional moment identity')
            physical_coeff = raw.re*tail_factor  # actual 3D moment = coefficient/pi^3
            row = {'column': j+1, 'order': k,
                   'kind': 'signed' if k == 1 else 'positive',
                   'pi_times_scalar_fiber_moment': str(raw.re),
                   'pi_cubed_times_3d_moment': str(physical_coeff),
                   'exact_expected': str(expected), 'passed': True}
            if k in (0, 2):
                Dinv = 1/F(basis_dets[str(j+1)])
                bound_coeff = Dinv*t**(k-3)
                require(0 <= physical_coeff <= bound_coeff, 'Positive conditional moment bound')
                row['pi_cubed_times_Mjk_majorant'] = str(bound_coeff)
            results.append(row)
    z = [s*basis[0][i]+y2*basis[1][i]+y3*basis[2][i] for i in range(3)]
    output = {'schema': 'B076-parallel-quotient-results-v1', 'status': 'PASS',
              'actual_configurations': 1, 'columns': cfg['columns'],
              'point': list(map(str, z)), 't': str(t), 'basis_determinants': basis_dets,
              'moment_results': results, 'python': sys.version,
              'platform': platform.platform(), 'dependencies': ['Python standard library only'],
              'scope': 'Exact fixed 3D coarea moments with parallel outside columns; no Green integral, no semiconcavity constant, no body32/full6 kernel.'}
    (root/'results.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'actual_configurations': 1, 'exact_moment_identities': 6}))


if __name__ == '__main__':
    main()
