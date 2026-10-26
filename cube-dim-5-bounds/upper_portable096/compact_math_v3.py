"""Exact compact-interval routines shared by the 096 packet builder only.

The independent portable checker reimplements the arithmetic separately.
"""
from fractions import Fraction as F
from itertools import product

SCALE_BITS = 896
SERIES_TERMS = 224
SCALE = 1 << SCALE_BITS
G = F(21, 200000)
SUPPORTED_EPSILONS = (F(1, 64), F(1, 128))


def ceildiv(a: int, b: int) -> int:
    return -(-a // b)


def fraction(text: str) -> F:
    if not isinstance(text, str):
        raise TypeError("source rationals must be strings")
    value = F(text)
    if str(value) != text:
        raise ValueError(f"noncanonical source rational {text!r}")
    return value


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def ceil_fraction(x: F) -> int:
    return ceildiv(x.numerator, x.denominator)


def quantize_pair(pair: list[str] | tuple[str, str], bits: int) -> list[str]:
    scale = 1 << bits
    lower, upper = map(fraction, pair)
    if lower > upper:
        raise ValueError("source interval endpoints are reversed")
    low_int = floor_fraction(lower * scale)
    high_int = ceil_fraction(upper * scale)
    width_int = high_int - low_int
    if width_int < 0:
        raise ArithmeticError("outward interval quantization produced negative width")
    return [str(low_int), str(width_int)]


def quantize_values(cache_record: dict, bits: int) -> list[list[str]]:
    values = cache_record["values"]
    if len(values) != 5 or any(len(pair) != 2 for pair in values):
        raise ValueError("source cache record must contain five interval pairs")
    return [quantize_pair(pair, bits) for pair in values]


def decode_values(encoded: list[list[str]], bits: int) -> list[tuple[F, F]]:
    scale = 1 << bits
    if len(encoded) != 5:
        raise ValueError("compact record must contain five intervals")
    result = []
    for pair in encoded:
        if len(pair) != 2 or any(not isinstance(x, str) or str(int(x)) != x for x in pair):
            raise ValueError("quantized intervals use canonical decimal integers")
        lo_i, width_i = map(int, pair)
        if width_i < 0:
            raise ValueError("quantized interval width must be nonnegative")
        result.append((F(lo_i, scale), F(lo_i + width_i, scale)))
    return result


def log_atanh(z: F) -> tuple[int, int]:
    if not F(0) <= z <= F(1, 3):
        raise ValueError("atanh argument outside [0,1/3]")
    zl = z.numerator * SCALE // z.denominator
    zu = ceildiv(z.numerator * SCALE, z.denominator)
    z2l = zl * zl // SCALE
    z2u = ceildiv(zu * zu, SCALE)
    pl, pu = zl, zu
    lo = hi = 0
    for j in range(SERIES_TERMS):
        d = 2 * j + 1
        lo += 2 * pl // d
        hi += ceildiv(2 * pu, d)
        pl = pl * z2l // SCALE
        pu = ceildiv(pu * z2u, SCALE)
    hi += ceildiv(2 * pu * SCALE, (2 * SERIES_TERMS + 1) * (SCALE - z2u))
    return lo, hi


LOG2 = log_atanh(F(1, 3))


def log_fixed(x: F) -> tuple[int, int]:
    if x <= 0:
        raise ValueError("log argument is nonpositive")
    e = x.numerator.bit_length() - x.denominator.bit_length()
    y = x / (F(2) ** e)
    if y < 1:
        e -= 1
        y *= 2
    if y >= 2:
        e += 1
        y /= 2
    if not F(1) <= y < F(2):
        raise ArithmeticError("fixed-point log normalization failed")
    lo, hi = log_atanh((y - 1) / (y + 1))
    if e >= 0:
        return lo + e * LOG2[0], hi + e * LOG2[1]
    return lo + e * LOG2[1], hi + e * LOG2[0]


def root_catalog() -> list[tuple[int, int, int, int]]:
    return [(i, j, k, l)
            for i in range(32)
            for j in range(i + 1)
            for k in range(j + 1)
            for l in range(k + 1)]


def reconstruct_box(root: list[int] | tuple[int, ...], path: str) -> tuple[list[F], list[F]]:
    if len(root) != 4 or any(type(x) is not int for x in root):
        raise ValueError("root requires four integer coordinates")
    if any(x < 0 or x > 31 for x in root) or any(root[i] < root[i + 1] for i in range(3)):
        raise ValueError("root coordinates must be nonincreasing in [0,31]")
    if not isinstance(path, str) or any(c not in "01" for c in path):
        raise ValueError("node path must contain only binary split bits")
    lo = [F(x, 32) for x in root]
    hi = [x + F(1, 32) for x in lo]
    for depth, bit in enumerate(path):
        j = depth % 4
        mid = (lo[j] + hi[j]) / 2
        if bit == "0":
            hi[j] = mid
        else:
            lo[j] = mid
    return lo, hi


def mode_from_certificate(cert: dict) -> str:
    mechanism = cert.get("mechanism")
    quadratic = cert.get("quadratic")
    allowed_mechanisms = {"full_frame_16", "full_frame_16_new_target_g105", "075_full_frame_16"}
    allowed_quadratics = {"Q_full=I4/16", "I4/(16epsilon^2)"}
    if mechanism not in allowed_mechanisms or quadratic not in allowed_quadratics:
        raise ValueError(f"unsupported actual source certificate mode: mechanism={mechanism!r}, quadratic={quadratic!r}")
    return "full_frame_16"


def assess_compact(node: list, root: list[int], point_tail: list[str], epsilon: F, quant_bits: int, value_pairs: list[list[str]], mode: str) -> dict:
    """Return exact local proof status from outward intervals, with no source calls."""
    root_id, path = node
    if epsilon not in SUPPORTED_EPSILONS:
        raise ValueError(f"unsupported per-leaf epsilon: {epsilon}")
    lo, hi = reconstruct_box(root, path)
    if len(point_tail) != 4:
        raise ValueError("point_tail must have four exact coordinates")
    point = tuple(map(fraction, point_tail))
    if any(not a <= x <= b for a, x, b in zip(lo, point, hi)):
        raise ValueError("anchor point lies outside its root/path box")
    vals = decode_values(value_pairs, quant_bits)
    Hlo, Hhi = vals[0]
    if not F(0) < Hlo <= Hhi:
        return {"passed": False, "reason": "quantized H lower endpoint is not positive"}
    if Hlo <= G:
        return {"passed": False, "reason": "quantized H lower endpoint does not exceed g"}
    delta_i, _ = log_fixed(Hlo / G)
    delta = F(delta_i, SCALE)
    gl = []
    for dlo, dhi in vals[1:]:
        ratios = (dlo / Hlo, dlo / Hhi, dhi / Hlo, dhi / Hhi)
        gl.append((min(ratios), max(ratios)))
    lower = tuple(a - x for a, x in zip(lo, point))
    upper = tuple(b - x for b, x in zip(hi, point))
    def linear(v):
        return sum((min(a * x, b * x) for x, (a, b) in zip(v, gl)), F(0))
    if mode == "Q=I4/4-J4/36":
        offsets = []
        for bits in product((0, 1), repeat=4):
            v = tuple(upper[j] if bits[j] else lower[j] for j in range(4))
            q = sum((x*x for x in v), F(0))/4 - sum(v, F(0))**2/36
            offsets.append(delta + linear(v) - q/(epsilon*epsilon))
    elif mode == "signed_basis_81":
        axes = [sorted({lower[j], F(0), upper[j]}) for j in range(4)]
        if any(len(axis) != 3 for axis in axes):
            return {"passed": False, "reason": "signed-basis point grid is degenerate"}
        offsets = []
        for v in product(*axes):
            q = sum((x*x for x in v), F(0))/4 - sum((abs(x) for x in v), F(0))**2/36
            offsets.append(delta + linear(v) - q/(epsilon*epsilon))
    else:
        offsets = []
        for bits in product((0, 1), repeat=4):
            v = tuple(upper[j] if bits[j] else lower[j] for j in range(4))
            q = sum((x*x for x in v), F(0))/16
            offsets.append(delta + linear(v) - q/(epsilon*epsilon))
    minimum = min(offsets)
    return {"passed": minimum > 0, "minimum": minimum, "delta_lower": delta}
