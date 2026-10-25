#!/usr/bin/env python3
"""Independent scalar audit of the primary local H_epsilon box certificate.

This reads the locked Arb vertex intervals but does not rerun their 967,680
cone-term evaluations.  It checks geometry, symmetry multiplicities, ODE loss,
face-half-box integral direction, and the conditional global cost arithmetic.
"""
from __future__ import annotations

from pathlib import Path
from hashlib import sha256
from math import comb
import json
from flint import arb, ctx

ctx.prec = 384


def main() -> None:
    primary = Path(__file__).resolve().parents[1]
    source = primary / 'local_box' / 'arb384_results.json'
    data = json.loads(source.read_text())
    rows = sorted(data['outputs'], key=lambda row: row['index'])
    assert [row['index'] for row in rows] == list(range(5))
    h = arb(1)/1000
    epsilon = arb(1)/100
    q = arb(6)/5
    lower = arb(3)/100000
    loss = ((arb(129)/128).sqrt()*h/epsilon).cosh()**5
    mean_vertices = arb(0)
    for k, row in enumerate(rows):
        assert row['epsilon'] == '1/100' and row['denominator'] == 1000
        assert sorted(row['point_numerators']) == sorted([1000]+[751]*k+[749]*(4-k)+[1])
        value = arb(row['H_epsilon_arb_384'])
        assert value > lower*loss
        # There are C(4,k) x2 fifth-axis vertices.  The full-box average has
        # denominator 32, so fifth-axis reflection leaves C(4,k)/16.
        mean_vertices += comb(4, k)*value**q/16
    # Integrate the full [-h,h] fifth axis by the log-semiconvex vertex bound,
    # then divide by 2 by reflection to obtain x6 in [0,h].
    local_h_power_upper = (2*h)**4*h * (q*5*h*h/(32*epsilon*epsilon)).exp() * mean_vertices
    assert local_h_power_upper < arb(7)/10**20
    tail_epsilon = arb(1)/1000
    tail = (arb(6)/7)*(arb(32)/243)**(arb(1)/6)*arb.pi()**(-arb(11)/6)*tail_epsilon**(arb(7)/6)
    global_g_lower = arb(1)/60000
    global_h_norm_upper = arb(1)/16000
    conditional_cost = 6**q*((global_h_norm_upper+tail)/global_g_lower)**q
    assert conditional_cost < 64
    result = {
        'status': 'PASS_scalar_audit_local_only_and_conditional_global_only',
        'source_sha256': sha256(source.read_bytes()).hexdigest(),
        'five_representative_vertices': True,
        'full_box_vertex_count': 32,
        'positive_half_box_factor': '1/2 by x6 reflection',
        'ode_loss': str(loss),
        'local_H_epsilon_lower_implies_local_G_lower': 'H_1/100 >= 3/100000 on x2..x5 in [749,751]/1000, x6 in [0,1]/1000',
        'local_H_power_integral_upper': str(local_h_power_upper),
        'warning': 'This is an H_(1/100) local integral upper, not an I or G integral upper.',
        'conditional_only': {
            'requires_g_ge': '1/60000',
            'requires_H_1_1000_norm_le': '1/16000',
            'cost_upper_if_both_requirements_hold': str(conditional_cost),
            'requirements_proved_here': False,
        },
        'source_code_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name('local_box_independent_audit.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
