"""Random equal-cell atoms: autocorrelation and low-frequency cross diagnostics.

FFT values are numerical quadrature, not an interval certificate or lower bound.
The concentration and expectation identities are proved separately.
"""
from pathlib import Path
import numpy as np,json,math

OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(202610290345);rows=[]
for K in [15,63,255,1023]:
    Q=max(65536,1<<(K*K-1).bit_length())
    u=rng.random(K);atoms=(np.arange(K)+u)/K
    x=(np.arange(Q)+.5)/Q
    p=np.searchsorted(atoms,x,side='right')-K*x
    spectrum=np.fft.rfft(p)
    correlation=np.fft.irfft(abs(spectrum)**2,n=Q)/Q
    t=np.arange(Q)/Q
    expected=np.maximum(0.,1-K*np.minimum(t,1-t))**3/6
    error=float(np.max(abs(correlation-expected)))
    theorem_envelope=4*math.sqrt(math.log(K)/K)+3/K
    assert error<theorem_envelope
    m=np.arange(1,min(16,K-1)+1)
    b=4*(729/29)*((np.sin(np.pi*m*10/27)/(np.pi*m))**4
                  +3*(np.sin(np.pi*m/9)/(np.pi*m))**4)
    c=np.sqrt(2*b/3)
    z=np.exp(-2j*np.pi*atoms[:,None]*m).sum(axis=0)
    phat=z/(2j*np.pi*m)
    cross_envelope=float(4*np.sum(c*abs(phat)))
    markov_threshold=4*(20/9)*math.pi**2/(3*math.sqrt(K))
    rows.append(dict(K=K,Q=Q,seed=202610290345,
       P_abs_max=float(abs(p).max()),autocorrelation_error_grid=error,
       BV_quadrature_error_budget_without_FFT_rounding=4*K/Q,
       uniform_concentration_envelope=theorem_envelope,
       concentration_failure_upper=2/K**6,
       low_frequency_count=len(m),actual_fourier_cross_envelope=cross_envelope,
       four_times_expected_cross_envelope=markov_threshold,
       source_mean_P=float(p.mean())))
(OUT/'results.json').write_text(json.dumps(dict(numpy=np.__version__,rows=rows,
    scope='Deterministic-seed numerical diagnostics only; no new C_n lower bound.'),indent=2))
print(json.dumps(rows,indent=2))
