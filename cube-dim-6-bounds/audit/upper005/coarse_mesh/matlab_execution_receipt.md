# MATLAB execution receipt (one completed scalar cross-check)

This receipt records the single MATLAB call already made for this audit.  It
is historical metadata only: **do not rerun it**.  On subsequent instruction,
shared MATLAB is to be coordinated serially by the main agent; this audit will
not make further MATLAB calls.

## Tool metadata

- Tool: `mcp__matlab__evaluate_matlab_code`.
- Reported elapsed wall time: 17.1 seconds.
- Start and end timestamps: not supplied by the tool response, so unavailable.
- Project path passed to the tool:
  `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature`.
- Concurrency observation: this audit had no report of a simultaneous MATLAB
  job; this is an observation, not a locking guarantee.

## Exact executed code

```matlab
format long g; D = 25/192; eps = 1/100; L = 6e-8; U = 23e-6; q = 6/5; K = (6/7)*(1/243)^(1/6)*pi^(-11/6); threshold = (16/3)*U; coarse = L*exp(D/(32*eps^2)); tail = K*eps^(7/6); fprintf('D=%.17g\nK=%.17g\nthreshold=%.17g\ncoarse=%.17g\ntail=%.17g\nratios: coarse/threshold=%.17g tail/threshold=%.17g\n',D,K,threshold,coarse,tail,coarse/threshold,tail/threshold);
```

## Returned command-window output

```text
D=0.13020833333333334
K=0.04207370247458074
threshold=0.00012266666666666665
coarse=28160408990.533825
tail=0.00019528882759476022
ratios: coarse/threshold=229568551553264.91 tail/threshold=1.5920284858268499
```

## State effects

- The code created or overwrote workspace variables `D`, `eps`, `L`, `U`, `q`,
  `K`, `threshold`, `coarse`, and `tail` in the existing MATLAB session.
- It created no files, figures, plots, or graphics.
- The submitted MATLAB code did not run `clear`, `clearvars`, `close`,
  `restoredefaultpath`, `addpath`, `rmpath`, or `cd`.
- The MCP call supplied `project_path`, whose documented behavior is to set the
  current working folder to that path.  Thus the session working folder may
  have been changed by the tool wrapper; no MATLAB path search-path change was
  requested or performed.
