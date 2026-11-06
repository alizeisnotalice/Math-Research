# r029 formal audit — 2026-11-06

Read-only independent audit of `rounds/r029/report.md`, sections 2--6.

**Accepted in the stated `R_N` feedback-model scope.**  It is not a new
lower bound for the strong-maximal constant `C_N`.

The triangular signed-tilt Fourier estimate covers both signs and provides
the needed all-shift exponential factor.  At `r=tm+x`, the future argument
is exactly `r-m=(t-1)m+x`; (6) supplies the required density replacement.
The action cancellation in (15) is genuine: with
`p_c=(1+e^c)^-1` and `m(c)=1+cp_c`, one has
`m'(c)=-p_c(c-m(c))`, so the first parameter derivative is controlled by
`|b_N|U' + U''`.

For (14), `j<=t/2` is controlled by a common signed tilt, whose factors
cancel and leave `exp(C|u|/sqrt(t))`; for `j>t/2`, differentiating the
smooth source and convolving against the prefix density gives the same
bound for every real shift.  The window `T=N/L^2` makes the resulting
action error negligible compared with the source.

The exact sign check is
`W_t-K_{c_t}W_{t-1}=G_t+E_t-g_t`.
The reserves in (10), (11), and (15) make it nonnegative, so the induction
does not require an occupation estimate or prior positivity of `W`.
Finally, the common-tilt conditional-density cancellation in (18), the
`t-1/t` bridge estimates, and the harmonic source sum are consistent.
