# Pre-patch interval diagnostic

Six exact samples were collected for the originally numbered left interval 0.
The sixth-point quartic interpolation check failed, so no polynomial or
exclusion was emitted.  The cause was an interval-selection bug: when four
left events already existed, the code also retained the far `-1/16` window
endpoint, creating a first interval that crossed omitted events.  This was a
structural selection failure, not evidence for or against a bound.

The selector now uses exactly the four nearest intervals on each side (or the
window endpoint only if fewer events are available).  The pre-patch samples
are retained as diagnostics and are not used by the corrected audit.
