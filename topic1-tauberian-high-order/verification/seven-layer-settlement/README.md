# Seven-layer exact algebra

Run locally with Mathematica/WolframScript 15 ARM: wolframscript -file seven_layers.wl. The retained log is the actual successful run. This script proves finite identities and series under the supplied scalar interfaces; the geometric interfaces and all-input remainder control are established by the full proof and its independent review, not by this program.

Checks: telescoping residual is 0; first five layer weights start at v^2; normalized initial endpoint is O(v^8 S); the exact T=V+M weighted error agrees with Q(v) through v^7; solving the true ratio denominator gives coefficients {4/27,-176/243}. Both symbolic constants b and cD remain arbitrary positive fixed constants, so their cancellation is tested.

The early V-normalized intermediate series and the final exact T-normalized section are recorded together; only the latter is the interface used in the proof. Some non-ASCII printed label characters have encoding noise in the original log; rational expressions and the ASCII result labels are intact. No numerical fitting, interval-rounding certification, or geometric optimization is claimed.

seven_layers.wl SHA-256 876c0d8191c7c100a9402e652d2aefddf31a0fd3756a5df81d4402a4bdc8c37a

seven_layers.log SHA-256 62acf028e9467ec971eb1bb5e1b6b37cbfb80238eaaa8ad10542c84be80059db
