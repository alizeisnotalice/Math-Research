# F063: geometry retained, implementation not ready for release

Reviewed seal: C063 preparation manifest `45f2defd70681f5d4cc826c58d2062f688bded2aca1d16dba7b487f22f4d287c`. All thirteen pinned preparation files were unchanged during the independent run. The actual prepared snapshot has the correct full root, no production release, no inflight job and no production call. This review does not invalidate that factual zero-production record or assert that a false global g result has already been produced.

The result is **REJECT_IMPLEMENTATION_FOR_RELEASE** because the public append/recovery path accepts forbidden state transitions. The underlying closed-region geometry and the documented conditional analytical certificate remain valid under their proper hypotheses.

## 1. Independent geometry proof

Let D={1>=y1>=...>=y5>=0}, and B=product_i[l_i,u_i] contained in [0,1]^5. Define L_i=max_(j>=i) l_j and U_i=min_(j<=i) u_j. If y is in B intersect D, then y_i>=l_j for every j>=i and y_i<=u_j for every j<=i. Hence L_i<=y_i<=U_i. Conversely the tightened rectangle B* is contained in B. Thus B intersect D=B* intersect D exactly.

Both L and U are descending. If each L_i<=U_i, the vector L itself is in B* intersect D; otherwise the offending coordinate inequality makes the intersection empty. The emptiness test is therefore necessary and sufficient, with strict L_i>U_i only. Equality, lower-dimensional ties and zero coordinates must remain. Applying the same maxima/minima to the already descending bounds leaves them unchanged. Their midpoint is descending and is an eligible face anchor.

When U_k>L_k, splitting its interval into the two closed halves at m=(L_k+U_k)/2 has union B*. Intersecting with D and tightening each child preserves that union. In particular every point with y_k=m remains in both applicable closed children; no seam is discarded. A truly empty child may be marked empty because its D intersection is empty, not because its volume is zero.

Starting from the exact root [0,1]^5, replacing one current leaf by the two registered children preserves the root union. Merely changing a leaf's processing status does not change its region. Therefore passed, queued, inflight and blocked nonempty leaves cover D after every **legal** transition, with exact empty witnesses for the remaining leaves. This induction requires the original root and the already-created topology to persist; it does not justify arbitrary states that only satisfy local node equations.

The independent asymmetric-box configuration verifies one different exact instance, including a seam point with repeated and zero coordinates. Its finite checks support implementation conformance, not the general proof above.

## 2. Analytical and same-window contracts

For each positive registered a and fixed T=1/2, the A057 inequality applies on the complete tightened rectangle, including S2/S3 and projection collisions. At every corner d the robust quantity is f_d exp(b_d/f_d-|d|^2/(32a^2)), with b_d the endpoint-minimized gradient dot product and f_d=clamp(b_d,[L,U]). The quadratic is one fixed matrix, so corner minimization follows from concavity before interval minimization. No tiny sub-box may replace the stored node.

The ladder (1/32,1/64,1/128,1/250,1/1000) is strictly decreasing and stays in (0,T). A strict center failure F_upper<target permits the next smaller cutoff but cannot reuse the old window's value or gradient. A strict center lower endpoint above target permits an actual whole-box test; failure of that test may subdivide at the same cutoff. Target overlap/equality and terminal ladder failure stop unresolved. These policies are conservative and are implemented on the ordinary reserve/ingest path. The independent wrong-window configuration is rejected by the adapter.

The external authentication boundary is stated clearly in theory.md and state_format.md. `check_saved_box` verifies rational corner data, the shared-value clamp and parsed-certificate digest; it deliberately does not prove that arbitrary endpoint numbers enclose an exponential. `validate_jet` checks metadata and hash-shaped identifiers; that is not native-result authentication. Main must inspect the actual pinned producer, native output, original receipt, A057 checker execution and relevant source/audit files. This limitation is not itself an undisclosed mathematical defect. In particular, a fabricated hash string alone is never evidence of a genuine gradient or an Arb exponential enclosure.

The parent reports final F061 acceptance, with manifest `91b26762f71229fbaa4c011be12b5a678ce968afc298e23a2c5320403c51fbda`, for two actual cases and all twelve full-gradient fields. The frozen C063 preparation predates that gate and correctly records no release. A future release must bind the final accepted evidence separately. Neither those finite cases nor the source proof means every future numerical input has already been independently rerun. This review performs no F061 rerun.

## 3. Blocking implementation findings: valid public API chains

These findings use `append_snapshot` itself with the correct current predecessor hash. They do not edit published bytes, invent predecessor digests, bypass canonical serialization, or rely merely on the obvious possibility of arbitrary external file forgery. Each resulting chain is accepted by `recover`.

### R1 — the validator does not fix the root it claims to cover

Location: author `cover_core.py` lines 338–350 and 388–390; cross-snapshot append lines 427–430.

`validate_state` requires an `r` node with no parent, but does not require raw root [0,1]^5 or depth zero. In configuration 2, the initial correct state is published, then a state whose unresolved root is rebuilt with `make_node('r',None,0,[[0,1/2]]*5,0)` is appended. Both append and recovery accept it, and the returned structural claim still has `root_partition_preserved=True`. The point y=(1,1,1,1,1) belongs to D but no longer belongs to the root region.

The correct initial snapshot is not disputed. The defect is that the advertised recovery verifier does not preserve that initial coverage invariant when an unresolved node is changed. A root-shrink state is not a legal geometric transition, yet the public persistence boundary accepts it.

Required repair: enforce the exact registered root geometry and identity on every state, and preserve existing node geometry/topology across snapshots according to legal split transitions. A full-root local check is necessary even when no prior snapshot is available; history validation is necessary for later topology.

### R2 — a hash-valid rollback drops an orphan and refunds its reservation

Location: author `cover_core.py` lines 393–414 and 417–430. Only the input-manifest identity and previously passed nodes are compared across snapshots.

Configuration 3 performs the following public sequence:

1. Publish initial state s0.
2. `reserve(s0)` creates inflight token `r.a0.n000001`, exact input digest and reservation count 1; publish this s1.
3. Call `append_snapshot(directory,s0,current_hash_of_s1)` using the stale initial state object but the correct current predecessor hash.
4. Recovery accepts this third snapshot, returns no inflight job, resets the reservation count to 0 and recommends continuing from the unique queue.
5. `reserve(recovered_state)` issues exactly `r.a0.n000001` with exactly the original input digest again.

This models the forbidden recovery window after a reservation checkpoint but before an original handle/receipt was attached. No actual process was launched in the test. The persisted state must treat an unknown original launch as unresolved and must not make the same token dispatchable again. The current append/recover path cannot distinguish the stale rollback from an admissible next state.

Ordinary reserve/ingest transitions enforce the eight-reservation cap, as independently checked with eight distinct synthetic nodes and a rejected ninth reservation. That local cap does not protect the cumulative batch when the history can be rolled back. The same gap can erase unresolved topology, original attempts, failure receipts, and a consumed native budget.

Required repair: validate legal transitions between snapshots or persist an operation/event that deterministically derives the next state from the previous one. At minimum preserve node topology, attempt-history prefixes, unique tokens, original handles/receipts, nondecreasing serials/reservations and started-node history; an inflight attempt may leave the active slot only through a validated completion or retained failure transition. Monotone counts alone are insufficient if an old token/input can be reissued under a new ledger.

### R3 — a missing release becomes dispatchable through snapshot metadata change

Location: author `cover_core.py` lines 120–134, 149–152, 338–387 and 406–409.

The `authorize` helper explicitly rejects TEST_ONLY state, but recovery does not validate the environment transition or require a release object when environment is MAIN_RELEASED. Configuration 7 publishes a TEST_ONLY initial state, then appends a copy with `environment='MAIN_RELEASED'` and `release=None`. The correctly linked chain recovers, and `reserve` succeeds without any release or F061 gate object.

This is not a successful native launch; no launcher exists. Main's independent external-authentication contract could catch it before launch. Nonetheless the state API currently promises a preparation/test/production distinction that its own public persistence-plus-reservation path does not enforce. That distinction must not rely solely on callers remembering to use `authorize`.

Required repair: validate environment values and release consistency in every state and in `reserve`; allow only the explicitly reviewed preparation-to-release transition; forbid TEST_ONLY promotion across all public state paths. Bind immutable release identity to the state history and retain the external-file verification requirement.

## 4. Safeguards that did pass

The independent legal-path configurations confirm the eight distinct-node/eight-attempt initial batch cap, same-window rejection after a cutoff decrease, and no revisit of a passed leaf. A previously passed node is protected against structural changes by `append_snapshot`. Synthetic completion leaves global_g_certificate and new_C6_bound false.

Canonical serialization, sequence/hash links, exclusive append and atomic publication are useful integrity mechanisms. They establish the provenance of the submitted records, not the legality of state transitions between those records. The three blockers arise after valid canonical publication rather than through breaking the file hash chain.

## 5. Release decision and repair boundary

Retain the pure geometry theorem, the A057 analytic interface, the documented external-proof boundary, the current zero-production preparation and the successful local safeguards. Do not release production on this core version. A repaired version should bind root and immutable metadata, validate transition legality, preserve attempts/reservations/topology and reject test-to-production promotion. It needs a new immutable seal and a bounded independent continuation; this audit does not edit or silently repair the frozen source.

No synthetic tree or successful state test is an actual body32 cover. A later complete tree still needs genuine same-window jets, authenticated exponential certificates, settled original receipts, an exact final partition audit, and main's separate symmetry/homogeneity integration. No global g or C6 conclusion is made here.
