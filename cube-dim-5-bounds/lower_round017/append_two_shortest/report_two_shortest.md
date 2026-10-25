# 017 append: two shortest remaining streams, starting from2/10

Exactly two new streams were registered and completed sequentially. The
starting checkpoint had2/10 streams complete:0<->1 and1/2<->2. Neither was
rerun. All433 raw affine functions and1024 original cases remain in scope.
No union, parameter grid, fixed-ID partition, angle search or plane-capture
case was launched in this batch.

## Stream1: slopes1/2 and1

This was the shortest remaining stream,1x101=101 pairs. Exact classification
gives4 events left of the window,17 inside[3/2,2] and80 to its right.
The17 window positions are exactly a=N^(1/5), with

    N=[8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 26, 27, 28, 30, 32].

Their outward rational isolating intervals are saved individually. They
range from8^(1/5) to the exact upper endpoint2. Every one is formally
identical, after fifth-power normalization, to an event already in the
completed0/1 stream. Thus there are0 new unique window events. Both the
within-stream isolation and the cross-completed-event check have0 unresolved
clusters. This is an exact expression comparison plus interval separation,
not a comparison of rounded decimals.

These are raw endpoint equalities, not necessarily active breakpoints of
the physical union. Some derive from geometrically inactive cases, which
were deliberately retained in the complete raw catalog.

Actual execution: PID28680, nice10,
1 thread, RSS25133056 bytes,
0.012825250 seconds,101 nodes, exit0. No live
session handle remained when the next stream was registered.

## Stream2: slopes-1 and1/2

After stream1 completed, this was the new shortest remaining stream,
102x1=102 pairs. All102 events lie strictly left of[3/2,2]. There are no
window events and no unresolved cluster. Analytically, r-a=a/2 implies

    a=2r/3=N^(1/5)/3 <=4/3 <3/2.

Equivalently, throughout this window r-a<=1/2<a/2, so this raw upper endpoint
is never clipped by the domain endpoint a/2. This simplifies that branch
without asserting that its underlying family is otherwise nonempty.

Actual execution: PID28762, nice10,
1 thread, RSS25247744 bytes,
0.013553250 seconds,102 nodes, exit0. The process
exited. Both jobs used the six numerical-library thread variables fixed to1,
with no pool or application window.

## Updated structure and untouched work

Now4/10 streams are complete. The completed streams collectively still have
710 isolated unique formal potential window events; the new batch adds none.
This is not the complete event set: six streams remain

    (-1,0):16218 pairs; (-1,1):10302; (-1,2):7140;
    (0,1/2):159; (0,2):11130; (1,2):7070.

Their events were not evaluated. No ordering is asserted inside any future
unresolved cluster; the runner saves those clusters rather than forcing
separation.

Capacity diagnosis is unchanged: a*s^4 exceeds the complete frozen C_cert
throughout[3/2,2], so capacity alone excludes no parameter subinterval.
The prior uniform side12 cutoff remains recorded, but the fixed-reference
side4 completeness proof is not extended to this parameter window.

The original18 materialized-plane-packet list remains explicit. The earlier
analytic whole-plane qualification of(8,1) was read, not rerun, and no fixed-ID
partition was created. There are17 still-unresolved capture cases:

    (8,2), (27,10..16), (36,10..14), (48,8..11).

A downstream fixed-ID-only workflow may still materialize(8,1); its existing
whole-plane certificate uses a moving nearest-point selector and must not be
misrepresented as one fixed atom over the entire period.

All old017 checkpoint files, the previous append and the15-file007 package
were hash-verified unchanged. All new files live in this append directory.
No commit or formal lower-bound update was made. The batch stops here with
no live computational handle and no further stream or qualification launched.
