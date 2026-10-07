import numpy as np
import math


def maximal_from_q(q, weights, n, t):
    """q_j = full side needed for atom j / source radial side.

    Evaluate only the END of each exact equal-distance group. Choosing the
    smallest radius in an exact response tie is deterministic.
    """
    order = np.argsort(q, axis=1, kind="stable")
    qs = np.take_along_axis(q, order, axis=1)
    mass = np.cumsum(weights[order], axis=1)
    last = np.concatenate([qs[:, :-1] != qs[:, 1:], np.ones((len(q), 1), dtype=bool)], axis=1)
    with np.errstate(divide="ignore"):
        logq = np.log(qs)
        response = np.log(mass) - t - n * logq
    response = np.where(last, response, -np.inf)
    winner = np.argmax(response, axis=1)
    rows = np.arange(len(q))
    return qs, logq, mass, last, response, winner, response[rows, winner]

def window_maximal(q, weights, n, t):
    """True continuous-window maximal function, full sides in [1,2]."""
    qs, logq, mass, last, _, _, _ = maximal_from_q(q, weights, n, t)
    r = math.exp(t/n)
    qa, qb = 1/r, 2/r
    inside = (qs >= qa) & (qs <= qb)
    endpoints = np.broadcast_to([qa, qb], (len(q), 2))
    endpoint_mass = np.stack([np.sum((q <= qa)*weights, axis=1),
                              np.sum((q <= qb)*weights, axis=1)], axis=1)
    qs = np.concatenate([qs, endpoints], axis=1)
    mass = np.concatenate([mass, endpoint_mass], axis=1)
    valid = np.concatenate([last & inside, np.ones((len(q), 2), dtype=bool)], axis=1)
    perm = np.argsort(qs, axis=1, kind='stable')
    qs = np.take_along_axis(qs, perm, axis=1)
    mass = np.take_along_axis(mass, perm, axis=1)
    valid = np.take_along_axis(valid, perm, axis=1)
    with np.errstate(divide='ignore'):
        logq = np.log(qs)
        resp = np.log(mass)-t-n*logq
    resp = np.where(valid, resp, -np.inf)
    ix = np.argmax(resp, axis=1)
    return qs, logq, mass, valid, resp, ix, resp[np.arange(len(q)), ix]
