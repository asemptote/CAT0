"""Task functions for validation/workqueue.py ({"fn": "tasks:<name>", "args": {...}}).

Each takes JSON arguments and returns a JSON-serialisable result.
"""

import json
import os
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VLINKCHECK = os.path.join(ROOT, 'lean', '.lake', 'build', 'bin', 'vlinkcheck')


def run_check(impl, isosig, sqlens, **meta):
    """One input through an implementation of validation/run_check.py."""
    from run_check import IMPLS
    out = IMPLS[impl](isosig, sqlens)
    out.update(impl=impl, isosig=isosig, sqlens=sqlens, **meta)
    return out


def false_certificate(isosig, sqlens):
    from vlink.check import certificate
    return certificate(isosig, list(sqlens))['certificate']


def true_certificate(isosig, sqlens, pi_beads=True):
    from vlink.certify import true_certificate as tc
    return tc(isosig, list(sqlens), pi_beads=pi_beads)


def kernel_verdict(certificate, exact_only=True, **meta):
    """The Python kernel's verdict on a certificate of either kind."""
    from vlink import kernel
    if certificate.get('format') == kernel.TRUE_FORMAT:
        ok, msg = kernel.check_true_certificate(certificate, exact_only=exact_only)
    else:
        ok, msg = kernel.check_certificate(certificate)
    return dict(meta, ok=ok, msg=msg)


def lean_verdict(certificate, **meta):
    """vlinkcheck's verdict on one certificate."""
    with tempfile.NamedTemporaryFile('w', suffix='.jsonl', delete=False) as f:
        f.write(json.dumps(certificate) + '\n')
        path = f.name
    try:
        p = subprocess.run([VLINKCHECK, path], capture_output=True, text=True, check=True)
    finally:
        os.unlink(path)
    verdict, claim, msg = (p.stdout.strip().split(' ', 2) + ['', ''])[:3]
    return dict(meta, ok=verdict == 'ACCEPT', claim=claim, msg=msg)


def shape_metrics(isosig, shape, limit=100000, timeout=600.0, **meta):
    """All metrics of one shape satisfying the edge condition (or as many as
    found before limit / timeout): counts and the metrics themselves."""
    import sys
    sys.path.insert(0, os.path.join(HERE, 'controls'))
    from shape_metrics import enumerate_metrics
    metrics, complete = enumerate_metrics(isosig, shape, limit=limit, timeout=timeout)
    return dict(meta, isosig=isosig, shape=shape, n=len(metrics), complete=complete,
                metrics=metrics[:1000])


def shape_metrics_batch(isosigs, shape, limit=100000, timeout=3600.0, **meta):
    """shape_metrics for several triangulations (one queue task per batch)."""
    out = [shape_metrics(sig, shape, limit=limit, timeout=timeout) for sig in isosigs]
    return dict(meta, shape=shape, n_triangulations=len(out),
                incomplete=[r['isosig'] for r in out if not r['complete']],
                found=[r for r in out if r['n']])
