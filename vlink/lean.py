"""The Lean checker (lean/: vlinkcheck) as vlink's reference checker.

vlinkcheck re-checks certificates of both kinds with exact arithmetic whose
soundness is proved in Lean (lean/README.md). It is about ten times faster
than the Python kernel with exact_only=True, so check.linkcond runs it on
every certificate when it is available (built, or pointed to by the
environment variable VLINKCHECK) and records its verdict as 'lean'.
Nothing here is trusted: it only runs the binary and reads its output.
"""

import json
import os
import subprocess
import tempfile

_DEFAULT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'lean', '.lake', 'build', 'bin', 'vlinkcheck')


def binary():
    """Path of vlinkcheck, or None if it is not built."""
    path = os.environ.get('VLINKCHECK', _DEFAULT)
    return path if os.path.isfile(path) and os.access(path, os.X_OK) else None


def check(certificates):
    """[(ok, claim, message)] for a list of certificates (one vlinkcheck run)."""
    exe = binary()
    if exe is None:
        raise FileNotFoundError('vlinkcheck is not built (lean/README.md)')
    with tempfile.NamedTemporaryFile('w', suffix='.jsonl', delete=False) as f:
        for c in certificates:
            f.write(json.dumps(c) + '\n')
        path = f.name
    try:
        out = subprocess.run([exe, path], capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(path)
    rows = []
    for line in out.splitlines():
        verdict, claim, msg = (line.split(' ', 2) + ['', ''])[:3]
        rows.append((verdict == 'ACCEPT', claim, msg))
    if len(rows) != len(certificates):
        raise RuntimeError('vlinkcheck returned %d verdicts for %d certificates'
                           % (len(rows), len(certificates)))
    return rows
