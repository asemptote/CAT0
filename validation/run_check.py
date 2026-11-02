"""Run a checker implementation over the candidate examples and record witnesses.

    python validation/run_check.py --impl current --workers 4
    python validation/run_check.py --impl legacy2020 --files examples_1

Each input is run in its own process; one JSON line per input is appended to
validation/results/<impl>.jsonl with the result, the vertex, the necklace
(beads as link-triangle indices plus start/end link-triangle vertex indices),
and the wall time. Identical (isosig, sqlens) inputs are run once.

impl = current:   gallery.py (LinkTriangles / build_galleries / find_necklaces),
                  driven exactly as vertex_linkcond(mode='cert') drives them,
                  but keeping the galleries so the witness can be recorded.
impl = legacy2020: validation/legacy/vlinkcond_2020.py, a mechanical Python 3 /
                  Regina 7 port of the code that produced the Sept 2020 results.
"""

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, ROOT, os.path.join(HERE, 'legacy')]

from examples import FILES, load_all


def _float(x):
    try:
        return float(x)
    except Exception:
        return repr(x)


def run_current(isosig, sqlens):
    import gallery
    gallery.DEBUG = False
    T = gallery.fromIsoSig(isosig)
    for i_v, v in enumerate(T.vertices()):
        L = v.buildLink()
        lt = gallery.LinkTriangles(T, v, sqlens)
        G = []
        for t in range(L.triangles().size()):
            G.extend(gallery.build_galleries(lt, t))
        iss_G = tuple([] for _ in range(L.vertices().size()))
        for i_G, g in enumerate(G):
            iss_G[L.triangle(g.ts[0]).vertex(g.ds[0]).index()].append(i_G)
        for i in range(len(G)):
            necklaces = gallery.find_necklaces(G, lt, iss_G, False, None, None, (i,))
            if necklaces:
                beads = [dict(ts=[int(t) for t in G[j].ts], start=int(G[j].ds[0]),
                              end=int(G[j].rs[-1]),
                              length=_float(G[j].geodesiclength))
                         for j in necklaces[0]]
                return dict(result=False, vertex=i_v, necklace=beads, n_galleries=len(G))
    return dict(result=True)


def run_legacy2020(isosig, sqlens):
    import vlinkcond_2020 as vl
    v, loop, L, T = vl.vertexLinkCondition(isosig, sqlens)
    if v is True:
        return dict(result=True)
    beads = [dict(ts=[t.index() for t in g[0]], start=int(g[1]), end=int(g[2]))
             for g in loop]
    return dict(result=False, vertex=v.index(), necklace=beads)


def run_vlink(isosig, sqlens):
    """vlink with its defaults (link triangles may be reused)."""
    from vlink.check import vertex_linkcond
    return vertex_linkcond(isosig, sqlens)


def run_vlink_norepeat(isosig, sqlens):
    """vlink restricted to gallery.py's search space (no reused triangles)."""
    from vlink.check import vertex_linkcond
    return vertex_linkcond(isosig, sqlens, allow_repeats=False)


def run_vlink_repeats(isosig, sqlens):
    from vlink.check import vertex_linkcond
    return vertex_linkcond(isosig, sqlens, allow_repeats=True)


def run_vlink_stopstart(isosig, sqlens):
    """vlink with gallery.py's rule of not extending paths past their start."""
    from vlink.check import vertex_linkcond
    return vertex_linkcond(isosig, sqlens, stop_at_start=True)


def run_linkcond(isosig, sqlens):
    """vlink's full test (validity, edges, vertex links), every False answer
    verified by the kernel; certificates kept in the output."""
    from vlink.check import linkcond
    r = linkcond(isosig, sqlens)
    r.pop('necklace', None)
    return r


def run_linkcond_nopi(isosig, sqlens):
    """linkcond without beads of length exactly pi (Lemma 4.9 of EMM 2003)."""
    from vlink.check import linkcond
    r = linkcond(isosig, sqlens, pi_beads=False)
    r.pop('necklace', None)
    return r


def _certified(isosig, sqlens, **kw):
    """linkcond with certificates for True answers too, checked by the Python
    kernel and (when built) by vlinkcheck, the Lean reference checker, timed
    separately. Certificates are summarised (size); with VLINK_CERT_DIR set
    they are also stored there (gzipped JSON)."""
    import gzip
    import hashlib
    from vlink import lean
    from vlink.check import linkcond
    r = linkcond(isosig, sqlens, certify=True, lean=False, **kw)
    r.pop('necklace', None)
    cert = r.pop('certificate', None)
    r['certificate_bytes'] = len(json.dumps(cert)) if cert else 0
    if cert is not None and lean.binary():
        t0 = time.time()
        ok, _, msg = lean.check([cert])[0]
        r['lean'] = dict(ok=ok, msg=msg, seconds=round(time.time() - t0, 2))
        if ok != r['verified']:
            r['lean_disagrees'] = True
    d = os.environ.get('VLINK_CERT_DIR')
    if d and cert is not None:
        key = hashlib.sha1(json.dumps([isosig, list(sqlens), kw], sort_keys=True).encode())
        path = os.path.join(d, key.hexdigest()[:16] + '.json.gz')
        with gzip.open(path, 'wt') as f:
            json.dump(cert, f)
        r['certificate_file'] = path
    return r


def run_linkcond_cert(isosig, sqlens):
    return _certified(isosig, sqlens)


def run_linkcond_cert_nopi(isosig, sqlens):
    return _certified(isosig, sqlens, pi_beads=False)


IMPLS = dict(current=run_current, legacy2020=run_legacy2020, vlink=run_vlink,
             linkcond=run_linkcond, linkcond_nopi=run_linkcond_nopi,
             linkcond_cert=run_linkcond_cert, linkcond_cert_nopi=run_linkcond_cert_nopi,
             vlink_norepeat=run_vlink_norepeat, vlink_repeats=run_vlink_repeats,
             vlink_stopstart=run_vlink_stopstart)


def job(args):
    impl, ex = args
    t0 = time.time()
    try:
        out = IMPLS[impl](ex['isosig'], ex['sqlens'])
    except Exception as e:  # record, don't abort the batch
        out = dict(result=None, error='%s: %s' % (type(e).__name__, e))
    out.update(impl=impl, seconds=round(time.time() - t0, 1), **ex)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--impl', choices=sorted(IMPLS), required=True)
    ap.add_argument('--files', nargs='*', default=list(FILES))
    ap.add_argument('--jsonl', default=None,
                    help='take inputs from a JSONL file of {isosig, sqlens, ...} instead of --files')
    ap.add_argument('--only-nonflat', action='store_true',
                    help='with --jsonl: skip records marked flat')
    ap.add_argument('--chunk', type=int, default=0,
                    help='run only inputs chunk, chunk+nchunks, ... (PBS array index)')
    ap.add_argument('--nchunks', type=int, default=1)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(HERE, 'results', a.impl + '.jsonl')
    os.makedirs(os.path.dirname(out), exist_ok=True)

    # resume: skip inputs that already have a (non-error) result in `out`
    seen, jobs = set(), []
    if os.path.exists(out):
        for line in open(out):
            r = json.loads(line)
            if r.get('result') is not None:
                seen.add((r['isosig'], tuple(r['sqlens'])))
    n_done = len(seen)
    if a.jsonl:
        exs = []
        for i, line in enumerate(open(a.jsonl)):
            r = json.loads(line)
            if a.only_nonflat and r.get('flat'):
                continue
            exs.append(dict(file=os.path.basename(a.jsonl), idx=i, isosig=r['isosig'],
                            sqlens=r['sqlens'], raw_sqlens=r['sqlens'],
                            name=r.get('name')))
    else:
        exs = load_all(a.files)
    for ex in exs:
        key = (ex['isosig'], tuple(ex['sqlens']))
        if key not in seen:
            seen.add(key)
            jobs.append((a.impl, ex))
    # longest (most tetrahedra) first so the big ones don't finish last
    jobs.sort(key=lambda j: -len(j[1]['sqlens']))
    jobs = jobs[a.chunk::a.nchunks]
    print('%d distinct inputs to run (%d already done)' % (len(jobs), n_done), flush=True)
    with Pool(a.workers, maxtasksperchild=1) as p, open(out, 'a') as f:
        for r in p.imap_unordered(job, jobs):
            f.write(json.dumps(r) + '\n')
            f.flush()
            print(r['file'], r['idx'], r['isosig'], r['result'],
                  r.get('error', ''), '%.0fs' % r['seconds'], flush=True)


if __name__ == '__main__':
    main()
