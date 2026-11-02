"""Collect every metric checked so far into one table of results.

    python validation/build_results_table.py

Inputs (validation/results/): the metric list all_metrics2.jsonl with names,
census_shape_metrics.jsonl and census_walks.*.jsonl for geometry and shape,
and the linkcond runs (allm2.linkcond*.jsonl). Outputs:
  results/examples.jsonl  one record per metric
  results/examples.md     tables grouped by geometry and result
"""

import glob
import json
import os
import re
from collections import Counter, defaultdict

import regina

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')


def load(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            r = json.loads(line)
            out[(r['isosig'], tuple(r['sqlens']))] = r
    return out


def shape_of(sqlens):
    s = set(sqlens)
    if s <= {1, 2, 3}:
        return 'cube'
    if s <= {3, 4}:
        return 'disphenoid'
    return 'other'


def volume_of(name):
    m = re.match(r'(?:Hyp|HW)_([0-9.]+)', name or '')
    return float(m.group(1)) if m else None


def canonical(name, isosig):
    """A key identifying the manifold across census naming schemes
    (Hyp_v (closed census) vs HW_v #k): volume plus first homology for
    hyperbolic manifolds (equal volumes do occur), the name otherwise."""
    v = volume_of(name)
    if v is None:
        return ('other', name)
    h1 = str(regina.Triangulation3.fromIsoSig(isosig).homology())
    return ('hyperbolic', round(v, 6), h1)


def census_names():
    """isosig -> manifold name, from the closed census and Hodgson-Weeks."""
    import sys
    sys.path.insert(0, os.path.join(HERE, 'controls'))
    from topology import load_census
    cdir = os.path.join(os.path.dirname(HERE), 'CAT0_Candidates')
    names = {sig: name for sig, name, _ in load_census(os.path.join(cdir, 'closed-or-census-11.dat'))}
    for line in open(os.path.join(cdir, 'hodgson-weeks.dat')):
        vol_k, _, sig = line.strip().partition(' ')
        vol, _, k = vol_k.partition('_')
        names.setdefault(sig, 'HW_%s #%s' % (vol, k))
    # the 2016 candidate files record retriangulations as 'isosig (from source)'
    for path in glob.glob(os.path.join(cdir, '*.txt')):
        for line in open(path):
            m = re.match(r'\s*(\S+) \(from (\S+)\)', line)
            if m and m.group(2) in names:
                names.setdefault(m.group(1), names[m.group(2)])
    return names


def main():
    meta = {}
    for path in [os.path.join(R, 'census_shape_metrics.jsonl')] + sorted(
            glob.glob(os.path.join(R, 'census_walks.*.jsonl'))):
        for k, r in load(path).items():
            meta.setdefault(k, r)
    metrics = load(os.path.join(R, 'all_metrics2.jsonl'))
    runs = {impl: load(os.path.join(R, 'allm2.%s.jsonl' % impl))
            for impl in ('linkcond', 'linkcond_nopi', 'linkcond_cert', 'linkcond_cert_nopi')}
    import sys
    sys.path.insert(0, os.path.join(HERE, 'controls'))
    from topology import classify
    cnames = census_names()
    rows = []
    for k, m in metrics.items():
        info = dict(meta.get(k, {}))
        if not info.get('name') and k[0] in cnames:
            info['name'] = cnames[k[0]]
        if info.get('name') and not info.get('geometry'):
            n = info['name']
            info['geometry'] = 'hyperbolic' if n.startswith(('Hyp_', 'HW_')) else classify(n)[0]
        name = info.get('name') or m.get('name')
        r = runs['linkcond'].get(k, {})
        cert = runs['linkcond_cert'].get(k) or runs['linkcond_cert_nopi'].get(k) or {}
        how = None
        if r.get('result') is False:
            how = 'witness (kernel-verified)' if r.get('verified') else 'witness (unverified)'
        elif r.get('result') is True:
            if runs['linkcond_cert'].get(k, {}).get('verified'):
                how = 'certificate (kernel-verified)'
            elif runs['linkcond_cert_nopi'].get(k, {}).get('verified'):
                how = 'certificate (kernel-verified, assuming EMM Lemma 4.9)'
            else:
                how = 'search only'
        rows.append(dict(
            isosig=k[0], sqlens=list(k[1]), name=name, geometry=info.get('geometry'),
            volume=volume_of(name), shape=info.get('shape') or shape_of(k[1]),
            n_tet=regina.Triangulation3.fromIsoSig(k[0]).size(), flat=info.get('flat', m.get('flat')),
            result=r.get('result'), established_by=how,
            certificate_bytes=cert.get('certificate_bytes')))
    rows.sort(key=lambda x: (x['geometry'] or '', x['volume'] or 0, x['name'] or '', x['n_tet']))
    with open(os.path.join(R, 'examples.jsonl'), 'w') as f:
        for row in rows:
            f.write(json.dumps(row) + '\n')

    # markdown summary
    by = defaultdict(list)
    for row in rows:
        by[(row['geometry'] or 'unknown', row['result'])].append(row)
    lines = ['# Metrics checked by vlink', '',
             '%d metrics; results %s; established by %s.' % (
                 len(rows), dict(Counter(r['result'] for r in rows)),
                 dict(Counter(r['established_by'] for r in rows))), '']
    manifolds = defaultdict(set)
    for row in rows:
        manifolds[(row['geometry'], row['result'])].add(canonical(row['name'], row['isosig']))
    lines += ['| geometry | locally CAT(0) (metrics / manifolds) | not CAT(0) (metrics / manifolds) |',
              '|---|---|---|']
    for g in sorted({r['geometry'] or 'unknown' for r in rows}):
        t, f = by.get((g, True), []), by.get((g, False), [])
        lines.append('| %s | %d / %d | %d / %d |' % (g, len(t), len(manifolds[(g, True)]),
                                                     len(f), len(manifolds[(g, False)])))
    for (g, res), items in sorted(by.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        lines += ['', '## %s — %s' % (g, {True: 'locally CAT(0)', False: 'not CAT(0)'}.get(res, res)), '',
                  '| manifold | volume | shape | tetrahedra | flat | established by | isosig | squared lengths |',
                  '|---|---|---|---|---|---|---|---|']
        for r in items:
            lines.append('| %s | %s | %s | %d | %s | %s | `%s` | %s |' % (
                r['name'], '%.6f' % r['volume'] if r['volume'] else '', r['shape'], r['n_tet'],
                'yes' if r['flat'] else '', r['established_by'] or '', r['isosig'], r['sqlens']))
    open(os.path.join(R, 'examples.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines[:12]))


if __name__ == '__main__':
    main()
