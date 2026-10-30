"""The final dataset: every metric checked, with provenance and certified answers.

    python validation/build_dataset.py

Reads validation/results/ and writes dataset.jsonl (one record per distinct
metric) and dataset.md (summary). Each record has:

  isosig, sqlens, n_tet
  sources      where the metric comes from: the 2020 candidate files
               (examples_*, hyperbolic_candidate), the full census scan
               (census_full_metrics), Pachner walks (census_walks*), controls
  manifold     census name and geometry where known
  flat         every edge's dihedral sum is exactly 2*pi
  answer       locally CAT(0) or not, with the reason for False
  pi, nopi     per setting (pi_beads True / False): answer, Python kernel
               verdict, Lean (vlinkcheck) verdict, certificate file on taiji1
"""

import glob
import json
import os
import sys
from collections import Counter, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')


def load(name):
    path = os.path.join(RES, name)
    return [json.loads(line) for line in open(path)] if os.path.exists(path) else []


def key(r):
    return r['isosig'], tuple(r['sqlens'])


def main():
    data = OrderedDict()
    for r in load('all_metrics2.jsonl'):
        data[key(r)] = dict(isosig=r['isosig'], sqlens=list(r['sqlens']),
                            n_tet=len(r['isosig']) and None, sources=[r['name']], manifold={})
    for r in load('census_full_metrics.jsonl'):
        d = data.setdefault(key(r), dict(isosig=r['isosig'], sqlens=r['sqlens'], sources=[], manifold={}))
        d['sources'].append('census %s (%s)' % (r['census'], r['shape']))
        d['manifold'].setdefault('name', r.get('name'))
    for path in sorted(glob.glob(os.path.join(RES, 'census_walks.*.jsonl'))) + \
            [os.path.join(RES, 'controls.jsonl')]:
        for r in map(json.loads, open(path)):
            if key(r) not in data:
                continue
            d = data[key(r)]
            src = os.path.basename(path)[:-len('.jsonl')]
            if r.get('from_isosig'):
                src += ' (walk from %s)' % r['from_isosig']
            if src not in d['sources']:
                d['sources'].append(src)
            if r.get('geometry'):
                d['manifold'].setdefault('geometry', r['geometry'])
            if r.get('name') and not r['name'].startswith('examples'):
                d['manifold'].setdefault('name', r['name'])
    for impl, tag in (('linkcond_cert', 'pi'), ('linkcond_cert_nopi', 'nopi')):
        for r in load('allm2.%s.jsonl' % impl):
            d = data[key(r)]
            lean = r.get('lean') or {}
            d[tag] = dict(result=r['result'], reason=r.get('reason'), python_verified=r['verified'],
                          lean_verified=lean.get('ok'), certificate=r.get('certificate_file'),
                          seconds=r.get('seconds'))
            d['flat'] = r.get('flat')
    from vlink.tri import Triangulation
    sys.path.insert(0, os.path.join(HERE, 'controls'))
    from topology import classify
    for d in data.values():
        name = d['manifold'].get('name')
        if name and 'geometry' not in d['manifold']:
            base = name.split(' : #')[0]
            if base.startswith('Hyp_') or base.startswith('HW_'):
                d['manifold']['geometry'] = 'hyperbolic'
            elif ' U/m ' in base or '] / [' in base:
                d['manifold']['geometry'] = 'graph manifold'
            else:
                d['manifold']['geometry'] = classify(base)[0]
    for impl in ('linkcond_cert',):
        for r in load('allm2.%s.jsonl' % impl):
            d = data[key(r)]
            # a closed manifold with a flat metric is a flat (Euclidean) manifold
            if r.get('flat') and d['manifold'].get('geometry') in (None, 'other'):
                d['manifold']['geometry'] = 'flat'

    for d in data.values():
        d['n_tet'] = Triangulation(d['isosig']).n
        pi, nopi = d.get('pi'), d.get('nopi')
        assert pi and nopi and pi['result'] == nopi['result'], d['isosig']
        assert pi['python_verified'] and pi['lean_verified'], d['isosig']
        assert nopi['python_verified'] and nopi['lean_verified'], d['isosig']
        d['answer'] = pi['result']
        d['reason'] = pi['reason']
    with open(os.path.join(RES, 'dataset.jsonl'), 'w') as f:
        for d in data.values():
            f.write(json.dumps(d) + '\n')

    rows = list(data.values())
    c = Counter((d['answer'], d['reason'], bool(d.get('flat'))) for d in rows)
    geo = Counter(d['manifold'].get('geometry', 'unknown') for d in rows)
    lines = ['# Dataset (validation/results/dataset.jsonl)', '',
             'Every metric examined: %d distinct (isosig, squared edge lengths), from the 2020'
             % len(rows),
             'candidate files, the full census scan (`census_full.md`), Pachner walks and the',
             'controls. Every answer is certified in both settings (`pi_beads` True and False):',
             'the certificate is accepted by the Python kernel and by `vlinkcheck` (Lean);',
             'certificates are stored on `/import/taiji1/wardak/certificates/allm2`.', '',
             '| answer | reason | flat | metrics |', '|---|---|---|---|']
    for (ans, reason, flat), n in sorted(c.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]), kv[0][2])):
        lines.append('| %s | %s | %s | %d |' % ('locally CAT(0)' if ans else 'not CAT(0)',
                                              reason or '', 'yes' if flat else 'no', n))
    lines += ['', 'Geometry of the manifold (where recorded): ' +
              ', '.join('%s %d' % kv for kv in geo.most_common()) + '.', '',
              'Tetrahedra per triangulation: %d to %d.' % (min(d['n_tet'] for d in rows),
                                                          max(d['n_tet'] for d in rows))]
    open(os.path.join(RES, 'dataset.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(HERE))
    main()
