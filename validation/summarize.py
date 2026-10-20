"""Compare implementations and verify every witness independently.

    python validation/summarize.py > validation/results/summary.md

Reads validation/results/*.jsonl (one file per implementation / input set),
runs verify_cert.metric_report on every input and verify_cert.verify_necklace
on every False result, and prints a markdown table plus a verdict per input:

  expected True  -- the metric is flat (locally Euclidean, hence CAT(0))
  witness valid  -- a False result whose necklace is a genuine short geodesic
  witness INVALID-- a False result whose necklace is not (implementation bug)
  invalid input  -- the input is not a Euclidean metric at all
"""

import glob
import json
import os
from collections import defaultdict

import verify_cert as V

HERE = os.path.dirname(os.path.abspath(__file__))


def short(r):
    if r.get('result') is None:
        return 'error'
    return 'True' if r['result'] else 'False'


def main():
    runs = defaultdict(dict)  # (file, idx) -> impl -> record
    meta = {}
    for path in sorted(glob.glob(os.path.join(HERE, 'results', '*.jsonl'))):
        for line in open(path):
            r = json.loads(line)
            if 'impl' not in r:
                continue  # control-set files, not checker results
            key = (r['file'], r['idx'])
            runs[key][r['impl']] = r
            meta[key] = r
    impls = sorted({i for v in runs.values() for i in v})

    print('| input | tets | metric | ' + ' | '.join(impls) + ' | verdict |')
    print('|---|---|---|' + '---|' * len(impls) + '---|')
    for key in sorted(runs):
        m = meta[key]
        rep = V.metric_report(m['isosig'], m['sqlens'])
        if not rep['ok']:
            metric = 'invalid input'
        elif rep['flat']:
            metric = 'flat'
        else:
            metric = 'edge cond. ' + ('ok' if rep['edge_condition'] else 'FAILS')
        cells, verdicts = [], []
        for impl in impls:
            r = runs[key].get(impl)
            if r is None:
                cells.append('—')
                continue
            c = short(r) + ' (%ds)' % r['seconds']
            if r.get('result') is False and rep['ok']:
                v = V.verify_necklace(m['isosig'], m['sqlens'], r['vertex'], r['necklace'])
                tag = 'witness valid' if v['valid'] else 'witness INVALID'
                if v['valid'] and any(j.get('borderline') for j in v['junctions']):
                    tag += ' (angle = pi at a junction)'
                c += ', ' + tag
                verdicts.append('%s: %s' % (impl, tag))
            if r.get('result') is True and metric == 'flat':
                verdicts.append('%s: correct (flat)' % impl)
            if r.get('result') is False and metric == 'flat':
                verdicts.append('%s: WRONG (flat metric reported not CAT(0))' % impl)
            if r.get('error'):
                c += ': ' + r['error'][:50]
            cells.append(c)
        print('| %s %d `%s` | %d | %s | %s | %s |' % (
            key[0], key[1], m['isosig'], rep['n_tet'], metric,
            ' | '.join(cells), '; '.join(verdicts)))


if __name__ == '__main__':
    main()
