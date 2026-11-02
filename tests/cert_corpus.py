"""Certificates for cross-checking the Lean checker (lean/: vlinkcheck) against
the Python kernel.

    python tests/cert_corpus.py OUT.jsonl [--workers N]   # build the corpus
    python tests/cert_corpus.py --compare OUT.jsonl LEAN_OUTPUT.txt

The corpus has one line per certificate: {certificate, label, python_ok,
python_msg}. It holds genuine False certificates (every False metric in the
regression data and walks), wrong ones (the 2020 code's witnesses, the
junction-bug witnesses, corrupted certificates), and True certificates for
the regression metrics answered True, in both pi settings, with corrupted
variants. --compare checks that vlinkcheck (one line per certificate:
ACCEPT/REJECT ...) reaches the Python kernel's verdict on every one.
"""

import copy
import json
import os
import sys
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation'), os.path.join(ROOT, 'tests')]

from vlink import kernel  # noqa: E402
from vlink.certify import true_certificate  # noqa: E402
from vlink.check import CERT_FORMAT, certificate  # noqa: E402

RESULTS = os.path.join(ROOT, 'validation', 'results')
DATA = json.load(open(os.path.join(ROOT, 'tests', 'data', 'regression_cases.json')))
FLAT = ('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])


def _false_cert(args):
    sig, sq = args
    return certificate(sig, list(sq))['certificate']


def _true_certs(args):
    sig, sq = args
    return [true_certificate(sig, list(sq), pi_beads=pi) for pi in (True, False)]


def build(workers):
    import test_certificates as TC
    import verify_cert as V
    out = []
    with Pool(workers) as pool:
        falses = TC.false_metrics()
        for c in pool.map(_false_cert, falses):
            out.append(('false genuine', c))
        trues = sorted({(c['isosig'], tuple(c['sqlens'])) for c in DATA['cases']
                        if c['expected'] is True})
        for cs in pool.map(_true_certs, trues):
            out += [('true genuine', c) for c in cs]
    # wrong False certificates
    for path in ('legacy2020.jsonl', 'legacy2020.disphenoid.jsonl'):
        for line in open(os.path.join(RESULTS, path)):
            r = json.loads(line)
            if r.get('result') is not False or r['file'] == 'examples_3':
                continue
            if any(len(b['ts']) > 1 for b in r['necklace']):
                continue
            beads = V.regina_to_canonical(r['isosig'], r['vertex'], r['necklace'])
            out.append(('false legacy2020', TC._thin_cert(r['isosig'], r['sqlens'], r['vertex'], beads)))
    for line in open(os.path.join(RESULTS, 'vlink_repeats.buggy.jsonl')):
        r = json.loads(line)
        beads = V.regina_to_canonical(r['isosig'], r['vertex'], r['necklace'])
        out.append(('false junction-bug', TC._thin_cert(r['isosig'], r['sqlens'], r['vertex'], beads)))
    for _, good in [o for o in out if o[0] == 'false genuine' and o[1]['reason'] == 'vertex'][:5]:
        c = copy.deepcopy(good); c['necklace'] = c['necklace'][1:]; out.append(('false corrupt', c))
        c = copy.deepcopy(good); b = c['necklace'][0]; b['start'], b['end'] = b['end'], b['start']
        out.append(('false corrupt', c))
        c = copy.deepcopy(good); c['necklace'] = c['necklace'] + c['necklace'][:1]; out.append(('false corrupt', c))
        c = copy.deepcopy(good); c['vertex'] = 99; out.append(('false corrupt', c))
        c = copy.deepcopy(good); c['sqlens'] = list(FLAT[1]); c['isosig'] = FLAT[0]
        out.append(('false corrupt', c))
    for reason, key in (('edge', 'edge'), ('invalid', 'tet')):
        out.append(('false corrupt', dict(format=CERT_FORMAT, isosig=FLAT[0], sqlens=FLAT[1],
                                          reason=reason, **{key: 0})))
    out.append(('false genuine', certificate(FLAT[0], [1] * 7)['certificate']))
    out.append(('false genuine', certificate('gLLPQbefefefhhxxqxh', [1, 2, 1, 2, 1, 1, 2])['certificate']))
    # wrong True certificates
    good = true_certificate(*FLAT)
    s, tree = next((s, t) for s, t in good['vertices'][0]['trees'].items() if t)
    c = copy.deepcopy(good); del c['vertices'][0]['trees'][s][next(iter(tree))]; out.append(('true corrupt', c))
    c = copy.deepcopy(good); t = c['vertices'][0]['trees'][s]; x = next(iter(t)); t[x] = int(x)
    out.append(('true corrupt', c))
    c = copy.deepcopy(good); c['vertices'][0]['n_beads'] += 1; out.append(('true corrupt', c))
    c = copy.deepcopy(good); c['vertices'] = []; out.append(('true corrupt', c))
    c = copy.deepcopy(good); c['options']['pi_beads'] = False; out.append(('true corrupt', c))
    c = copy.deepcopy(good); c['options']['allow_repeats'] = False; out.append(('true corrupt', c))
    c = copy.deepcopy(good); c['sqlens'] = [1] * 7; out.append(('true corrupt', c))
    out.append(('true not-cat0', true_certificate('oLLvzLLQQQcbehlmmklnknmnklxxahaqhqxhaqaqq',
                                                  [1, 2, 1, 3, 2, 1, 3, 1, 2, 1, 2, 2, 1, 2, 2])))
    return out


def _verdict(lc):
    label, c = lc
    if c.get('format') == kernel.TRUE_FORMAT:
        ok, msg = kernel.check_true_certificate(c, exact_only=True)
    else:
        ok, msg = kernel.check_certificate(c)
    return dict(label=label, python_ok=ok, python_msg=msg, certificate=c)


def main():
    if sys.argv[1] == '--compare':
        corpus = [json.loads(line) for line in open(sys.argv[2])]
        lean = [line.split(' ', 2) for line in open(sys.argv[3]) if line.strip()]
        assert len(corpus) == len(lean), (len(corpus), len(lean))
        bad = 0
        counts = {}
        for r, (verdict, claim, msg) in zip(corpus, lean):
            ok = verdict == 'ACCEPT'
            key = (r['label'], r['python_ok'])
            counts[key] = counts.get(key, 0) + 1
            if ok != r['python_ok']:
                bad += 1
                print('DISAGREE', r['label'], 'python', r['python_ok'], r['python_msg'][:80],
                      '| lean', verdict, msg.strip()[:80])
        for (label, ok), n in sorted(counts.items()):
            print('%-20s python %-5s %4d' % (label, ok, n))
        print('AGREE' if not bad else 'DISAGREEMENTS: %d' % bad)
        return
    workers = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 4
    items = build(workers)
    with Pool(workers) as pool:
        rows = pool.map(_verdict, items)
    with open(sys.argv[1], 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')
    print('corpus:', len(rows), 'certificates;',
          sum(r['python_ok'] for r in rows), 'accepted by the Python kernel')


if __name__ == '__main__':
    main()
