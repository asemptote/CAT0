"""certcheck must accept every certificate vlink produces for a False answer
and reject wrong ones.

    python tests/test_certificates.py
"""

import copy
import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation')]

from certcheck import check_certificate  # noqa: E402
from vlink.check import CERT_FORMAT, certificate  # noqa: E402

RESULTS = os.path.join(ROOT, 'validation', 'results')
FLAT = ('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])


def false_metrics():
    """Metrics vlink answered False (validation inputs and walks)."""
    seen = set()
    data = json.load(open(os.path.join(ROOT, 'tests', 'data', 'regression_cases.json')))
    for c in data['cases']:
        if c['expected'] is False:
            seen.add((c['isosig'], tuple(c['sqlens'])))
    for path in glob.glob(os.path.join(RESULTS, 'census_walks_vlink.*.jsonl')):
        for line in open(path):
            r = json.loads(line)
            if r.get('result') is False:
                seen.add((r['isosig'], tuple(r['sqlens'])))
    return sorted(seen)


def test_vlink_certificates_accepted():
    n = 0
    for sig, sq in false_metrics():
        r = certificate(sig, list(sq))
        assert r['result'] is False, sig
        ok, msg = check_certificate(r['certificate'])
        assert ok, (sig, msg)
        n += 1
    # an edge failure (regular tetrahedra) and an invalid input
    r = certificate(FLAT[0], [1] * 7)
    assert r['reason'] == 'edge' and check_certificate(r['certificate'])[0]
    r = certificate('gLLPQbefefefhhxxqxh', [1, 2, 1, 2, 1, 1, 2])   # examples_3 labels
    assert r['reason'] == 'invalid' and check_certificate(r['certificate'])[0]
    return n


def _thin_cert(sig, sq, vertex, beads):
    return dict(format=CERT_FORMAT, isosig=sig, sqlens=sq, reason='vertex', vertex=vertex,
                necklace=[dict(ts=b['ts'], start=b['start'], end=b['end'], steps=[],
                               pi_point=None) for b in beads])


def test_wrong_certificates_rejected():
    rejected = 0
    # 1. the September 2020 code's witnesses (all shown invalid by verify_cert);
    #    the one-triangle-bead ones are complete certificates as they stand
    for path in ('legacy2020.jsonl', 'legacy2020.disphenoid.jsonl'):
        for line in open(os.path.join(RESULTS, path)):
            r = json.loads(line)
            if r.get('result') is not False or r['file'] == 'examples_3':
                continue
            if any(len(b['ts']) > 1 for b in r['necklace']):
                continue
            ok, _ = check_certificate(_thin_cert(r['isosig'], r['sqlens'], r['vertex'], r['necklace']))
            assert not ok, (r['file'], r['idx'])
            rejected += 1
    # 2. the 'there and back' loop accepted by the junction bug fixed earlier
    for line in open(os.path.join(RESULTS, 'vlink_repeats.buggy.jsonl')):
        r = json.loads(line)
        ok, msg = check_certificate(_thin_cert(r['isosig'], r['sqlens'], r['vertex'], r['necklace']))
        assert not ok, (r['file'], r['idx'])
        rejected += 1
    # 3. corrupted genuine certificates
    sig, sq = false_metrics()[0]
    good = certificate(sig, list(sq))['certificate']
    assert check_certificate(good)[0]
    bad = []
    c = copy.deepcopy(good); c['necklace'] = c['necklace'][1:]; bad.append(c)
    c = copy.deepcopy(good); b = c['necklace'][0]; b['start'], b['end'] = b['end'], b['start']; bad.append(c)
    c = copy.deepcopy(good); c['necklace'] = c['necklace'] + c['necklace'][:1]; bad.append(c)
    c = copy.deepcopy(good); c['isosig'], c['sqlens'] = FLAT[0], FLAT[1]; bad.append(c)
    c = copy.deepcopy(good); c['vertex'] = 99; bad.append(c)
    bad.append(dict(format=CERT_FORMAT, isosig=FLAT[0], sqlens=FLAT[1], reason='edge', edge=0))
    bad.append(dict(format=CERT_FORMAT, isosig=FLAT[0], sqlens=FLAT[1], reason='invalid', tet=0))
    bad.append(dict(format='other', isosig=FLAT[0], sqlens=FLAT[1], reason='edge', edge=0))
    for c in bad:
        ok, _ = check_certificate(c)
        assert not ok, c.get('reason')
        rejected += 1
    return rejected


if __name__ == '__main__':
    print('ok, accepted vlink certificates:', test_vlink_certificates_accepted())
    print('ok, rejected wrong certificates:', test_wrong_certificates_rejected())
