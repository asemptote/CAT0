"""Load the candidate (isosig, squared edge lengths) pairs from CAT0_Candidates/.

The parsing mirrors __main__.py (and the 2020 vLinkCond.py) exactly, including
the replacement of zero squared lengths by 2 in examples_2.txt, so that every
implementation is run on the same inputs that produced the historical results.
"""

import os
from ast import literal_eval
from string import ascii_lowercase

CANDIDATES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          'CAT0_Candidates')
FILES = ('examples_1', 'examples_2', 'examples_3', 'hyperbolic_candidate')


# examples_3.txt holds edge *labels* from search_in_Pachner_graph__symmetric_tet_shape.py
# (2 = the pair of opposite edges of each tetrahedron with dihedral angle pi/2,
# 1 = the other four edges, dihedral angle pi/3), not squared lengths. Read as
# squared lengths every tetrahedron is degenerate. The Sommerville disphenoid with
# these angles has edge lengths 2 and sqrt(3), i.e. squared lengths 4 and 3.
DISPHENOID = {1: 3, 2: 4}


def load(fname):
    """Returns a list of dicts with keys file, idx, isosig, sqlens, raw_sqlens."""
    if fname == 'examples_3_disphenoid':
        out = load('examples_3')
        for e in out:
            e['file'] = fname
            e['sqlens'] = [DISPHENOID[x] for x in e['raw_sqlens']]
        return out
    with open(os.path.join(CANDIDATES, fname + '.txt')) as f:
        lines = f.readlines()
    if fname == 'examples_1':
        isosigs = [l.strip() for l in lines[1::5]]
        raw = [literal_eval(l.strip()) for l in lines[2::5]]
    elif fname == 'examples_2':
        pairs = [(literal_eval(l[l.index('['):l.index(']')+1]),
                  l[l.index(']')+2:].strip())
                 for l in lines if l.startswith('candidate found')]
        raw, isosigs = zip(*pairs)
    elif fname == 'examples_3':
        isosigs = [l.split()[0] for l in lines if l[0] in ascii_lowercase]
        raw = [literal_eval(l) for l in lines if l[0] == '[']
    elif fname == 'hyperbolic_candidate':
        isosigs = [lines[0].split()[0]]
        raw = [literal_eval(lines[1])]
    else:
        raise ValueError(fname)
    return [dict(file=fname, idx=i, isosig=s, raw_sqlens=list(r),
                 sqlens=[2 if x == 0 else x for x in r])
            for i, (s, r) in enumerate(zip(isosigs, raw))]


def load_all(fnames=FILES):
    return [e for f in fnames for e in load(f)]
