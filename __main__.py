"""



"""


from ast import literal_eval

import os

from gallery import vertex_linkcond

path = os.path.dirname(os.path.abspath(__file__)) + "/CAT0_Candidates/"

for fname in [
              'examples_1',
              'examples_2',
              'examples_3',
              'hyperbolic_candidate',
              ]:

    if fname:
        with open(path+fname+'.txt', 'r') as f: lines = f.readlines()

    if fname == 'examples_1':
        isosigs = tuple(map(str.strip, lines[1::5]))
        edgelengths_sqs = [literal_eval(l.strip()) for l in lines[2::5]]
    elif fname == 'examples_2':
    # this file has edges of length 0, leading to ZeroDivisionErrors. Replace with length 2
        edgelengths_sqs, isosigs \
        = zip(*[(literal_eval(l[l.index('['):l.index(']')+1]),
                 l[l.index(']')+2:].strip())
                for l in lines if l.startswith('candidate found')])
        edgelengths_sqs = [[2 if x==0 else x for x in y]
                            for y in edgelengths_sqs]
    elif fname == 'examples_3':
        from string import ascii_lowercase
        isosigs = [l.split()[0] for l in lines if l[0] in ascii_lowercase]
        edgelengths_sqs = [literal_eval(l) for l in lines if l[0]=='[']
    elif fname == 'hyperbolic_candidate':
        isosigs = [lines[0].split()[0]]
        edgelengths_sqs = [literal_eval(lines[1])]
    else:
        continue
    
    print(isosigs)
    print(edgelengths_sqs)
    print()
        
    
    results = []
    for isosig, edgelengths_sq in zip(isosigs, edgelengths_sqs):
        print(isosig, edgelengths_sq)
        print(vertex_linkcond(isosig, edgelengths_sq, parallel=False, progress=True, mode='findall'))#mode='cert'))
        #quit()
        #if loop:
        #    result = ('False: loop in link of vertex %s with necklace triangles %s'
        #              % (v.index(), [[t.index() for t in g[0]] for g in loop]))
        #else: result = 'True'
        #print(result)
        #results.append(result)
