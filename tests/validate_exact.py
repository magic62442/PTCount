#!/usr/bin/env python
"""Independent injective enumeration for post-release ISO regressions.
No network, external data, NetworkX, or historical result oracle required.
"""
import itertools,random,subprocess,sys,tempfile
from pathlib import Path
from collections import Counter
from validate_release import find_executable,run_query
ROOT=Path(__file__).resolve().parents[1]
def write_graph(path,n,edges,mark=()):
    path.write_text(f'{n} {len(edges)}\n'+''.join(f'{u} {v}\n' for u,v in edges)+(' '.join(map(str,mark))+'\n' if mark else ''))
def read_query(path):
    lines=path.read_text().splitlines();n,m=map(int,lines[0].split())
    edges=[tuple(map(int,x.split())) for x in lines[1:m+1]]
    mark=list(map(int,lines[m+1].split())) if len(lines)>m+1 else []
    return n,edges,mark

def oracle(n,E,q,edges,mark):
    # Divide fixed-root mappings by the automorphisms fixing those root roles.
    qe={tuple(sorted(e)) for e in edges}; roots=mark[1:] if mark else []
    aut=sum(all(p[u]==u for u in roots) and all(tuple(sorted((p[u],p[v]))) in qe for u,v in edges)
            for p in itertools.permutations(range(q)))
    counts=Counter()
    for p in itertools.permutations(range(n),q):
        if all(tuple(sorted((p[u],p[v]))) in E for u,v in edges):
            counts[tuple(p[u] for u in roots)]+=1
    assert all(x%aut==0 for x in counts.values())
    return {k:v//aut for k,v in counts.items()}

def main():
    build=Path(sys.argv[1]).resolve();checks=0;failures=[]
    with tempfile.TemporaryDirectory(prefix='ptcount-exact-') as td:
        td=Path(td);queries=[]
        patterns=[('cycle5',5,[(i,(i+1)%5) for i in range(5)]),
                  ('book4',6,[(u,v) for u in range(4) for v in (4,5)]+[(4,5)]),
                  ('path6',6,[(i,i+1) for i in range(5)]),
                  ('nested5',5,[(0,4),(1,3),(1,4),(2,3),(2,4),(3,4)])]
        for name,q,edges in patterns:
            for mode,mark in [('gc',()),('voc',(1,0)),('eoc',(2,*edges[0]))]:
                path=td/f'{name}_{mode}.graph';write_graph(path,q,edges,mark);queries.append(path)
        queries+=sorted((ROOT/'tests/queries').glob('eoc6_*.graph'))
        if len(sys.argv)>2:
            queries=sorted(Path(sys.argv[2]).glob('*.txt'),key=lambda p:int(p.stem))
            if not queries: raise ValueError('external query directory is empty')
        rng=random.Random(20261106);n=7
        complete=list(itertools.combinations(range(n),2))
        sparse=[e for e in complete if rng.random()<0.63]
        for name,edges in [('complete',complete),('sparse',sparse)]:
            data=td/f'{name}.graph';write_graph(data,n,edges);E=set(edges)
            for query in queries:
                q,qe,mark=read_query(query);want=oracle(n,E,q,qe,mark)
                for prog in ('ptcount','piso','tiso'):
                    output=td/'output.txt'
                    try: stdout=run_query(find_executable(build,prog),data,query,output,60)
                    except Exception as e:
                        failures.append((name,query.name,prog,str(e)[-250:]));continue
                    lines=output.read_text().splitlines()
                    if not mark:got={():int(lines[0])};expected={():want.get((),0)}
                    elif mark[0]==1:
                        assert len(lines)==n
                        got={(i,):int(x) for i,x in enumerate(lines)};expected={(i,):want.get((i,),0) for i in range(n)}
                    else:
                        got={(int(a),int(b)):int(c) for a,b,c in map(str.split,lines)}
                        expected={(a,b):want.get((a,b),0) for u,v in edges for a,b in [(u,v),(v,u)]}
                    if got!=expected:
                        mismatch=next((k,got.get(k),expected.get(k)) for k in got.keys()|expected.keys() if got.get(k)!=expected.get(k))
                        failures.append((name,query.name,prog,mismatch))
                    if prog=='ptcount' and query.name.startswith('cycle5'):
                        header,row=stdout.strip().splitlines()[-2:]
                        record=dict(zip(header.split('\t'),row.split('\t')))
                        if record['method']!='tiso': failures.append((name,query.name,prog,'cycle route is not tISO'))
                    checks+=1
    if failures:
        for f in failures: print('FAIL:',f)
        raise SystemExit(f'{len(failures)} failures across {checks} completed runs')
    print(f'PASS: {checks} independent GC/VOC/EOC exact-output checks, across {len(queries)} queries')
if __name__=='__main__':main()
