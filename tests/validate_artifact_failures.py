#!/usr/bin/env python
"""Regress the reported representable binomial overflow and prism EOC failure."""
import math
import sys
import tempfile
from pathlib import Path
from validate_release import find_executable, run_query
from validate_exact import write_graph, oracle


def main():
    build = Path(sys.argv[1]).resolve()
    checks = 0
    with tempfile.TemporaryDirectory(prefix='ptcount-artifact-') as tmp:
        tmp = Path(tmp)
        data, query, output = (tmp / name for name in ('data.graph', 'query.graph', 'out.txt'))
        write_graph(data, 10001, [(0, i) for i in range(1, 10001)])
        write_graph(query, 6, [(0, i) for i in range(1, 6)])
        for program in ('piso', 'ptcount'):
            run_query(find_executable(build, program), data, query, output, 60)
            assert int(output.read_text()) == math.comb(10000, 5), program
            checks += 1
        edges = [(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3),(1,4),(2,5)]
        write_graph(data, 6, edges)
        # Both edge orbits, both root directions, and GC/VOC protect the
        # normalization convention rather than just suppressing its error.
        for mark in [(), (1,0), (2,0,3), (2,3,0), (2,0,1), (2,1,0)]:
            write_graph(query, 6, edges, mark)
            want = oracle(6, {tuple(sorted(e)) for e in edges}, 6, edges, mark)
            for program in ('piso', 'tiso', 'ptcount'):
                run_query(find_executable(build, program), data, query, output, 60)
                lines = output.read_text().splitlines()
                if not mark:
                    got = {(): int(lines[0])}; expected = {(): want.get((), 0)}
                elif mark[0] == 1:
                    got = {(i,): int(v) for i,v in enumerate(lines)}
                    expected = {(i,): want.get((i,), 0) for i in range(6)}
                else:
                    got = {(int(a),int(b)):int(c) for a,b,c in map(str.split, lines)}
                    expected = {(a,b):want.get((a,b),0) for u,v in edges for a,b in [(u,v),(v,u)]}
                assert got == expected, (program, mark, got, expected)
                checks += 1
    print(f'PASS: {checks} artifact regression checks')

if __name__ == '__main__':
    main()
