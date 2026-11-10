# GraphSet IEP: count discrepancies

## Small graph verification

Unmodified upstream IEP-on code returns 1 and 2 for two isomorphic six-node hosts, each with reference count 1. Across all ten labeled K3,3 edge sets, IEP-on has seven mismatches; an IEP-off driver returns 1 on all ten.

## Large graph verification

The comparison using the lastfm data graph covers six-node queries 70, 91, 92 and 106, shared with the GraphPi report.

| Query | GraphSet IEP on | Verified reference | Historical comparison |
|---|---:|---:|---|
| 70 | 8275161 | 8119142 | mismatch |
| 91 | 19690294 | 47273619 | mismatch |
| 92 | 4855834 | 48826660 | mismatch |
| 106 | 50105792 | 37417280 | mismatch |

## Queries

The pattern edge lists retain their original query IDs: [70](queries/70.txt) (K3,3), [91](queries/91.txt), [92](queries/92.txt), and [106](queries/106.txt).
