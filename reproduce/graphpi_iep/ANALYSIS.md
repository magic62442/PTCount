# GraphPi IEP: count discrepancies

## Small-graph verification

Six K3,3 host labelings were checked with the existing local GraphPi executable, both with IEP enabled and disabled. Exhaustive enumeration gives 72 injective mappings and one distinct matched subgraph for each host. With IEP disabled, all six counts are 1; with IEP enabled, two counts are 1 and four are 0. These outputs match GraphSet on the same six labelings. The input edge order preserves node IDs under GraphPi's first-seen-node loader. The one-versus-two pair is verified separately in the GraphSet report.

## Large graph verification

The comparison using the lastfm data graph covers six-node queries 70, 91, 92 and 106, shared with the GraphSet report. The data graph has 7,624 nodes and 27,806 undirected edges.

## Observations

| Six-node query ID | IEP on | IEP off | Reference count |
|---|---:|---:|---:|
| 70 | 8,275,161 | 8,119,142 | 8,119,142 |
| 91 | 52,535,639 | 47,273,619 | 47,273,619 |
| 92 | 73,213,475 | 48,826,660 | 48,826,660 |
| 106 | 50,105,792 | 37,417,280 | 37,417,280 |

## Queries

The pattern edge lists retain their original query IDs: [70](queries/70.txt) (K3,3), [91](queries/91.txt), [92](queries/92.txt), and [106](queries/106.txt).
