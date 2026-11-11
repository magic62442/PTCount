# PTCount

PTCount is a tree-decomposition-based framework for exact, non-induced subgraph isomorphism counting. It supports global counts, vertex-orbit counts (VOC), and edge-orbit counts (EOC). For each query, PTCount chooses between two execution methods:

- **pISO** evaluates the original pattern by constructing partial pattern matches. Its grouped plans can aggregate valid completions directly without materializing a list of complete matches, and it does not evaluate tISO correction patterns.
- **tISO** evaluates a tree-isomorphism equation over the original pattern and its correction patterns.

The optimizer combines structural rules with cardinality-based cost estimates to select the method. The release interface provides exact, non-induced ISO counting.

## Technical report

[Revisiting Exact Subgraph Counting — technical report](paper-techreport.pdf) describes the methods, analysis, and experiments, including the supplementary appendix.

## Included programs

- `ptcount`: the main executable. It plans the query and automatically selects pISO or tISO.
- `piso`: runs the pISO method. The compressed ConNode case is integrated as an internal pISO execution path rather than exposed as a separate program.
- `tiso`: runs the tISO method.

All three programs accept the same single-query interface and produce the same result-file format. Their compact standard-output summaries differ by program.

## Requirements

- CMake 3.14 or newer
- A GCC- or Clang-compatible C++17 compiler
- GLPK development files (`glpk.h` and the GLPK library)
- Python 3 (validation only)

The required nauty sources are included in `third_party/nauty`. This release has been validated on macOS arm64 with AppleClang.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executables are written to `build/bin`.

On Apple Silicon, if the installed CMake executable itself runs under Rosetta, configure a native build explicitly:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64
cmake --build build -j
```

## Input format

Data and query graphs are undirected, unlabelled simple graphs. Their text format begins with the number of vertices and undirected edges, followed by one edge per line:

```text
<number_of_vertices> <number_of_edges>
<endpoint_1> <endpoint_2>
...
```

Vertex IDs must be `0, ..., |V|-1`. Self-loops and parallel edges are not supported. Each undirected edge must appear exactly once. Query patterns must be connected and contain between 3 and 8 vertices. The upper bound follows the release's collision-free 64-bit canonical graph representation.

A data graph ends after its edge list. A query optionally appends one orbit specification:

```text
1 <query_vertex>
2 <query_endpoint_1> <query_endpoint_2>
```

The first form requests a vertex-orbit count (VOC), and the second requests an edge-orbit count (EOC). Omit the orbit specification for a global count. The two vertices in an EOC specification must be adjacent in the query.

PTCount uses non-induced semantics: a match must preserve every query edge, while additional edges among the corresponding data vertices are allowed. Counts are normalized by the applicable pattern or orbit automorphisms.

The data loader also recognizes the existing `BINGRPH` binary data-graph format. It is a host-endian compatibility format; the text format above is the documented portable interchange format. Query graphs are always text files.

## Run

Run the optimizer-selected plan:

```sh
./build/bin/ptcount -d DATA.graph -q QUERY.graph -r RESULT.txt
```

Run either method directly:

```sh
./build/bin/piso -d DATA.graph -q QUERY.graph -r RESULT.txt
./build/bin/tiso -d DATA.graph -q QUERY.graph -r RESULT.txt
```

`ptcount` and `piso` also accept a directory of `.txt` or `.graph` queries:

```sh
./build/bin/ptcount -d DATA.graph -q QUERY_DIR -b -r RESULT_DIR
./build/bin/piso -d DATA.graph -q QUERY_DIR -b -r RESULT_DIR
```

In batch mode, each result is written as `RESULT_DIR/<query-filename>`. `RESULT_DIR` must therefore differ from `QUERY_DIR`; PTCount rejects any result path that would overwrite a query or data graph. File-open and write failures produce an error.

The result file depends on the query type:

- **Global:** one integer, the number of subgraph occurrences.
- **VOC:** one integer per data vertex, in data-vertex ID order.
- **EOC:** one whitespace-separated `data_u data_v count` record for each directed data-graph arc. Endpoint roles follow the marked query edge; reverse arcs can therefore have different counts for a role-sensitive edge orbit.

If `-r` is omitted, the executables still print a compact execution summary, but do not write a count file. Its `count` is the global count or the sum of all emitted VOC/EOC entries; for same-orbit EOC this includes both directed data-graph arcs.

## Validation

The self-contained validation uses a five-clique as the data graph and checks triangle global, vertex-orbit, and edge-orbit queries. A four-node star checks the compressed pISO path integrated from ConNode, and a sparse graph checks VOC when isolated vertices outnumber directed edges. After building, run:

```sh
python tests/validate_release.py build
ctest --test-dir build --output-on-failure
```

The test runs `ptcount`, `piso`, and `tiso`, covering five counting cases per program. It verifies both result files and their aggregate summary counts. It also checks rejection of too-small, too-large, and disconnected queries, query-file overwrite attempts, and unwritable result targets (30 checks in total).

Additional CTest targets validate the production suffix-counting kernels, structural optimizer guards, and 114 exact-output comparisons against independent injective enumeration. These cover global, VOC, and EOC queries on complete and sparse graphs, including distinct edge orbits whose endpoint vertex-orbit pairs coincide.

For additional datasets and query workloads, use the text format above and the same execution commands.

## Third-party software

The bundled nauty sources retain their upstream copyright and Apache-2.0 notice in `third_party/nauty/COPYRIGHT`. GLPK is detected as an external system dependency and is not redistributed here.
