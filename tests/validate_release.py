#!/usr/bin/env python
"""Run PTCount's self-contained global/VOC/EOC release checks."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TESTS = Path(__file__).resolve().parent
PROGRAMS = ("ptcount", "piso", "tiso")


def find_executable(build_dir: Path, program: str) -> Path:
    candidates = (build_dir / "bin" / program, build_dir / program)
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    checked = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"cannot find {program}; checked: {checked}")


def run_query(
    executable: Path,
    data: Path,
    query: Path,
    output: Path,
    timeout: float,
) -> str:
    command = [
        str(executable),
        "-d",
        str(data),
        "-q",
        str(query),
        "-r",
        str(output),
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )
    if completed.returncode != 0:
        detail = (
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
        raise RuntimeError(detail)
    if not output.is_file():
        raise RuntimeError(
            f"{executable.name} returned success but did not create {output}"
        )
    return completed.stdout


def expect_rejection(
    executable: Path,
    data: Path,
    query: Path,
    expected_error: str,
    timeout: float,
) -> None:
    command = [str(executable), "-d", str(data), "-q", str(query)]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )
    if completed.returncode == 0:
        raise AssertionError(
            f"unsupported query was accepted: {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}"
        )
    if completed.returncode < 0:
        raise AssertionError(
            f"unsupported query terminated by signal {-completed.returncode}: "
            f"{' '.join(command)}"
        )
    if expected_error not in completed.stderr:
        raise AssertionError(
            f"rejection should mention {expected_error!r}; stderr was: "
            f"{completed.stderr!r}"
        )


def expect_result_rejection(
    executable: Path,
    data: Path,
    query: Path,
    output: Path,
    expected_error: str,
    timeout: float,
) -> None:
    command = [
        str(executable),
        "-d",
        str(data),
        "-q",
        str(query),
        "-r",
        str(output),
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )
    if completed.returncode == 0:
        raise AssertionError(
            f"unsafe result path was accepted: {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}"
        )
    if completed.returncode < 0:
        raise AssertionError(
            f"unsafe result path terminated by signal {-completed.returncode}: "
            f"{' '.join(command)}"
        )
    if expected_error not in completed.stderr:
        raise AssertionError(
            f"rejection should mention {expected_error!r}; stderr was: "
            f"{completed.stderr!r}"
        )


def summary_count(program: str, stdout: str) -> int:
    lines = [line for line in stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise AssertionError(f"{program} did not print a result row: {stdout!r}")
    fields = lines[-1].split("\t")
    index = {"ptcount": 3, "piso": 2, "tiso": 1}[program]
    try:
        return int(fields[index])
    except (IndexError, ValueError) as error:
        raise AssertionError(
            f"cannot parse aggregate count from {program} output: {lines[-1]!r}"
        ) from error


def parse_ints(path: Path) -> list[int]:
    try:
        return [int(token) for token in path.read_text().split()]
    except ValueError as error:
        raise AssertionError(f"non-integer result in {path}: {error}") from error


def check_global(path: Path) -> None:
    values = parse_ints(path)
    if values != [10]:
        raise AssertionError(f"global result should be [10], got {values}")


def check_star_global(path: Path) -> None:
    values = parse_ints(path)
    if values != [20]:
        raise AssertionError(f"star global result should be [20], got {values}")


def check_voc(path: Path) -> None:
    values = parse_ints(path)
    expected = [6] * 5
    if values != expected:
        raise AssertionError(f"VOC result should be {expected}, got {values}")


def check_eoc(path: Path) -> None:
    records: dict[tuple[int, int], int] = {}
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 3:
            raise AssertionError(
                f"EOC line {line_number} should have three fields, got: {line!r}"
            )
        try:
            source, target, count = map(int, fields)
        except ValueError as error:
            raise AssertionError(
                f"EOC line {line_number} contains a non-integer: {line!r}"
            ) from error
        key = (source, target)
        if key in records:
            raise AssertionError(f"duplicate EOC record for directed edge {key}")
        records[key] = count

    expected_edges = {(u, v) for u in range(5) for v in range(5) if u != v}
    if set(records) != expected_edges:
        missing = sorted(expected_edges - set(records))
        extra = sorted(set(records) - expected_edges)
        raise AssertionError(f"unexpected EOC edge set; missing={missing}, extra={extra}")
    wrong = {edge: count for edge, count in records.items() if count != 3}
    if wrong:
        raise AssertionError(f"each K5 edge should have EOC 3, got {wrong}")


def check_sparse_voc(path: Path) -> None:
    values = parse_ints(path)
    expected = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
    if values != expected:
        raise AssertionError(f"sparse VOC result should be {expected}, got {values}")


def validate_program(executable: Path, output_dir: Path, timeout: float) -> None:
    cases = (
        ("global", TESTS / "data_k5.graph", TESTS / "query_triangle_global.graph", check_global, 10),
        ("compressed-global", TESTS / "data_k5.graph", TESTS / "query_star_global.graph", check_star_global, 20),
        ("voc", TESTS / "data_k5.graph", TESTS / "query_triangle_voc.graph", check_voc, 30),
        ("eoc", TESTS / "data_k5.graph", TESTS / "query_triangle_eoc.graph", check_eoc, 60),
        ("sparse-voc", TESTS / "data_sparse_triangle.graph", TESTS / "query_triangle_voc.graph", check_sparse_voc, 3),
    )
    for name, data, query, checker, expected_count in cases:
        output = output_dir / f"{executable.name}_{name}.txt"
        stdout = run_query(executable, data, query, output, timeout)
        checker(output)
        actual_count = summary_count(executable.name, stdout)
        if actual_count != expected_count:
            raise AssertionError(
                f"{executable.name} {name} summary should be {expected_count}, "
                f"got {actual_count}"
            )
        if executable.name == "ptcount" and name == "compressed-global":
            fields = [line for line in stdout.splitlines() if line.strip()][-1].split("\t")
            if len(fields) < 3 or fields[2] != "compressed-piso":
                raise AssertionError(
                    "ptcount star query should exercise the compressed-piso "
                    f"implementation, got: {stdout!r}"
                )
        print(f"PASS {executable.name}: {name}")

    invalid_cases = (
        ("too-small", "2 1\n0 1\n", "between 3 and 8"),
        (
            "too-large",
            "9 8\n0 1\n1 2\n2 3\n3 4\n4 5\n5 6\n6 7\n7 8\n",
            "between 3 and 8",
        ),
        ("disconnected", "3 1\n0 1\n", "must be connected"),
    )
    for name, contents, expected_error in invalid_cases:
        query = output_dir / f"invalid_{name}.graph"
        query.write_text(contents)
        expect_rejection(
            executable,
            TESTS / "data_k5.graph",
            query,
            expected_error,
            timeout,
        )
        print(f"PASS {executable.name}: reject-{name}")

    overwrite_query = output_dir / f"{executable.name}_overwrite_query.graph"
    shutil.copyfile(TESTS / "query_triangle_global.graph", overwrite_query)
    original_query = overwrite_query.read_bytes()
    expect_result_rejection(
        executable,
        TESTS / "data_k5.graph",
        overwrite_query,
        overwrite_query,
        "result path must differ from query path",
        timeout,
    )
    if overwrite_query.read_bytes() != original_query:
        raise AssertionError(f"{executable.name} modified its query on rejection")
    print(f"PASS {executable.name}: reject-query-overwrite")

    unwritable_target = output_dir / f"{executable.name}_result_directory"
    unwritable_target.mkdir()
    expect_result_rejection(
        executable,
        TESTS / "data_k5.graph",
        TESTS / "query_triangle_global.graph",
        unwritable_target,
        "cannot open result path",
        timeout,
    )
    print(f"PASS {executable.name}: reject-unwritable-result")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "build_dir",
        nargs="?",
        type=Path,
        default=ROOT / "build",
        help="CMake build directory (default: %(default)s)",
    )
    parser.add_argument(
        "--program",
        action="append",
        choices=PROGRAMS,
        help="validate only this program; repeat to select several",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=120.0,
        help="timeout in seconds for each query (default: %(default)s)",
    )
    arguments = parser.parse_args()

    programs = tuple(arguments.program or PROGRAMS)
    build_dir = arguments.build_dir.expanduser().resolve()
    try:
        executables = [find_executable(build_dir, program) for program in programs]
        with tempfile.TemporaryDirectory(prefix="ptcount-release-") as temporary:
            output_dir = Path(temporary)
            for executable in executables:
                validate_program(executable, output_dir, arguments.timeout)
    except (FileNotFoundError, RuntimeError, AssertionError, subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1

    print(f"PASS: validated {', '.join(programs)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
