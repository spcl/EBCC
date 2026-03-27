#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Sequence


DEFAULT_RELATIVE_ERROR_TARGETS = (0.001, 0.01, 0.1)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def format_number(value: float) -> str:
    return format(value, "g")


def make_target_label(value: float) -> str:
    return f"rel{format_number(value).replace('.', 'p')}"


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(
        description=(
            "Run perf profiling for EBCC CDO compression across a sweep of "
            "relative error targets."
        )
    )
    parser.add_argument(
        "--plugin-path",
        type=Path,
        required=True,
        help="Path exported as HDF5_PLUGIN_PATH for the CDO command.",
    )
    parser.add_argument(
        "--input-data",
        type=Path,
        required=True,
        help="Input NetCDF file passed to cdo copy.",
    )
    parser.add_argument(
        "--output-path",
        type=Path,
        required=True,
        help="Directory where perf data files and compressed outputs are written.",
    )
    parser.add_argument(
        "--python",
        type=Path,
        default=root / ".venv" / "bin" / "python",
        help="Python interpreter used to run filter_wrapper.py.",
    )
    parser.add_argument(
        "--filter-wrapper",
        type=Path,
        default=root / "ebcc" / "filter_wrapper.py",
        help="Path to ebcc/filter_wrapper.py.",
    )
    parser.add_argument(
        "--perf",
        default="perf",
        help="Path to the perf executable.",
    )
    parser.add_argument(
        "--cdo",
        default="cdo",
        help="Path to the cdo executable.",
    )
    parser.add_argument(
        "--base-cr",
        type=float,
        default=200.0,
        help="Base compression ratio passed to filter_wrapper.py.",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=721,
        help="Tile height passed to filter_wrapper.py.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1440,
        help="Tile width passed to filter_wrapper.py.",
    )
    parser.add_argument(
        "--relative-error-targets",
        nargs="+",
        type=float,
        default=list(DEFAULT_RELATIVE_ERROR_TARGETS),
        help=(
            "Relative error targets to sweep. Defaults to 0.001 0.01 0.1 "
            "(0.1%%, 1%%, 10%%)."
        ),
    )
    parser.add_argument(
        "--call-graph",
        default="dwarf",
        help="perf record --call-graph mode.",
    )
    parser.add_argument(
        "--mmap-pages",
        default="16M",
        help="perf record --mmap-pages value.",
    )
    parser.add_argument(
        "--event",
        default=None,
        help="Optional perf event passed to perf record -e.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing perf data and compressed outputs.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing perf record.",
    )
    return parser.parse_args(list(argv))


def ensure_inputs(args: argparse.Namespace) -> None:
    if not args.plugin_path.exists():
        raise SystemExit(f"Plugin path does not exist: {args.plugin_path}")
    if not args.input_data.is_file():
        raise SystemExit(f"Input data does not exist: {args.input_data}")
    if not args.python.is_file():
        raise SystemExit(f"Python executable does not exist: {args.python}")
    if not args.filter_wrapper.is_file():
        raise SystemExit(f"filter_wrapper.py does not exist: {args.filter_wrapper}")
    args.output_path.mkdir(parents=True, exist_ok=True)


def run_filter_wrapper(args: argparse.Namespace, relative_error_target: float) -> str:
    command = [
        str(args.python),
        str(args.filter_wrapper),
        "--base_cr",
        format_number(args.base_cr),
        "--height",
        str(args.height),
        "--width",
        str(args.width),
        "--relative_error_target",
        format_number(relative_error_target),
    ]
    completed = subprocess.run(
        command,
        cwd=repo_root(),
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        if completed.stdout:
            sys.stderr.write(completed.stdout)
        if completed.stderr:
            sys.stderr.write(completed.stderr)
        raise SystemExit(completed.returncode)

    filter_opts = completed.stdout.strip().splitlines()
    if not filter_opts:
        raise SystemExit(
            f"filter_wrapper.py produced no filter options for target {relative_error_target}"
        )
    return filter_opts[-1].strip()


def run_command(command: Sequence[str], env: dict[str, str], dry_run: bool) -> None:
    print(shlex.join(command), file=sys.stderr)
    if dry_run:
        return

    completed = subprocess.run(command, env=env, check=False)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


def build_perf_command(
    args: argparse.Namespace,
    perf_data_path: Path,
    filter_opts: str,
    output_nc_path: Path,
) -> list[str]:
    command = [
        args.perf,
        "record",
        "--call-graph",
        args.call_graph,
        "--mmap-pages",
        args.mmap_pages,
        "--stat",
    ]
    if args.event:
        command.extend(["-e", args.event])
    command.extend(
        [
            "-o",
            str(perf_data_path),
            "--",
            args.cdo,
            "-b",
            "F32",
            "--filter",
            filter_opts,
            "copy",
            str(args.input_data),
            str(output_nc_path),
        ]
    )
    return command


def maybe_check_overwrite(paths: Sequence[Path], force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        joined = "\n".join(str(path) for path in existing)
        raise SystemExit(
            "Refusing to overwrite existing outputs without --force:\n"
            f"{joined}"
        )


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    ensure_inputs(args)

    env = os.environ.copy()
    env["HDF5_PLUGIN_PATH"] = str(args.plugin_path.resolve())

    input_stem = args.input_data.stem
    for relative_error_target in args.relative_error_targets:
        label = make_target_label(relative_error_target)
        perf_data_path = args.output_path / f"perf_cpu_{input_stem}_{label}.data"
        output_nc_path = args.output_path / f"{input_stem}_ebcc_perf_{label}.nc"
        maybe_check_overwrite((perf_data_path, output_nc_path), args.force)

        filter_opts = run_filter_wrapper(args, relative_error_target)
        print(
            f"Running profiling sweep for relative_error_target={format_number(relative_error_target)}",
            file=sys.stderr,
        )
        command = build_perf_command(args, perf_data_path, filter_opts, output_nc_path)
        run_command(command, env=env, dry_run=args.dry_run)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
