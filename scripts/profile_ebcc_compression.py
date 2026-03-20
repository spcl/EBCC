#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


HEADER_RE = re.compile(r":\s*(?:(?P<period>\d+)\s+)?(?P<event>.+):\s*$")
HEX_RE = re.compile(r"^(?:0x)?[0-9a-fA-F]+$")

CATEGORY_PREFIXES = {
    "base_encode": ("j2k_encode_internal",),
    "base_decode": ("j2k_decode_internal",),
    "residual_encode": (
        "spiht_encode",
        "ZSTD_compress",
        "ZSTD_compress2",
        "ZSTD_compressCCtx",
        "ZSTD_compressStream",
        "ZSTD_compressStream2",
    ),
    "residual_decode": (
        "spiht_decode",
        "ZSTD_decompress",
        "ZSTD_decompressDCtx",
        "ZSTD_decompressStream",
    ),
}

CATEGORY_ORDER = (
    "base_encode",
    "base_decode",
    "residual_encode",
    "residual_decode",
)

PHASE_DEFINITIONS = (
    ("base", ("base_encode", "base_decode")),
    ("residual", ("residual_encode", "residual_decode")),
)


@dataclass(frozen=True)
class AnalysisResult:
    total_weight: float
    classified_weight: float
    category_weights: Counter[str]
    processed_blocks: int
    weighted_by_period: bool
    unclassified_leaves: Counter[str]


def parse_cli(argv: Sequence[str]) -> tuple[argparse.Namespace, list[str]]:
    perf_command: list[str] = []
    if "--" in argv:
        split = argv.index("--")
        perf_command = list(argv[split + 1 :])
        argv = argv[:split]

    parser = argparse.ArgumentParser(
        description=(
            "Profile EBCC compression time from perf samples. "
            "Either analyze an existing perf.data file or record one first."
        )
    )
    parser.add_argument(
        "perf_data",
        nargs="?",
        default="perf_cpu.data",
        help="Path to an existing perf.data file, or the output path when recording.",
    )
    parser.add_argument(
        "--perf",
        default="perf",
        help="Path to the perf executable (default: %(default)s).",
    )
    parser.add_argument(
        "--call-graph",
        default="dwarf",
        help="perf record --call-graph mode (default: %(default)s).",
    )
    parser.add_argument(
        "--mmap-pages",
        default="64M",
        help="perf record --mmap-pages value (default: %(default)s).",
    )
    parser.add_argument(
        "--event",
        default=None,
        help="Optional perf event selector passed to perf record -e.",
    )
    parser.add_argument(
        "--show-unclassified",
        type=int,
        default=10,
        help="How many unclassified leaf symbols to print (default: %(default)s).",
    )

    args = parser.parse_args(list(argv))
    return args, perf_command


def run_command(command: Sequence[str]) -> None:
    completed = subprocess.run(command)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


def record_perf(
    perf_exe: str,
    output_path: Path,
    command: Sequence[str],
    call_graph: str,
    mmap_pages: str,
    event: str | None,
) -> None:
    perf_cmd = [
        perf_exe,
        "record",
        "--call-graph",
        call_graph,
        "--mmap-pages",
        mmap_pages,
        "--stat",
    ]
    if event:
        perf_cmd.extend(["-e", event])
    perf_cmd.extend(["-o", str(output_path), "--", *command])

    print(f"Recording perf data to {output_path}...", file=sys.stderr)
    run_command(perf_cmd)


def extract_period(header_line: str) -> tuple[float, bool]:
    match = HEADER_RE.search(header_line)
    if not match:
        return 1.0, False

    period = match.group("period")
    if period is None:
        return 1.0, False
    return float(period), True


def extract_symbol(frame_line: str) -> str | None:
    stripped = frame_line.strip()
    if not stripped:
        return None

    parts = stripped.split()
    if not parts:
        return None

    index = 0
    if HEX_RE.fullmatch(parts[0]):
        index = 1
    if index >= len(parts):
        return None

    symbol = parts[index]
    if symbol in {"(", ")", "[unknown]"}:
        return None

    symbol = symbol.split("+", 1)[0]
    symbol = symbol.split("(", 1)[0]
    return symbol or None


def classify_symbols(symbols: Sequence[str]) -> str | None:
    best_category = None
    best_depth = None

    for depth, symbol in enumerate(symbols):
        for category in CATEGORY_ORDER:
            prefixes = CATEGORY_PREFIXES[category]
            if any(symbol.startswith(prefix) for prefix in prefixes):
                if best_depth is None or depth < best_depth:
                    best_category = category
                    best_depth = depth

    return best_category


def analyze_perf_script(
    perf_exe: str,
    perf_data_path: Path,
    show_unclassified: int,
) -> AnalysisResult:
    command = [perf_exe, "script", "-i", str(perf_data_path)]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
    )

    assert process.stdout is not None
    assert process.stderr is not None

    total_weight = 0.0
    classified_weight = 0.0
    category_weights: Counter[str] = Counter()
    unclassified_leaves: Counter[str] = Counter()
    processed_blocks = 0
    weighted_by_period = False
    block: list[str] = []

    def consume_block(lines: list[str]) -> None:
        nonlocal total_weight
        nonlocal classified_weight
        nonlocal processed_blocks
        nonlocal weighted_by_period

        if not lines:
            return

        symbols = [symbol for symbol in (extract_symbol(line) for line in lines[1:]) if symbol]
        if not symbols:
            return

        weight, has_period = extract_period(lines[0])
        weighted_by_period = weighted_by_period or has_period
        total_weight += weight
        processed_blocks += 1

        category = classify_symbols(symbols)
        if category is None:
            unclassified_leaves[symbols[0]] += weight
            return

        category_weights[category] += weight
        classified_weight += weight

    for raw_line in process.stdout:
        line = raw_line.rstrip("\n")
        if line.strip():
            block.append(line)
            continue

        consume_block(block)
        block = []

    consume_block(block)

    stderr = process.stderr.read()
    return_code = process.wait()
    if return_code != 0:
        print(stderr, file=sys.stderr, end="")
        raise SystemExit(return_code)

    if stderr.strip():
        print(stderr, file=sys.stderr, end="")

    if processed_blocks == 0:
        raise SystemExit(f"No perf sample blocks were parsed from {perf_data_path}.")

    if show_unclassified < 0:
        raise SystemExit("--show-unclassified must be >= 0.")

    return AnalysisResult(
        total_weight=total_weight,
        classified_weight=classified_weight,
        category_weights=category_weights,
        processed_blocks=processed_blocks,
        weighted_by_period=weighted_by_period,
        unclassified_leaves=unclassified_leaves,
    )


def percent(numerator: float, denominator: float) -> float:
    if denominator <= 0.0:
        return 0.0
    return 100.0 * numerator / denominator


def print_report(result: AnalysisResult, show_unclassified: int) -> None:
    basis = "event periods" if result.weighted_by_period else "raw sample counts"
    print(f"Measurement basis: {basis}")
    print(f"Processed perf sample blocks: {result.processed_blocks}")
    print(f"Classified sample weight: {percent(result.classified_weight, result.total_weight):.2f}%")
    print(f"Unclassified sample weight: {percent(result.total_weight - result.classified_weight, result.total_weight):.2f}%")
    print()
    print("Phase fractions of total sampled time:")

    for phase_name, categories in PHASE_DEFINITIONS:
        phase_weight = sum(result.category_weights[category] for category in categories)
        print(f"  {phase_name}: {percent(phase_weight, result.total_weight):6.2f}% of total")
        for category in categories:
            category_weight = result.category_weights[category]
            label = category.replace(f"{phase_name}_", "", 1)
            print(
                f"    {label}: {percent(category_weight, result.total_weight):6.2f}% of total, "
                f"{percent(category_weight, phase_weight):6.2f}% of {phase_name}"
            )

    if show_unclassified == 0 or not result.unclassified_leaves:
        return

    print()
    print("Top unclassified leaf symbols:")
    for symbol, weight in result.unclassified_leaves.most_common(show_unclassified):
        print(f"  {symbol}: {percent(weight, result.total_weight):6.2f}% of total")


def main(argv: Sequence[str]) -> int:
    args, perf_command = parse_cli(argv)
    perf_data_path = Path(args.perf_data)

    if perf_command:
        record_perf(
            perf_exe=args.perf,
            output_path=perf_data_path,
            command=perf_command,
            call_graph=args.call_graph,
            mmap_pages=args.mmap_pages,
            event=args.event,
        )
    elif not perf_data_path.exists():
        raise SystemExit(
            f"{perf_data_path} does not exist. Pass an existing perf.data file, "
            "or provide a command after '--' to record one first."
        )

    result = analyze_perf_script(
        perf_exe=args.perf,
        perf_data_path=perf_data_path,
        show_unclassified=args.show_unclassified,
    )
    print_report(result, show_unclassified=args.show_unclassified)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
