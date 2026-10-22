"""Bounded DNA triage signals. No sequence or decoded text is returned."""
from __future__ import annotations

import re
from typing import Iterable

import numpy as np

XOR_KEY = 0x3C
ASCII_WINDOW_BYTES = 16
ASCII_WINDOW_BASES = ASCII_WINDOW_BYTES * 4
MAX_CANDIDATES = 100
MIN_PERSISTENT_STRONG_WINDOWS = 20


def _candidate_scores(values: np.ndarray) -> np.ndarray:
    printable = ((values >= 32) & (values <= 126)).mean(axis=1)
    letters = (((values >= 65) & (values <= 90)) | ((values >= 97) & (values <= 122))).mean(axis=1)
    digit = ((values >= 48) & (values <= 57)).any(axis=1)
    punctuation = np.isin(values, np.frombuffer(b" _-.,:;!?/", dtype=np.uint8)).any(axis=1)
    classes = (np.ones(len(values), dtype=np.uint8) + digit + punctuation) / 3.0
    return 0.65 * printable + 0.25 * letters + 0.10 * classes


def _scan_profiles(sequence: str) -> tuple[list[dict[str, float | int]], list[dict[str, float | int]], dict[str, dict[str, float | int]]]:
    # Encode each base only as its 2-bit value for scoring; never convert bytes to text.
    codes = np.fromiter(("ACGT".find(base) for base in sequence), dtype=np.int8, count=len(sequence))
    valid = codes >= 0
    raw_output: list[dict[str, float | int]] = []
    xor_output: list[dict[str, float | int]] = []
    profile_stats = {
        "raw": {"max_score": 0.0, "review_windows": 0, "strong_windows": 0, "longest_strong_run": 0},
        "xor": {"max_score": 0.0, "review_windows": 0, "strong_windows": 0, "longest_strong_run": 0},
    }
    # A non-ACGT character breaks each frame so candidates cannot cross it.
    edges = np.flatnonzero(np.diff(np.r_[False, valid, False].astype(np.int8)))
    runs = zip(edges[::2], edges[1::2])
    for run_start, run_end in runs:
        run = codes[run_start:run_end].astype(np.uint8)
        if len(run) < ASCII_WINDOW_BASES:
            continue
        for frame in range(4):
            framed = run[frame:]
            framed = framed[: (len(framed) // 4) * 4]
            if len(framed) < ASCII_WINDOW_BASES:
                continue
            byte_stream = ((framed[0::4] << 6) | (framed[1::4] << 4) |
                           (framed[2::4] << 2) | framed[3::4]).astype(np.uint8)
            count = len(byte_stream) - ASCII_WINDOW_BYTES + 1
            streaks = {"raw": 0, "xor": 0}
            # Process small batches; stride views avoid allocating a whole window matrix.
            for offset in range(0, count, 250_000):
                stop = min(count, offset + 250_000)
                view = np.lib.stride_tricks.sliding_window_view(
                    byte_stream[offset:stop + ASCII_WINDOW_BYTES - 1], ASCII_WINDOW_BYTES)
                for windows, output, profile in ((view, raw_output, "raw"),
                                                 (np.bitwise_xor(view, XOR_KEY), xor_output, "xor")):
                    scores = _candidate_scores(windows)
                    stats = profile_stats[profile]
                    stats["max_score"] = max(float(stats["max_score"]), float(scores.max(initial=0.0)))
                    stats["review_windows"] = int(stats["review_windows"]) + int(np.count_nonzero(scores >= .58))
                    strong = scores >= .78
                    stats["strong_windows"] = int(stats["strong_windows"]) + int(np.count_nonzero(strong))
                    # Track contiguous high-scoring byte windows, including across batches.
                    edges = np.flatnonzero(np.diff(np.r_[False, strong, False].astype(np.int8)))
                    chunk_run = int(np.max(edges[1::2] - edges[::2], initial=0))
                    all_strong = bool(np.all(strong))
                    prefix = (len(strong) if all_strong else int(np.argmax(~strong))) if strong[0] else 0
                    suffix = (len(strong) if all_strong else int(np.argmax(~strong[::-1]))) if strong[-1] else 0
                    current_run = streaks[profile] + prefix if strong[0] else 0
                    longest = max(chunk_run, current_run, suffix)
                    stats["longest_strong_run"] = max(int(stats["longest_strong_run"]), longest)
                    streaks[profile] = (streaks[profile] + len(strong) if all_strong else suffix) if strong[-1] else 0
                    candidate_mask = scores >= .58
                    candidate_edges = np.flatnonzero(np.diff(np.r_[False, candidate_mask, False].astype(np.int8)))
                    batch_candidates = []
                    for first, after_last in zip(candidate_edges[::2], candidate_edges[1::2]):
                        start = int(run_start + frame + (offset + first) * 4)
                        end = start + ASCII_WINDOW_BASES + (int(after_last) - int(first) - 1) * 4
                        batch_candidates.append({"start": start, "end": end,
                                                 "score": float(scores[first:after_last].max())})
                    batch_candidates.sort(key=lambda row: (-float(row["score"]), int(row["start"]), int(row["end"])))
                    output.extend(batch_candidates[:MAX_CANDIDATES])
    for output in (raw_output, xor_output):
        output.sort(key=lambda row: (-float(row["score"]), int(row["start"]), int(row["end"])))
    return raw_output[:MAX_CANDIDATES], xor_output[:MAX_CANDIDATES], profile_stats


def scan_ascii_encoding(sequence: str) -> dict[str, object]:
    """Evaluate raw and fixed XOR 0x3C profiles and deduplicate coordinates."""
    raw, masked, profile_stats = _scan_profiles(sequence)
    best_by_region: dict[tuple[int, int], float] = {}
    for row in raw + masked:
        key = (int(row["start"]), int(row["end"]))
        best_by_region[key] = max(best_by_region.get(key, 0.0), float(row["score"]))
    window_candidates = [{"start": start, "end": end, "score": score}
                         for (start, end), score in best_by_region.items()]
    window_candidates.sort(key=lambda row: (row["start"], row["end"]))
    candidates: list[dict[str, float | int]] = []
    for row in window_candidates:
        if candidates and row["start"] < candidates[-1]["end"]:
            candidates[-1]["end"] = max(candidates[-1]["end"], row["end"])
            candidates[-1]["score"] = max(candidates[-1]["score"], row["score"])
        else:
            candidates.append(dict(row))
    candidates.sort(key=lambda row: (-row["score"], row["start"], row["end"]))
    candidates = candidates[:MAX_CANDIDATES]
    best = candidates[0]["score"] if candidates else 0.0
    longest_strong_run = max(int(stats["longest_strong_run"]) for stats in profile_stats.values())
    # A single high window is expected when millions of correlated windows are tested.
    # Keep it visible as a candidate, but reserve the strongest triage level for sustained evidence.
    evidence_score = aggregate_evidence_score(best, longest_strong_run)
    return {"score": best, "evidence_score": evidence_score,
            "longest_strong_run": longest_strong_run,
            "strong_windows": sum(int(stats["strong_windows"]) for stats in profile_stats.values()),
            "review_windows": sum(int(stats["review_windows"]) for stats in profile_stats.values()),
            "candidates": candidates,
            "profile": "Raw and fixed XOR 0x3C", "windows_examined": max(0, len(sequence) - ASCII_WINDOW_BASES + 1)}


def aggregate_evidence_score(best_window_score: float, longest_strong_run: int) -> float:
    """Down-weight isolated high windows; persistence supports a stronger triage level."""
    return float(best_window_score if longest_strong_run >= MIN_PERSISTENT_STRONG_WINDOWS
                 else min(best_window_score, .77))


def composition_score(sequence: str) -> float:
    """Heuristic composition irregularity score on a 0–1 scale, not a probability."""
    runs = [run for run in re.split(r"[^ACGT]+", sequence.upper()) if run]
    valid = "".join(runs)
    if len(valid) < 64:
        return 0.0
    counts = np.array([valid.count(base) for base in "ACGT"], dtype=float)
    p = counts / counts.sum()
    entropy = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
    gc = float((counts[1] + counts[2]) / counts.sum())
    longest = max((len(m.group(0)) for run in runs for m in re.finditer(r"A+|C+|G+|T+", run)), default=1)
    gc_extreme = min(1.0, abs(gc - .5) / .5)
    entropy_low = min(1.0, max(0.0, (2.0 - entropy) / 2.0))
    repeat = min(1.0, longest / max(1, len(valid)) * 12)
    return float(.45 * entropy_low + .40 * gc_extreme + .15 * repeat)


def host_specific_threshold(reference_sequences: Iterable[str]) -> float:
    """Estimate a conservative composition review threshold from local reference windows."""
    scores = []
    for sequence in reference_sequences:
        for start in range(0, len(sequence), 4096):
            chunk = sequence[start:start + 4096]
            if len(chunk) >= 64:
                scores.append(composition_score(chunk))
    if not scores:
        return 0.42
    # A high reference percentile avoids flagging composition that is ordinary for this host.
    return max(0.25, min(0.90, float(np.quantile(scores, 0.99))))


def matched_reference_ascii_baseline(query: str, reference_sequences: Iterable[str],
                                     sample_count: int = 128) -> dict[str, float | int] | None:
    """Compare a query's peak score with deterministic same-length host reference windows.

    This is a descriptive matched-background comparison, not a calibrated probability or
    a verdict threshold. Sampling is evenly spaced over valid A/C/G/T runs.
    """
    length = len(query)
    if length < ASCII_WINDOW_BASES or length > 4096 or sample_count < 1:
        return None
    runs: list[tuple[str, int, int]] = []
    total_starts = 0
    for sequence in reference_sequences:
        normalized = sequence.upper()
        for match in re.finditer(r"[ACGT]+", normalized):
            start, end = match.span()
            starts = end - start - length + 1
            if starts > 0:
                runs.append((normalized, start, starts))
                total_starts += starts
    if not total_starts:
        return None
    sample_n = min(sample_count, total_starts)
    targets = np.linspace(0, total_starts - 1, num=sample_n, dtype=np.int64)
    offsets = np.cumsum([0, *(starts for _, _, starts in runs)])
    query_scan = scan_ascii_encoding(query)
    query_score = float(query_scan["score"])
    query_run = int(query_scan["longest_strong_run"])
    controls = []
    run_index = 0
    for target in targets:
        while run_index + 1 < len(offsets) - 1 and target >= offsets[run_index + 1]:
            run_index += 1
        sequence, run_start, _ = runs[run_index]
        start = run_start + int(target - offsets[run_index])
        control_scan = scan_ascii_encoding(sequence[start:start + length])
        controls.append((float(control_scan["score"]), int(control_scan["longest_strong_run"])))
    joint_exceedances = sum(score >= query_score and run >= query_run for score, run in controls)
    return {
        "sample_count": sample_n,
        "query_length": length,
        "query_peak_score": query_score,
        "control_q95": float(np.quantile([score for score, _ in controls], .95)),
        "controls_at_or_above_query": int(sum(score >= query_score for score, _ in controls)),
        "joint_controls_at_or_above_query": int(joint_exceedances),
        "joint_tail_fraction_with_add_one": float((joint_exceedances + 1) / (sample_n + 1)),
    }


def marker_match(sequence: str, markers: Iterable[str] = ()) -> bool:
    """Exact local match against explicitly supplied markers; empty set means unavailable."""
    normalized = sequence.upper()
    return any(marker and marker.upper() in normalized for marker in markers)


def verdict(composition: float, ascii_score: float, marker: bool,
            composition_threshold: float = .42) -> str:
    if marker or ascii_score >= .78 or composition >= max(.72, composition_threshold):
        return "Needs manual recheck"
    if ascii_score >= .58 or composition >= composition_threshold:
        return "Worth a closer look"
    return "Looks OK"
