"""Strict, in-memory sequence file parsing for Polycheck."""
from __future__ import annotations

import io
import json
import re
from dataclasses import dataclass
from typing import Any

MAX_UPLOAD_BYTES = 10 * 1024 * 1024
MAX_TOTAL_BASES = 10_000_000
MAX_WINDOWS_PER_RECORD = 10_000_000
MAX_RECORDS = 500


class ParseError(ValueError):
    """Raised when an upload is malformed or exceeds a fixed safety limit."""


@dataclass(frozen=True)
class SequenceRecord:
    name: str
    sequence: str


def _check_size(data: bytes) -> None:
    if len(data) > MAX_UPLOAD_BYTES:
        raise ParseError("File exceeds the 10 MB upload limit.")


def clean_sequence(value: str) -> str:
    """Normalize whitespace/case while retaining ambiguity as a separator."""
    normalized = re.sub(r"\s+", "", value).upper()
    if re.search(r"[^ACGTURYSWKMBDHVN.-]", normalized):
        raise ParseError("Sequence contains characters outside the DNA IUPAC alphabet.")
    return normalized


def _validate_records(records: list[SequenceRecord]) -> list[SequenceRecord]:
    if not records:
        raise ParseError("Input contains no sequence records.")
    if len(records) > MAX_RECORDS:
        raise ParseError(f"Record limit reached ({MAX_RECORDS:,}).")
    total = sum(len(record.sequence) for record in records)
    if total > MAX_TOTAL_BASES:
        raise ParseError(f"Total-base limit reached ({MAX_TOTAL_BASES:,} bases).")
    if not any(base in "ACGT" for record in records for base in record.sequence):
        raise ParseError("Input is empty or contains no valid A/C/G/T bases.")
    for record in records:
        # For 16-byte ASCII windows, each window spans 64 bases; advance by one base.
        windows = max(0, len(record.sequence) - 64 + 1)
        if windows > MAX_WINDOWS_PER_RECORD:
            raise ParseError(f"Record {record.name} exceeds the windows-per-record limit.")
    return records


def parse_fasta(data: bytes) -> list[SequenceRecord]:
    _check_size(data)
    try:
        text = data.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise ParseError("Input is not valid UTF-8 text.") from exc
    records: list[SequenceRecord] = []
    name: str | None = None
    parts: list[str] = []
    for line_no, line in enumerate(text.splitlines(), 1):
        if line.startswith(">"):
            if name is not None:
                if not parts:
                    raise ParseError(f"Invalid FASTA input: record {name!r} has no sequence (line {line_no}).")
                records.append(SequenceRecord(name, clean_sequence("".join(parts))))
            name = line[1:].strip()
            if not name:
                raise ParseError(f"Invalid FASTA input: empty header on line {line_no}.")
            parts = []
        elif line.strip():
            if name is None:
                raise ParseError(f"Invalid FASTA input: sequence before first header on line {line_no}.")
            if line.lstrip().startswith(">"):
                raise ParseError(f"Invalid FASTA input: malformed header on line {line_no}.")
            parts.append(line.strip())
    if name is not None:
        if not parts:
            raise ParseError(f"Invalid FASTA input: record {name!r} has no sequence.")
        records.append(SequenceRecord(name, clean_sequence("".join(parts))))
    return _validate_records(records)


def parse_fastq(data: bytes) -> list[SequenceRecord]:
    _check_size(data)
    try:
        text = data.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise ParseError("Input is not valid UTF-8 text.") from exc
    lines = text.splitlines()
    records: list[SequenceRecord] = []
    i = 0
    while i < len(lines):
        if not lines[i].startswith("@") or not lines[i][1:].strip():
            raise ParseError(f"Invalid FASTQ input: expected a non-empty @header on line {i + 1}.")
        name = lines[i][1:].strip()
        i += 1
        seq_parts: list[str] = []
        while i < len(lines) and not lines[i].startswith("+"):
            seq_parts.append(lines[i].strip())
            i += 1
        if i >= len(lines):
            raise ParseError(f"Invalid FASTQ input: missing + separator for record {name!r}.")
        sequence = clean_sequence("".join(seq_parts))
        if not sequence:
            raise ParseError(f"Invalid FASTQ input: record {name!r} has no sequence.")
        i += 1
        quality_length = 0
        while i < len(lines) and quality_length < len(sequence):
            if any(ord(char) < 33 or ord(char) > 126 for char in lines[i]):
                raise ParseError(f"Invalid FASTQ input: quality line contains invalid characters for record {name!r}.")
            quality_length += len(lines[i])
            i += 1
        if quality_length != len(sequence):
            raise ParseError(f"Invalid FASTQ input: sequence/quality length mismatch for record {name!r}.")
        records.append(SequenceRecord(name, sequence))
    return _validate_records(records)


def parse_json(data: bytes) -> list[SequenceRecord]:
    _check_size(data)
    try:
        obj: Any = json.loads(data.decode("utf-8", errors="strict"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ParseError(f"Invalid JSON input: {exc}") from exc
    if isinstance(obj, dict):
        obj = obj.get("records")
    if not isinstance(obj, list) or not obj:
        raise ParseError('Invalid JSON input: expected a non-empty array or an object with a "records" array.')
    records = []
    for index, item in enumerate(obj, 1):
        if not isinstance(item, dict) or not isinstance(item.get("sequence"), str):
            raise ParseError(f"Invalid JSON input: record {index} must have a string sequence field.")
        name = item.get("id", item.get("name", f"Record {index}"))
        if not isinstance(name, str) or not name.strip():
            raise ParseError(f"Invalid JSON input: record {index} has an invalid id/name.")
        records.append(SequenceRecord(name.strip(), clean_sequence(item["sequence"])))
    return _validate_records(records)


def parse_upload(name: str, data: bytes) -> list[SequenceRecord]:
    _check_size(data)
    ext = name.lower().rsplit(".", 1)[-1] if "." in name else ""
    if ext in {"fa", "fna", "fasta"}:
        return parse_fasta(data)
    if ext in {"fq", "fastq"}:
        return parse_fastq(data)
    if ext == "json":
        return parse_json(data)
    raise ParseError("Unsupported file type. Choose FASTA, FASTQ, or JSON.")
