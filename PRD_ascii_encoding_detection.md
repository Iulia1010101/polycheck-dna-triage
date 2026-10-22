# Polycheck PRD reference

The active product requirements are maintained at:

<https://prd.32dots.de/p/eXnv5EgM8f63k2uLCHQGiQ>

This local note replaces the earlier ASCII-in-DNA PRD, which described synthetic encoded fixtures and selectable profiles that do not match the current product requirements.

## Implemented scope

- Local Streamlit workflow bound to `127.0.0.1:8501`.
- Strict FASTA, FASTQ, and JSON parsing; 10 MB upload, 10 million total bases, and 10 million 64-base windows per record.
- Pre-scan near-limit warnings, four base frames, raw and fixed XOR `0x3C` scoring, coordinate/score-only candidates, and deduplication.
- Overlapping window hits are combined into regions; an isolated high local score is capped to review-level evidence, while the strongest ASCII verdict needs 20 adjacent strong windows. This is a provisional threshold, not calibrated evidence.
- Optional host reference from local FASTA or NCBI, with an offline FASTA-plus-JSON cache and generalized fallback.
- Verdict wording, a visible “Looks OK” caveat, collapsed technical details, and limitations at the bottom.
- Parsing, engine, and Streamlit UI split across `parsers.py`, `engine.py`, and `app.py`.

## Current gap

The PRD requires known-marker matching but does not specify a marker set. The app tells researchers that known-marker screening is inactive. Do not enable it until a marker panel and matching policy are approved and provided.

The PRD refers to 38 pre-existing tests, but this project contained 13 tests. The current suite has 21 tests; the original tests tied to the old fixture encoder were replaced to honor the current no-encode requirement.
