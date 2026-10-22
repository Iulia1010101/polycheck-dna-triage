# Polycheck PRD reference

The active product requirements are maintained at:

<https://prd.32dots.de/p/eXnv5EgM8f63k2uLCHQGiQ>

This local note replaces the earlier ASCII-in-DNA PRD, which described synthetic encoded fixtures and selectable profiles that do not match the current product requirements.

## Implemented scope

- Local Streamlit workflow bound to `127.0.0.1:8501`.
- Strict FASTA, FASTQ, and JSON parsing; 10 MB upload, 10 million total bases, and 10 million 64-base windows per record.
- Pre-scan near-limit warnings, four base frames, raw and fixed XOR `0x3C` scoring, coordinate/score-only candidates, and deduplication.
- Overlapping window hits are combined into regions. Without a matched host reference, an isolated high local score is capped to review-level evidence and the strongest ASCII verdict needs 20 adjacent strong windows.
- Optional host reference from local FASTA or NCBI, with accession-based EFetch, an offline FASTA-plus-JSON cache, specific fallback notices, and a generalized model when lookup fails.
- For records 64–4,096 bases long, up to 128 same-length natural-reference windows are compared with the query. Common host-background signals are downweighted; rare strong short signals may reach manual recheck. This provisional comparison is not calibrated evidence.
- Verdict wording, a visible “Looks OK” caveat, collapsed technical details, and limitations at the bottom.
- Parsing, engine, and Streamlit UI split across `parsers.py`, `engine.py`, and `app.py`.

## Current gap

The PRD requires known-marker matching but does not specify a marker set. The app tells researchers that known-marker screening is inactive. Do not enable it until a marker panel and matching policy are approved and provided.

The PRD refers to pre-existing tests and the earlier local implementation had a different fixture set. Treat `tests/test_app.py` as the current software test suite; passing it does not validate biological performance. The accession lookup and matched-reference calibration changes should receive additional automated coverage before they are relied on operationally.
