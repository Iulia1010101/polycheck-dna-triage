# Security and scientific review notes

## Data flow

- Sequence uploads are parsed and scanned in process memory.
- Sequence content, candidate byte values, and upload identifiers are not written to logs or sent to NCBI.
- NCBI is used only when requested for a host reference; the organism name and, when available, a public nucleotide accession from the uploaded header are sent. Sequence content is not sent. The optional `NCBI_API_KEY` is read from the environment.
- Downloaded NCBI reference FASTA and metadata are cached under `~/.polycheck/references`.
- Streamlit is configured for localhost only and the upload limit is 10 MB.

## Detection boundaries

The ASCII triage checks only the documented A=00, C=01, G=10, T=11 mapping and fixed XOR profile `0x3C`, in four base frames. It evaluates byte statistics and returns only region coordinates and scores. It does not display or reconstruct text, search arbitrary keys, or claim to identify arbitrary encryption.

Composition and ASCII scores are heuristics and are not calibrated probabilities. No supplied biological evaluation set or approved marker panel is present. Known-marker screening is inactive. The app must not claim that a sequence is safe, harmful, or screened against markers.

## Review before enabling marker screening

Obtain an approved marker panel and provenance, define exact and partial-match rules, include versioning and false-positive handling, and add tests using access-approved data. Do not place sensitive sequence data in logs, screenshots, public issues, or source control. A match should be described as a review signal, not a conclusion about intent or safety.
