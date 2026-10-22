# Polycheck engineering handover

**Project:** Polycheck — local DNA sequence triage scanner  
**Handoff updated:** 22 October 2026  
**Workspace:** `DNA_Trojan_Scanner_PoC`  
**Product requirements:** [PRD](https://prd.32dots.de/p/eXnv5EgM8f63k2uLCHQGiQ)

## 1. Executive summary

Polycheck is a local Streamlit prototype for first-pass triage of FASTA, FASTQ, and JSON sequence files. It checks simple sequence-composition properties and text-like statistics under a narrow, fixed DNA-to-byte hypothesis. The UI reports candidate coordinates and a cautious review label. It does not reveal candidate bytes or decoded text, infer arbitrary keys, or establish that a sequence is safe, harmful, or intentionally encoded.

The app runs in the user's workspace at `http://127.0.0.1:8502/`. The checked-in Streamlit config defaults to port 8501; that port was already occupied, so the current instance was started on 8502. Stop only the Polycheck process when restarting it.

The most recent targeted case is `trigger.fna`, an 84-base FASTA record with header accession `NC_002516.2` (Pseudomonas aeruginosa PAO1). The app previously fell back to a generalized model and labeled it “Worth a closer look.” Changes now use that exact accession for NCBI retrieval and compare short records against same-length windows in the fetched natural host reference. In the observed scan, none of 128 sampled windows equaled or exceeded both the query's peak score and persistence; its provisional add-one tail fraction was 0.008. The current app labels it “Needs manual recheck.” This is a triage change, not a finding about intent or safety.

## 2. Start, stop, and operate

Use Python 3.10 or newer. From the project directory in PowerShell:

```powershell
python -m pip install -r requirements.txt
python -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501
```

If port 8501 is occupied, use 8502:

```powershell
python -m streamlit run app.py --server.address 127.0.0.1 --server.port 8502
```

Open the URL printed by Streamlit. Stop the process with `Ctrl+C` in the terminal where it was started. The current in-app browser uses port 8502. The checked-in `.streamlit/config.toml` binds the default server to localhost, sets the 10 MB upload limit, and disables Streamlit usage statistics.

The workflow is: upload a supported sequence file; review the organism name inferred from the FASTA header; optionally upload a local reference FASTA; select **Run scan**; inspect the verdict and expand **Technical details** for candidate coordinates and background comparisons. An NCBI lookup occurs only when the user requests a scan without supplying a local reference. The first uncached NCBI request can take several seconds.

## 3. Repository map

| File | Responsibility |
| --- | --- |
| `app.py` | Streamlit UI, upload flow, NCBI/local reference selection, per-record and file-level aggregation, verdict presentation. |
| `parsers.py` | Strict in-memory FASTA, FASTQ, and JSON parsing; IUPAC alphabet checks and input caps. |
| `engine.py` | Composition heuristic, byte-window score, candidate coordinates, persistence rule, host-background sample comparison, verdict thresholds. |
| `references.py` | Organism and accession extraction, NCBI E-utilities calls, rate spacing, FASTA cache, and reference validation. |
| `tests/test_app.py` | Parser, engine, UI requirement, and cache tests. Tests are software checks, not biological validation. |
| `README.md` | User-facing run, privacy, and interpretation notes. |
| `SECURITY_REVIEW.md` | Data-flow and detection-boundary summary. |
| `PRD_ascii_encoding_detection.md` | Local requirements/implementation note. Some historical statements in it may be stale; this handover describes current code. |
| `.streamlit/config.toml` | Local server defaults and upload limit. |

Dependencies are intentionally small: Streamlit, pandas, and NumPy. There is no database service or separate worker process.

## 4. Input and resource limits

- Upload limit: 10 MiB.
- Total sequence bases per upload: 10,000,000.
- Maximum parsed records: 500.
- Maximum 64-base candidate windows per record: 10,000,000.
- The UI warns when file size, total bases, or windows exceed 90% of their cap.
- FASTA is parsed as UTF-8 and requires a header before sequence. Sequence whitespace is removed and bases are uppercased. The accepted alphabet is DNA IUPAC; ambiguous symbols break ASCII scanning frames.
- FASTQ quality text is consumed for length validation and is never scored as a sequence.
- JSON accepts an array of records or an object containing a `records` array. Each record needs a string `sequence`; `id` or `name` is optional.
- Upload sequence data is held in process/session memory. It is not intentionally written to disk or included in logs.

## 5. Detection behavior

### 5.1 ASCII-like triage

`scan_ascii_encoding()` in `engine.py` maps A/C/G/T to two-bit values A=00, C=01, G=10, T=11. It packs groups of four bases into one byte, tests four base-frame offsets, and scores overlapping 16-byte windows. This corresponds to 64 bases per complete byte window. It evaluates two hypotheses only: raw bytes and fixed XOR `0x3C`.

The window score is a heuristic weighted from printable-byte fraction (0.65), alphabetic fraction (0.25), and presence of selected character classes (0.10). It is not a probability and does not verify that a meaningful message exists. The scanner returns region start/end coordinates and scores only; it does not return bytes or decoded content. The app displays coordinates as 1-based inclusive ranges.

Windows at or above 0.58 form candidate/review regions. Scores at or above 0.78 are considered “strong” for persistence tracking. When the longest adjacent strong run is shorter than 20 windows, the aggregate ASCII evidence is capped at 0.77. This guard was introduced because an isolated high-scoring window can occur by chance when scanning many correlated windows. The cutoff is provisional and has not been validated on a representative benchmark.

### 5.2 Composition heuristic

`composition_score()` combines base entropy, GC deviation from 50%, and the longest homopolymer relative to sequence length. It returns a heuristic from 0 to 1, not a pathogenicity or safety measure. The host-specific composition threshold is estimated from non-overlapping 4,096-base pieces of the selected reference: its 99th percentile is clamped to 0.25–0.90. Without a usable reference, the generalized threshold is 0.42.

### 5.3 Matched natural-reference background

For a record between 64 and 4,096 bases and an available reference, `matched_reference_ascii_baseline()` samples up to 128 same-length windows from contiguous A/C/G/T runs in the reference. Sampling positions are evenly spaced across all eligible start positions. It compares the query and controls by peak local score and longest strong run, and reports a 95th percentile and an add-one empirical tail fraction. The reference sequences themselves are not shown in the UI.

This background comparison currently changes ASCII evidence as follows:

- If the add-one joint tail fraction is at most 0.01 and query peak score is at least 0.78, the peak score can reach the manual-recheck threshold even when the record is too short to meet the 20-window persistence rule.
- If at least 5% of sampled controls equal or exceed both peak score and persistence, the ASCII contribution is set to zero for verdict calculation; the candidate and raw scores remain visible.
- In between those conditions, the original persistence-adjusted evidence score is used.
- Records shorter than 64 or longer than 4,096 bases have no matched ASCII background; they use the original persistence behavior.

The add-one fraction is a descriptive sample rank, **not a p-value, calibrated false-positive rate, or safety estimate**. With 128 controls its resolution is coarse, and one reference genome cannot represent within-species variation. The 1%/5% decision cutoffs are provisional and should not be treated as operational thresholds until evaluated on held-out positive and benign data. A naturally occurring reference genome is a genomic background control; “natural” does not mean harmless or biologically safe.

### 5.4 Verdicts and markers

The verdict function reports:

- **Needs manual recheck** for an enabled marker match, adjusted ASCII evidence at least 0.78, or composition above `max(0.72, composition threshold)`.
- **Worth a closer look** for adjusted ASCII evidence at least 0.58 or composition at/above the composition threshold.
- **Looks OK** otherwise, with explicit UI language that this does not establish safety.

No marker panel is configured. Marker matching is inactive; the UI says so. Do not interpret the absence of a marker match as a negative result or add a marker list without an approved source, version, matching policy, and review of false-positive behavior.

## 6. NCBI lookup and data handling

`references.py` uses NCBI E-utilities at `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`:

1. If every uploaded record header contains the same supported versioned nucleotide accession, the app extracts it and calls EFetch directly. The current accepted prefixes include NC, NG, NM, NR, NT, NW, NZ, AC, CP, CM, AE, and AL. Accession.version identifiers are supported by NCBI EFetch for sequence records ([NCBI E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/)). Assembly accessions such as GCA/GCF are not currently resolved directly; the FASTA header may still provide a chromosome accession such as NC_… .
2. Otherwise, if an organism name is available, the app searches `nuccore` for a RefSeq chromosome or complete-genome title while excluding plasmid titles, then fetches one result.
3. Downloads are capped at 10 MiB and must parse as FASTA. Requests have a 12-second timeout and are spaced by at least 0.34 seconds. An optional `NCBI_API_KEY` environment variable is sent only to NCBI.
4. Valid references are cached under `~/.polycheck/references` as FASTA plus `index.json`. Cache lookup is exact by accession when one is provided; organism-only lookups use a case-insensitive organism label.
5. A cache write failure is logged as a warning but does not discard a successfully fetched reference.

The app sends the organism label and, when selected, the public accession. It does not send uploaded sequence content. Uploaded reference FASTA stays local. Downloaded public reference sequence is written to the local cache. The app reports distinct user-facing fallback messages for network failure, rate limiting, HTTP failure, invalid response, excessive reference size, and no search result. Log messages avoid sequence content and query values.

The prior UI hid all NCBI exceptions behind one generic “unavailable” notice, which made the exact cause unrecoverable from that scan. On the latest trigger scan, direct retrieval of NC_002516.2 succeeded and returned one FASTA record (6,353,955 bytes). This confirms the direct accession path works from the current environment; it does not establish that all network environments can reach NCBI.

## 7. Recent trigger investigation

The inspected `trigger.fna` contains one 84-base record. Before the host-reference comparison, its peak local ASCII window score was 0.920, but only six adjacent strong windows were present. The persistence guard therefore capped its aggregate ASCII evidence at 0.770 and returned “Worth a closer look.”

After fetching the exact chromosome reference, the same-length controls had a 95th-percentile peak score of 0.693. None of 128 sampled controls equaled or exceeded both the query peak (0.920) and its persistence (six strong windows). The add-one sample fraction was 1/129 ≈ 0.00775. The host-aware short-record rule raised adjusted ASCII evidence to 0.920, and the UI displayed “Needs manual recheck.” It did not decode or display the sequence content.

This example is one observed case, not an estimate of sensitivity or precision. The accession identifies a natural P. aeruginosa PAO1 reference; the comparison is useful for sequence background, not a basis for calling either the input or organism benign.

## 8. Known limitations and technical debt

1. **No benchmark.** There is no representative, access-approved set of benign and positive cases, no held-out validation, and no calibrated operating point. The persistence and empirical-tail thresholds are manually selected.
2. **Small control sample.** Up to 128 deterministic windows from one reference can miss rare genomic contexts and does not model species/strain diversity. Adjacent genomic samples are correlated; evenly spaced positions do not make controls independent.
3. **Length coverage.** The matched null is only implemented for records 64–4,096 bases. Whole chromosomes and records above that limit fall back to the 20-window persistence rule; they do not receive a length-matched maximum-score null.
4. **Narrow encoding hypotheses.** Only four frame offsets, the documented 2-bit base mapping, raw output, and fixed XOR `0x3C` are evaluated. Other encodings and transforms may be missed. The app must not add arbitrary key search or render decoded payloads without a separate scope and safety review.
5. **Reference selection.** Organism-only lookup picks one ESearch result. It is better constrained than the former unconstrained `retmax=1` search, but it is not a curated pangenome or guaranteed canonical representative. The app currently uses one accession/reference per scan.
6. **Reference cache policy.** Cache entries have no expiry or explicit NCBI refresh path. They are parsed before use but cache provenance/version updates need a deliberate policy.
7. **Marker screening absent.** The PRD mentions markers but supplies no approved panel. This feature remains off.
8. **Local-only deployment.** There is no authentication, multi-user isolation, or hardened network deployment. Keep it bound to localhost unless deployment architecture and data controls are reviewed.
9. **Potential UI/API cleanup.** Streamlit logs currently warn that `use_container_width` will be removed after 2025-12-31; migrate to the supported `width` argument when convenient.
10. **Historical document drift.** `PRD_ascii_encoding_detection.md` and earlier handoff notes contain older test counts and implementation descriptions. Prefer the active linked PRD, source code, and this handover; reconcile archival notes during the next documentation pass.

## 9. Recommended next work

### First: validate software behavior

- Run the automated suite from the project directory: `python -m unittest discover -s tests`.
- Add coverage for accession extraction and the “all records share one accession” rule; mixed or missing accessions must fall back to organism lookup.
- Add mocked E-utilities cases for accession EFetch, organism ESearch/EFetch, HTTP 429, timeouts, malformed JSON/FASTA, oversized response, and cache-write failure. Tests must not require live NCBI access.
- Add engine cases for a host-common signal, a host-rare short strong signal, the middle tail band, no usable reference, and records outside 64–4,096 bases.
- Check that failed reference retrieval produces the correct notice and does not silently imply host-specific analysis.

### Then: establish biological performance

- Build a documented, permissioned validation set spanning multiple taxa, strains, sequence lengths, GC content, ambiguous bases, and natural genomic regions. Keep positive controls separately managed and access-approved.
- Split data by genome/strain for train, calibration, and held-out evaluation. Do not split adjacent windows randomly across sets; that leaks highly correlated sequence context.
- Measure record-level and window-level false-positive rates, sensitivity, precision/PPV, and review burden with confidence intervals. Stratify short fragments and whole genomes because their multiple-testing behavior differs substantially.
- Compare the current deterministic 128-window scheme with larger, block-aware reference sampling and an explicit length/multiple-testing null. Calibrate the 1%/5% and 20-window thresholds against the held-out data before describing rates or changing defaults.
- Expand references only from documented NCBI/approved sources and preserve accession/version, retrieval time, and checksum as local provenance. A curated, versioned multi-genome panel is more defensible than treating a single chromosome as representative of a species.
- Keep results as review signals. Do not make automatic safety, benignity, or intent decisions from these heuristics.

## 10. Validation record for this handoff

- The provided trigger was parsed and analyzed without emitting sequence content.
- Direct NCBI retrieval by NC_002516.2 succeeded from the current environment; the reference was 6,353,955 bytes and parsed as one FASTA record.
- The trigger was rerun in the live Streamlit UI after the code update. The UI showed “Needs manual recheck,” “host-specific NCBI reference,” and the 0/128 matched-control comparison.
- The automated suite was **not rerun after the latest accession/background-calibration changes**. Run it before merging or relying on the changes beyond this observed case.

