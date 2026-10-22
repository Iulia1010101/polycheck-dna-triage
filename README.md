# Polycheck - DNA triage scanner

Polycheck is a local, first-pass triage tool for FASTA, FASTQ, and JSON sequence files. It reports a cautious triage hint and candidate coordinates for regions with text-like byte statistics under two fixed hypotheses. It never displays decoded text, infers XOR keys, or claims that a sequence is safe or harmful. Without a matched host reference, the strongest ASCII-based verdict requires 20 adjacent high-scoring windows. With a reference, common same-length host-background signals are downweighted, while a rare strong signal in a short record can reach manual recheck. These are provisional heuristics, not validated statistical thresholds.

## Run locally

Use Python 3.10 or newer:

```powershell
python -m pip install -r requirements.txt
python -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501
```

The Streamlit configuration also binds to `127.0.0.1:8501` and caps uploads at 10 MB. The scanner accepts at most 10 million total bases and 10 million 64-base windows per record. Sequence uploads remain in memory for the session.

## References and privacy

The optional host reference can be an uploaded local FASTA or an NCBI reference. When NCBI is selected, the app sends the organism name and, when all uploaded headers share one versioned nucleotide accession, that public accession; sequence content is never sent. Accession lookups use EFetch directly, while organism-only lookups search for a RefSeq chromosome/complete genome and exclude plasmids. References are cached as FASTA plus a JSON index in `~/.polycheck/references`. An optional NCBI API key may be provided through the `NCBI_API_KEY` environment variable and is not stored by the app.

When a reference is available, records from 64 to 4,096 bases are compared with up to 128 deterministic, same-length windows sampled from that natural reference. Common text-like signals are downweighted; a very unusual short signal can reach manual recheck even when it is too short to satisfy the persistence rule. This is an exploratory control, not a calibrated probability model. Its thresholds must be validated against representative benign and positive datasets before operational use.

The cache stores reference sequences fetched from NCBI. Input sequences are not written to disk or logs. Keep the app local and do not place sequence files, identifiers, or screenshots in public issues or support requests.

## Interpretation limits

`Looks OK` is a triage hint only, not proof of safety. The scores are heuristic, not probabilities. Natural sequence composition can produce false alarms, and different encodings or transforms can be missed. A candidate is a coordinate range plus score; the application does not decode or display the underlying byte stream.

The PRD requires a known-marker match but does not provide a marker panel or approved marker source. No panel is bundled in this version, so marker matching is currently inactive and must not be interpreted as a negative finding. Supply an approved local marker definition before relying on that signal.

## Tests

```powershell
python -m unittest discover -s tests
```

The test suite covers parser strictness, FASTQ quality-line separation, resource caps, warnings, candidate coordinates, scoring, and cautious UI wording. These software tests do not validate biological performance.
