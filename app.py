"""Polycheck local DNA sequence triage scanner.

This application reports heuristic triage signals only. It does not decode or
display text, encode or reconstruct DNA, or establish biological safety.
"""
from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import re
from typing import Any

import pandas as pd
import streamlit as st

from engine import (composition_score, host_specific_threshold, marker_match,
                    matched_reference_ascii_baseline,
                    scan_ascii_encoding, verdict)
from parsers import (MAX_RECORDS, MAX_TOTAL_BASES, MAX_UPLOAD_BYTES,
                     MAX_WINDOWS_PER_RECORD, ParseError, SequenceRecord,
                     parse_fasta, parse_upload)
from references import extract_accession, extract_organism_name, fetch_ncbi_reference

APP_TITLE = "Polycheck - DNA triage scanner"
MIN_STREAMLIT_VERSION = "1.49.0"
MAX_UPLOAD_MB = 10
GENERALIZED_BACKGROUND_THRESHOLD = 0.42
LIMITATIONS = [
    "This scan is a triage hint and cannot confirm a file is safe or free of hidden content.",
    "It uses composition and ASCII-triage signals only. Known-marker screening is inactive because the PRD does not provide an approved marker panel.",
    "The ASCII scan checks four base frames under raw and fixed XOR 0x3C profiles; other mappings or transforms may be missed.",
    "A candidate can arise by chance. A high score does not establish intent or prove that text is present.",
    "Host-background calibration uses 128 deterministic same-length windows for records up to 4,096 bases. It is provisional and not validated as a probability model.",
]

st.set_page_config(page_title=APP_TITLE, page_icon="🧬", layout="centered")


def _version_tuple(value: str) -> tuple[int, ...]:
    numbers = re.match(r"\d+(?:\.\d+)*", value)
    return tuple(int(part) for part in numbers.group(0).split(".")) if numbers else (0,)


def _streamlit_version_error() -> str | None:
    try:
        installed = version("streamlit")
    except PackageNotFoundError:
        installed = "not installed"
    if installed == "not installed" or _version_tuple(installed) < _version_tuple(MIN_STREAMLIT_VERSION):
        return f"Polycheck requires Streamlit {MIN_STREAMLIT_VERSION} or newer; installed version: {installed}. Update the pinned dependency before scanning."
    return None


def _near_limit_warnings(filename: str, data: bytes, records: list[SequenceRecord]) -> list[str]:
    warnings = []
    if len(data) > .9 * MAX_UPLOAD_BYTES:
        warnings.append("File is above 90% of the 10 MB size limit")
    total_bases = sum(len(record.sequence) for record in records)
    if total_bases > .9 * MAX_TOTAL_BASES:
        warnings.append("File is above 90% of the 10 M base limit")
    for record in records:
        windows = max(0, len(record.sequence) - 64 + 1)
        if windows > .9 * MAX_WINDOWS_PER_RECORD:
            warnings.append(f"Record {record.name} is above 90% of the windows-per-record limit")
    return warnings


def _parse_uploaded(name: str, data: bytes) -> list[SequenceRecord]:
    return parse_upload(name, data)


def _reference_sequences(uploaded_reference: Any, organism: str,
                         accession: str | None = None) -> tuple[list[str], str, str | None]:
    if uploaded_reference is not None:
        try:
            records = parse_fasta(uploaded_reference.getvalue())
        except ParseError as exc:
            return [], "generalized background model only", f"Uploaded reference FASTA could not be used: {exc}"
        return [record.sequence for record in records], "uploaded reference FASTA", None
    if not organism.strip() and not accession:
        return [], "generalized background model only", "No organism name found in the FASTA header, so no host-specific reference was used"
    fasta, status = fetch_ncbi_reference(organism.strip(), accession)
    if status in {"cache", "ncbi"} and fasta:
        try:
            return [record.sequence for record in parse_fasta(fasta)], "host-specific NCBI reference", None
        except ParseError:
            return [], "generalized background model only", "Host-specific reference was unavailable, so the generalized background model was used"
    if status == "not_found":
        return [], "generalized background model only", "Organism not found in NCBI, so no host-specific reference was used"
    failure_notices = {
        "unavailable_network": "NCBI could not be reached before the request timed out or the network connection failed. The generalized background model was used.",
        "unavailable_rate_limited": "NCBI rate-limited the reference request. Retry later or configure an NCBI API key; the generalized background model was used.",
        "unavailable_http": "NCBI returned an HTTP error while retrieving the reference. The generalized background model was used.",
        "unavailable_too_large": "The NCBI reference exceeded the 10 MB reference limit. The generalized background model was used.",
        "unavailable_invalid_response": "NCBI returned data that could not be validated as FASTA. The generalized background model was used.",
    }
    if status in failure_notices:
        return [], "generalized background model only", failure_notices[status]
    return [], "generalized background model only", "Host-specific reference was unavailable, so the generalized background model was used"


def _analyze(records: list[SequenceRecord], reference_sequences: list[str], reference_label: str) -> dict[str, Any]:
    results = []
    composition_threshold = host_specific_threshold(reference_sequences) if reference_sequences else GENERALIZED_BACKGROUND_THRESHOLD
    # No marker panel is bundled in v1; avoid presenting an empty panel as a negative biological finding.
    markers: tuple[str, ...] = ()
    for index, record in enumerate(records, 1):
        composition = composition_score(record.sequence)
        ascii_result = scan_ascii_encoding(record.sequence)
        background = (matched_reference_ascii_baseline(record.sequence, reference_sequences)
                      if reference_sequences else None)
        decision_ascii = float(ascii_result["evidence_score"])
        if background:
            tail = float(background["joint_tail_fraction_with_add_one"])
            if tail <= .01 and float(ascii_result["score"]) >= .78:
                # A short, unusually strong signal can be important even when it cannot
                # meet the fixed persistence rule used for long sequences.
                decision_ascii = float(ascii_result["score"])
            elif tail >= .05:
                # Do not elevate a text-like pattern that is common in matched host windows.
                decision_ascii = 0.0
        marker = marker_match(record.sequence, markers)
        score = verdict(composition, decision_ascii, marker,
                        composition_threshold)
        results.append({"record": index, "name": record.name, "composition": composition,
                        "ascii": ascii_result, "reference_background": background,
                        "decision_ascii": decision_ascii, "marker": marker, "verdict": score})
    highest = max((row["decision_ascii"] for row in results), default=0.0)
    peak_window_score = max((row["ascii"]["score"] for row in results), default=0.0)
    longest_strong_run = max((row["ascii"]["longest_strong_run"] for row in results), default=0)
    highest_composition = max((row["composition"] for row in results), default=0.0)
    final_verdict = verdict(float(highest_composition), float(highest),
                            any(row["marker"] for row in results), composition_threshold)
    all_candidates = []
    for row in results:
        all_candidates.extend({"start": c["start"] + 1, "end": c["end"], "score": c["score"], "record": row["record"]}
                              for c in row["ascii"]["candidates"])
    all_candidates.sort(key=lambda c: (-c["score"], c["start"], c["end"]))
    return {"verdict": final_verdict, "records": results, "reference": reference_label,
            "candidates": all_candidates[:100], "composition": highest_composition,
            "ascii": highest, "peak_window_score": peak_window_score,
            "longest_strong_run": longest_strong_run,
            "composition_threshold": composition_threshold,
            "reference_sequences": reference_sequences}


def _render_results(result: dict[str, Any], notice: str | None) -> None:
    colors = {"Needs manual recheck": "#C62828", "Worth a closer look": "#F9A825", "Looks OK": "#2E7D32"}
    st.markdown(f"<div style='border:2px solid {colors[result['verdict']]};border-radius:8px;padding:16px;margin:8px 0 16px'>"
                f"<h2 style='margin:0;color:{colors[result['verdict']]}'>{result['verdict']}</h2>"
                f"{'<div>triage hint only, not proof of safety</div>' if result['verdict'] == 'Looks OK' else ''}"
                f"<div style='margin-top:8px'>Reference used: {result['reference']}</div></div>", unsafe_allow_html=True)
    if notice:
        st.info(notice)
    st.info("Known-marker screening was not performed because no approved marker panel is configured.")
    st.write("The scan found signals that may deserve review. Candidate scores summarize text-like byte characteristics under the displayed hypotheses; they are not probabilities or a safety determination." if result["verdict"] != "Looks OK" else "No strong triage signal was found under the tested checks. This result does not establish safety.")
    with st.expander("Technical details", expanded=False):
        if result["candidates"]:
            for record_id in sorted({row["record"] for row in result["candidates"]}):
                st.caption(f"Record {record_id}")
                candidates = [{"start": row["start"], "end": row["end"], "score": round(row["score"], 3)}
                              for row in result["candidates"] if row["record"] == record_id]
                st.dataframe(pd.DataFrame(candidates, columns=["start", "end", "score"]), hide_index=True, use_container_width=True)
        else:
            st.caption("No candidate regions met the review threshold.")
        st.caption(f"Composition signal: {result['composition']:.3f} · ASCII triage evidence: {result['ascii']:.3f} · peak local window: {result['peak_window_score']:.3f} · adjacent strong windows: {result['longest_strong_run']} · composition review threshold: {result['composition_threshold']:.3f}. Scores are heuristic, not calibrated probabilities.")
        backgrounds = [(row["record"], row["reference_background"]) for row in result["records"]
                       if row.get("reference_background")]
        if backgrounds:
            st.caption("Matched natural-reference controls use same-length windows sampled across the selected host reference. A common signal is downweighted; an unusually strong short signal can reach manual recheck. This provisional comparison is not a calibrated probability or a safety estimate.")
            for record_id, baseline in backgrounds:
                st.caption(f"Record {record_id}: query length {baseline['query_length']:,} bases; {baseline['controls_at_or_above_query']}/{baseline['sample_count']} sampled reference windows had an equal or higher peak score; {baseline['joint_controls_at_or_above_query']}/{baseline['sample_count']} equaled or exceeded both peak score and persistence; matched-reference 95th percentile {baseline['control_q95']:.3f}; empirical tail fraction with add-one {baseline['joint_tail_fraction_with_add_one']:.3f}.")
                if baseline["joint_tail_fraction_with_add_one"] >= .05:
                    st.caption(f"Record {record_id}: the ASCII signal was downweighted because a similar or stronger peak/persistence combination appeared in at least 5% of sampled host windows.")
                elif baseline["joint_tail_fraction_with_add_one"] <= .01 and baseline["query_peak_score"] >= .78:
                    st.caption(f"Record {record_id}: the short ASCII signal was allowed to reach manual recheck because it was rare in the sampled host background.")
        elif result["reference_sequences"]:
            st.caption("Matched natural-reference comparisons are calculated for records 64–4,096 bases long; longer records need a separate, length-matched control design.")
        st.caption("An isolated high-scoring window is kept as a review clue. The strongest ASCII-based verdict requires at least 20 adjacent high-scoring windows; this persistence threshold is provisional and needs validation on a representative positive/negative set.")
        if result["records"]:
            st.caption(f"Scanned {len(result['records'])} record(s).")
    st.subheader("Limitations")
    for limitation in LIMITATIONS:
        st.markdown(f"- {limitation}")


def main() -> None:
    st.title("🧬 " + APP_TITLE)
    version_error = _streamlit_version_error()
    if version_error:
        st.error(version_error)
        return
    st.write("A local first-pass triage tool for sequence files. Files stay in memory for this session; when you request an NCBI reference, the app sends the organism name and, when present in a FASTA header, its public nucleotide accession. Sequence content is not sent.")
    st.markdown("10 MB file size  \n10 M total bases  \n10 M windows per record")
    upload = st.file_uploader("Choose a FASTA, FASTQ, or JSON file", type=["fasta", "fa", "fna", "fastq", "fq", "json"])
    if upload is None:
        st.session_state.pop("upload_bytes", None)
        st.session_state.pop("records", None)
        return
    data = upload.getvalue()
    if st.session_state.get("upload_bytes") != data:
        st.session_state.pop("scan_result", None)
        st.session_state.pop("reference_notice", None)
    try:
        if len(data) > MAX_UPLOAD_BYTES:
            raise ParseError("File exceeds the 10 MB upload limit.")
        records = _parse_uploaded(upload.name, data)
    except ParseError as exc:
        st.error(str(exc))
        return
    st.session_state["upload_bytes"] = data
    st.session_state["upload_name"] = upload.name
    st.session_state["records"] = records
    warnings = _near_limit_warnings(upload.name, data, records)
    for warning in warnings:
        st.warning(warning)
    guessed_organism = next((extract_organism_name(record.name) for record in records
                             if extract_organism_name(record.name)), "")
    accessions = [extract_accession(record.name) for record in records]
    matching_accession = (accessions[0] if accessions and accessions[0] and
                          all(accession == accessions[0] for accession in accessions) else None)
    st.subheader("Reference (optional)")
    organism = st.text_input("NCBI organism name", value=guessed_organism,
                              help="The organism name is sent to NCBI when you start the scan. If all uploaded records share an NCBI nucleotide accession in their headers, that public accession is used for a more specific reference lookup.")
    reference_upload = st.file_uploader("Or use a local reference FASTA", type=["fasta", "fa", "fna"], key="reference")
    if st.button("Run scan", type="primary", use_container_width=True):
        reference_sequences, reference_label, notice = _reference_sequences(reference_upload, organism, matching_accession)
        result = _analyze(records, reference_sequences, reference_label)
        st.session_state["scan_result"] = result
        st.session_state["reference_notice"] = notice
    result = st.session_state.get("scan_result")
    if result is not None:
        st.fragment(_render_results)(result, st.session_state.get("reference_notice"))


if __name__ == "__main__":
    main()
