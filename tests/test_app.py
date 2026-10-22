import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import app
import engine
import parsers
import references
from parsers import ParseError, SequenceRecord


class ParserTests(unittest.TestCase):
    def test_fasta_multiline_and_header_preserved(self):
        records = parsers.parse_fasta(b">sample [Bacterium example]\nACGT\nNNAC\n")
        self.assertEqual(records, [SequenceRecord("sample [Bacterium example]", "ACGTNNAC")])

    def test_fasta_rejects_sequence_before_header(self):
        with self.assertRaisesRegex(ParseError, "before first header"):
            parsers.parse_fasta(b"ACGT\n>sample\nACGT\n")

    def test_fasta_rejects_empty_record(self):
        with self.assertRaisesRegex(ParseError, "has no sequence"):
            parsers.parse_fasta(b">empty\n>next\nACGT\n")

    def test_invalid_alphabet_is_rejected(self):
        with self.assertRaisesRegex(ParseError, "IUPAC alphabet"):
            parsers.parse_fasta(b">sample\nACGTxyz\n")

    def test_fastq_quality_is_not_in_sequence(self):
        records = parsers.parse_fastq(b"@read\nACGTAC\n+\nACGTAC\n")
        self.assertEqual(records[0].sequence, "ACGTAC")

    def test_fastq_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ParseError, "length mismatch"):
            parsers.parse_fastq(b"@read\nACGTAC\n+\nIIII\n")

    def test_fastq_quality_line_with_header_chars_stays_quality(self):
        records = parsers.parse_fastq(b"@read\nACGT\n+\n@+AC\n")
        self.assertEqual(records[0].sequence, "ACGT")

    def test_json_records_array(self):
        data = json.dumps({"records": [{"id": "r1", "sequence": "acgt"}]}).encode()
        self.assertEqual(parsers.parse_json(data), [SequenceRecord("r1", "ACGT")])

    def test_json_malformed_rejected(self):
        with self.assertRaisesRegex(ParseError, "Invalid JSON"):
            parsers.parse_json(b"{")

    def test_upload_size_checked_before_parser(self):
        with patch.object(parsers, "MAX_UPLOAD_BYTES", 3):
            with self.assertRaisesRegex(ParseError, "10 MB upload limit"):
                parsers.parse_upload("x.fasta", b">x\nACGT\n")

    def test_unsupported_extension_rejected(self):
        with self.assertRaisesRegex(ParseError, "Unsupported file type"):
            parsers.parse_upload("x.txt", b"ACGT")

    def test_record_and_total_base_caps(self):
        with patch.object(parsers, "MAX_TOTAL_BASES", 3):
            with self.assertRaisesRegex(ParseError, "Total-base limit"):
                parsers.parse_fasta(b">x\nACGT\n")
        with patch.object(parsers, "MAX_RECORDS", 1):
            with self.assertRaisesRegex(ParseError, "Record limit"):
                parsers.parse_fasta(b">x\nACGT\n>y\nACGT\n")

    def test_window_cap(self):
        with patch.object(parsers, "MAX_WINDOWS_PER_RECORD", 1):
            with self.assertRaisesRegex(ParseError, "windows-per-record"):
                parsers.parse_fasta(b">x\n" + b"A" * 66 + b"\n")


class EngineTests(unittest.TestCase):
    def test_scanner_returns_candidate_coordinates_and_score_only(self):
        result = engine.scan_ascii_encoding("A" * 256)
        self.assertTrue(result["candidates"])
        candidate = result["candidates"][0]
        self.assertEqual(set(candidate), {"start", "end", "score"})
        self.assertGreaterEqual(candidate["score"], .58)
        self.assertLessEqual(candidate["end"], 256)
        self.assertNotIn("text", result)
        self.assertNotIn("bytes", candidate)

    def test_ambiguous_base_splits_candidates(self):
        result = engine.scan_ascii_encoding("A" * 80 + "N" + "A" * 80)
        self.assertTrue(result["candidates"])
        self.assertTrue(all(row["end"] <= 80 or row["start"] >= 81 for row in result["candidates"]))

    def test_candidate_sort_and_deduplication(self):
        result = app._analyze([SequenceRecord("r", "A" * 256)], [], "generalized background model only")
        candidates = result["candidates"]
        self.assertEqual(candidates, sorted(candidates, key=lambda row: (-row["score"], row["start"], row["end"])))
        coords = [(row["start"], row["end"]) for row in candidates]
        self.assertEqual(len(coords), len(set(coords)))

    def test_overlapping_window_hits_are_merged_into_regions(self):
        result = engine.scan_ascii_encoding("A" * 256)
        self.assertEqual(len(result["candidates"]), 1)
        self.assertEqual(result["candidates"][0]["start"], 0)
        self.assertEqual(result["candidates"][0]["end"], 256)

    def test_verdict_labels(self):
        self.assertEqual(engine.verdict(.1, .1, False), "Looks OK")
        self.assertEqual(engine.verdict(.5, .1, False), "Worth a closer look")
        self.assertEqual(engine.verdict(.1, .8, False), "Needs manual recheck")

    def test_isolated_high_window_cannot_create_strong_verdict(self):
        evidence = engine.aggregate_evidence_score(.92, 11)
        self.assertEqual(evidence, .77)
        self.assertEqual(engine.verdict(.15, evidence, False), "Worth a closer look")
        self.assertEqual(engine.aggregate_evidence_score(.92, engine.MIN_PERSISTENT_STRONG_WINDOWS), .92)

    def test_host_threshold_is_reference_derived(self):
        self.assertEqual(engine.host_specific_threshold([]), .42)
        self.assertGreaterEqual(engine.host_specific_threshold(["ACGT" * 2048]), .25)

    def test_marker_match_requires_explicit_marker_panel(self):
        self.assertFalse(engine.marker_match("ACGT" * 20))
        self.assertTrue(engine.marker_match("AACCGGTT", ["CCGG"]))


class UIRequirementTests(unittest.TestCase):
    def test_start_page_contains_three_fixed_caps(self):
        import inspect
        source = inspect.getsource(app.main)
        self.assertIn("10 MB file size", source)
        self.assertIn("10 M total bases", source)
        self.assertIn("10 M windows per record", source)

    def test_near_limit_warnings_name_each_limit(self):
        records = [SequenceRecord("r1", "A" * 9_000_100)]
        with patch.object(app, "MAX_UPLOAD_BYTES", 10_000_000), patch.object(app, "MAX_TOTAL_BASES", 10_000_000), patch.object(app, "MAX_WINDOWS_PER_RECORD", 10_000_000):
            warnings = app._near_limit_warnings("large.fa", b"x" * 9_000_100, records)
        self.assertIn("File is above 90% of the 10 MB size limit", warnings)
        self.assertIn("File is above 90% of the 10 M base limit", warnings)
        self.assertIn("Record r1 is above 90% of the windows-per-record limit", warnings)

    def test_looks_ok_disclaimer_is_rendered(self):
        import inspect
        self.assertIn("triage hint only, not proof of safety", inspect.getsource(app._render_results))

    def test_no_organism_name_has_visible_fallback_notice(self):
        records, reference, notice = app._reference_sequences(None, "")
        self.assertEqual(records, [])
        self.assertEqual(reference, "generalized background model only")
        self.assertIn("No organism name found", notice)


class ReferenceCacheTests(unittest.TestCase):
    def test_organism_extracted_from_ncbi_accession_header(self):
        self.assertEqual(references.extract_organism_name(
            "AE004091.2 Pseudomonas aeruginosa PAO1, complete genome"),
            "Pseudomonas aeruginosa")

    def test_reference_cache_round_trip_avoids_network(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(references, "CACHE_DIR", Path(directory)):
            fasta = b">reference\nACGTACGT\n"
            references._write_cache("Test organism", fasta, "etag")
            loaded = references._read_cache("test organism")
        self.assertIsNotNone(loaded)
        self.assertEqual(loaded[0], fasta)

    def test_corrupt_cache_is_ignored(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(references, "CACHE_DIR", Path(directory)):
            path = Path(directory) / "index.json"
            path.write_text("{", encoding="utf-8")
            with self.assertLogs(references.LOGGER, level="WARNING"):
                self.assertIsNone(references._read_cache("Test organism"))


if __name__ == "__main__":
    unittest.main()
