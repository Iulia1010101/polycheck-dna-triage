"""Local FASTA reference cache and optional NCBI organism-name lookup."""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path
import re
import threading
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from parsers import ParseError, parse_fasta

CACHE_DIR = Path.home() / ".polycheck" / "references"
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TIMEOUT_SECONDS = 12
LOGGER = logging.getLogger(__name__)
_THROTTLE_LOCK = threading.Lock()
_LAST_REQUEST_AT = 0.0


def _open_ncbi(request: Request):
    global _LAST_REQUEST_AT
    with _THROTTLE_LOCK:
        delay = 0.34 - (time.monotonic() - _LAST_REQUEST_AT)
        if delay > 0:
            time.sleep(delay)
        _LAST_REQUEST_AT = time.monotonic()
    return urlopen(request, timeout=TIMEOUT_SECONDS)


def extract_organism_name(header: str | None) -> str | None:
    if not header:
        return None
    match = re.search(r"\[([A-Za-z][A-Za-z .'-]{1,100})\]", header)
    if match:
        return match.group(1).strip()
    # NCBI flat-file FASTA headers often place the organism after the accession.
    match = re.match(
        r"^(?:[A-Z]{1,4}_?\d+(?:\.\d+)?|[A-Z][A-Z0-9_]*\.\d+)\s+"
        r"([A-Z][a-zA-Z.-]+\s+[a-z][a-zA-Z.-]+)(?:\s|,|$)",
        header.strip(),
    )
    return match.group(1).strip() if match else None


def extract_accession(header: str | None) -> str | None:
    """Extract a versioned nucleotide accession from an NCBI-style FASTA header."""
    if not header:
        return None
    match = re.search(r"\b((?:NC|NG|NM|NR|NT|NW|NZ|AC|CP|CM|AE|AL)_\d+\.\d+)\b", header.upper())
    return match.group(1) if match else None


def _read_cache(organism: str, accession: str | None = None) -> tuple[bytes, dict] | None:
    index_path = CACHE_DIR / "index.json"
    try:
        raw_index = index_path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None
    except OSError:
        LOGGER.warning("Reference cache index ignored because it could not be read.")
        return None
    try:
        index = json.loads(raw_index)
        if not isinstance(index, dict) or not isinstance(index.get("references", []), list):
            raise ValueError("invalid cache index shape")
        entries = index.get("references", [])
        entry = next((item for item in entries if isinstance(item, dict) and accession and
                      isinstance(item.get("accession"), str) and
                      item["accession"].casefold() == accession.casefold()), None)
        if entry is None and not accession:
            entry = next((item for item in entries if isinstance(item, dict) and
                          isinstance(item.get("organism_name"), str) and
                          item["organism_name"].casefold() == organism.casefold()), None)
        if not entry:
            return None
        filename = entry["fasta_filename"]
        if Path(filename).name != filename:
            return None
        content = (CACHE_DIR / filename).read_bytes()
        parse_fasta(content)
        return content, entry
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        LOGGER.warning("Reference cache entry ignored because it could not be validated.")
        return None


def _write_cache(organism: str, fasta: bytes, etag: str | None,
                 accession: str | None = None) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    identity = accession or organism
    safe = re.sub(r"[^A-Za-z0-9_-]+", "_", identity).strip("_")[:60] or "reference"
    filename = f"{safe}.fasta"
    (CACHE_DIR / filename).write_bytes(fasta)
    index_path = CACHE_DIR / "index.json"
    try:
        index = json.loads(index_path.read_text(encoding="utf-8"))
        if not isinstance(index, dict) or not isinstance(index.get("references", []), list):
            index = {"references": []}
    except (OSError, ValueError, TypeError):
        index = {"references": []}
    entries = [item for item in index["references"] if not (
        isinstance(item, dict) and ((accession and item.get("accession", "").casefold() == accession.casefold()) or
                                   (not accession and item.get("organism_name", "").casefold() == organism.casefold())))]
    entries.append({"organism_name": organism, "accession": accession, "fasta_filename": filename,
                    "downloaded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "etag": etag, "size_bytes": len(fasta)})
    temp = index_path.with_suffix(".tmp")
    temp.write_text(json.dumps({"references": entries}, indent=2), encoding="utf-8")
    temp.replace(index_path)


def fetch_ncbi_reference(organism: str, accession: str | None = None) -> tuple[bytes | None, str]:
    """Return a cached or NCBI FASTA; request only the organism label/accession."""
    if accession and not re.fullmatch(r"(?:NC|NG|NM|NR|NT|NW|NZ|AC|CP|CM|AE|AL)_\d+\.\d+", accession.upper()):
        accession = None
    accession = accession.upper() if accession else None
    cached = _read_cache(organism, accession)
    if cached:
        return cached[0], "cache"
    key = os.environ.get("NCBI_API_KEY")
    params = {"db": "nuccore", "retmode": "text"}
    if key:
        params["api_key"] = key
    try:
        if accession:
            params.update({"id": accession, "rettype": "fasta"})
        else:
            search_params = {
                "db": "nuccore",
                "term": f'"{organism}"[Organism] AND ("complete genome"[Title] OR chromosome[Title]) AND srcdb_refseq[PROP] NOT plasmid[Title]',
                "retmax": "1", "retmode": "json",
            }
            if key:
                search_params["api_key"] = key
            request = Request(f"{NCBI_BASE}/esearch.fcgi?{urlencode(search_params)}",
                              headers={"User-Agent": "Polycheck-local-triage/1.0"})
            with _open_ncbi(request) as response:
                result = json.loads(response.read().decode("utf-8"))
            ids = result.get("esearchresult", {}).get("idlist", [])
            if not ids:
                return None, "not_found"
            params.update({"id": ids[0], "rettype": "fasta"})
        request = Request(f"{NCBI_BASE}/efetch.fcgi?{urlencode(params)}", headers={"User-Agent": "Polycheck-local-triage/1.0"})
        with _open_ncbi(request) as response:
            fasta = response.read(10 * 1024 * 1024 + 1)
            etag = response.headers.get("ETag")
        if len(fasta) > 10 * 1024 * 1024:
            return None, "unavailable_too_large"
        try:
            parse_fasta(fasta)
        except ParseError:
            return None, "unavailable_invalid_response"
        try:
            _write_cache(organism, fasta, etag, accession)
        except OSError:
            LOGGER.warning("NCBI reference retrieved but could not be cached.")
        return fasta, "ncbi"
    except HTTPError as exc:
        LOGGER.warning("NCBI request failed with HTTP status %s.", exc.code)
        return None, "unavailable_rate_limited" if exc.code == 429 else "unavailable_http"
    except (TimeoutError, URLError) as exc:
        LOGGER.warning("NCBI request could not complete: %s.", type(exc).__name__)
        return None, "unavailable_network"
    except (OSError, ValueError, json.JSONDecodeError):
        LOGGER.warning("NCBI returned an unreadable or invalid reference response.")
        return None, "unavailable_invalid_response"
