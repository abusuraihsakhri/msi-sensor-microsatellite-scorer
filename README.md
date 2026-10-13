# MSI Sensor Microsatellite Scorer

Five-marker Bethesda PCR classification and mismatch repair (MMR) immunohistochemistry (IHC) interpretation, with a Python CLI, batch CSV processing, an optional FastAPI service and a static browser interface.

**Scope:** Despite the repository name, this project does **not** implement the original MSIsensor or MSIsensor2 algorithms for BAM/NGS read data. It classifies already-interpreted locus instability observations. It has not been validated as a clinical diagnostic device.

## Browser application

The root `index.html` is a no-build static application deployable through GitHub Pages. Deployment is defined in `.github/workflows/pages.yml`; check the repository's Pages deployment for current availability. The legacy `web/` entry redirects to the root application.

The browser supports individual assessment of BAT-25, BAT-26, D2S123, D5S346 and D17S250, optional MLH1/MSH2/MSH6/PMS2 expression and molecular findings, as well as CSV import and download. All calculations and CSV handling take place in the browser. No backend, account, Pyodide or external libraries are needed. Browser JavaScript implements the same classification rules as the Python core.

## Classification

The five-marker Bethesda rule used here is:

| Unstable loci | PCR classification |
| --- | --- |
| 0 of 5 | MSS |
| 1 of 5 | MSI-L |
| 2–5 of 5 | MSI-H |

MMR IHC reports retained or lost expression and flags relevant patterns for further evaluation. With both assays supplied, the tool reports PCR/IHC concordance. **Discordant results are explicitly labelled and must be reviewed**, rather than being assigned a combined positive result. MSI-L is not silently collapsed into MSS when IHC is included.

These rules do not replace assay quality control, clinically validated interpretation algorithms, genetic counselling, or treatment eligibility assessment. The illustrative prevalence lookup and narrative treatment information in the Python CLI should not be used as current guidelines.

## Command-line usage

Python 3.10–3.12 is the CI-tested range. The core scoring and CSV functions use the Python standard library.

```bash
python cli.py msi --bat25 unstable --bat26 unstable --d2s123 stable --d5s346 stable --d17s250 stable
python cli.py mmr --mlh1 lost --msh2 retained --msh6 retained --pms2 lost
python cli.py assess --bat25 unstable --bat26 unstable --d2s123 stable --d5s346 stable --d17s250 stable --mlh1 lost --msh2 retained --msh6 retained --pms2 lost
python cli.py prevalence --cancer-type colorectal
python cli.py batch -i sample.csv -o results.csv
```

CSV must contain all five marker headings exactly as shown above and explicit `stable/unstable`, `true/false` or equivalent values. If IHC data are supplied, all four MMR protein columns are required. Optional columns are `BRAF_V600E` and `MLH1_PROMOTER_METHYLATION`. Invalid rows and duplicate headers are rejected. The CLI accepts caller-selected local file paths, which are **not** sandboxed.

For the optional API:

```bash
python -m pip install "fastapi>=0.110,<1" "uvicorn>=0.30,<1" "pydantic>=2,<3"
python cli.py serve --host 127.0.0.1 --port 8000
```

API routes in `agents/api.py` provide a separate demonstration audit/supervisor workflow; **they do not expose the Bethesda scoring functions**. The configured "mock" inference provider is a deterministic placeholder and does not perform guideline verification. Running the API on a public network requires an independent authentication, access-control, logging and privacy review.

## Development and tests

```bash
python -m pip install "pytest>=8,<10" "fastapi>=0.110,<1" "pydantic>=2,<3" "uvicorn>=0.30,<1"
python -m pytest -q
node --test tests/test_browser_scoring.cjs
node --check web/scoring.js
node --check web/app.js
```

The Node tests require Node.js 22; neither npm nor a JavaScript package installation is required. CI runs Python tests on 3.10–3.12, JavaScript tests and a Docker build. To serve the static app locally, run `python -m http.server 8080` and open `http://localhost:8080/`.

A container can be built with `docker build -t msi-scorer .` and run on port 8000. Configure `AUDIT_SECRET_KEY` for stable HMAC signing across restarts. Without it, a random in-memory key is generated; the demonstration audit log is not a persistent tamper-evident ledger.

## Privacy, limitations and compatibility

- The static page does not upload data, store it in browser storage or load third-party analytics. Avoid uploading identifiable patient information even for local use.
- Batch CSV exports escape spreadsheet formula prefixes in free-text input cells.
- The optional FastAPI service is separate, server-side and **cannot run on GitHub Pages**. It is not hardened for public patient-data processing.
- Modern desktop and mobile browsers supporting `File.text()` and `Blob` are expected to work. No browser compatibility matrix has been formally validated.
- This tool is **research/educational software**, not diagnostic software. The historical Bethesda panel is not equivalent to modern comprehensive MSI NGS assays or to tissue-specific validated diagnostic algorithms.

## License

[MIT](LICENSE).
